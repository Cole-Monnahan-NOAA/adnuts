#' Bayesian inference of a TMB model using the no-U-turn sampler.
#'
#' Draw Bayesian posterior samples from a Template Model Builder (TMB)
#' model using an MCMC algorithm. This function generates posterior samples
#' from which inference can be made. Adaptation schemes are used so
#' specification tuning parameters are not necessary, and parallel
#' execution reduces overall run time.
#'
#' @details This function implements algorithm 6 of Hoffman and Gelman (2014),
#' and loosely follows package \code{rstan}. The step size can be
#'   adapted or specified manually. The metric (i.e., mass matrix) can be
#'   unit diagonal, adapted diagonal (default and recommended), or a dense
#'   matrix specified by the user. Further control of algorithms can be
#'   specified with the \code{control} argument.  Elements are:
#' \describe{
#' \item{adapt_delta}{The target acceptance rate.}
#' \item{metric}{The mass metric to use. Options are: "unit" for a unit diagonal
#'   matrix; "diag" to estimate a diagonal matrix during warmup; a matrix
#'   to be used directly (in untransformed space).}
#' \item{adapt_engaged}{Whether adaptation of step size and metric is turned on.}
#' \item{max_treedepth}{Maximum treedepth for the NUTS algorithm.}
#' \item{stepsize}{The stepsize for the NUTS algorithm. If \code{NULL} it
#'   will be adapted during warmup.}
#' }
#'
#' @author Cole Monnahan
#' @param obj A TMB model object.
#' @param lower A vector of lower bounds for parameters. Allowed values are
#'   -Inf and numeric.
#' @param upper A vector of upper bounds for parameters. Allowed values are
#'   Inf and numeric.
#' @param algorithm The algorithm to use. NUTS is the default and
#'   recommended one, but "RWM" for the random walk Metropolis sampler and
#'   "HMC" for the static HMC sampler are available. These last two are
#'   deprecated but may be of use in some situations. These algorithms
#'   require different arguments; see their help files for more
#'   information.
#' @param laplace Whether to use the Laplace approximation if some
#'   parameters are declared as random. Default is to turn off this
#'   functionality and integrate across all parameters with MCMC.
#' @return A list containing the samples, and properties of the sampler
#'   useful for diagnosing behavior and efficiency.
#' @seealso \code{\link{extract_samples}} to extract samples and
#'   \code{\link{launch_shinytmb}} to explore the results graphically which
#'   is a wrapper for the \code{\link[shinystan]{launch_shinystan}} function.
#' @inheritParams sample_admb
#' @inheritSection sample_admb Warning
#' @examples
#' ## Build a fake TMB object with objective & gradient functions and some
#' ## other flags
#' \dontrun{
#' f <- function(x, order=0){
#'   if(order != 1) # negative log density
#'     -sum(dnorm(x=x, mean=0, sd=1, log=TRUE))
#'   else x # gradient of negative log density
#' }
#' init <- function() rnorm(2)
#' obj <- list(env=list(DLL='demo', last.par.best=c(x=init()), f=f,
#'   beSilent=function() NULL))
#' ## Run NUTS for this object
#' fit <- sample_tmb(obj, iter=1000, chains=3, init=init)
#' ## Check basic diagnostics
#' mon <- rstan::monitor(fit$samples, print=FALSE)
#' Rhat <- mon[,"Rhat"]
#' max(Rhat)
#' ess <- mon[, 'n_eff']
#' min(ess)
#' ## Or do it interactively with ShinyStan
#' launch_shinytmb(fit)
#' }
#'
sample_tmb <- function(obj, iter=2000, init, chains=3, seeds=NULL,
                       warmup=floor(iter/2), lower=NULL,
                       upper=NULL, thin=1, parallel=FALSE,
                       cores=NULL, path=NULL, algorithm="NUTS",
                       laplace=FALSE, control=NULL, ...){
  .Deprecated("tmbstan", package='tmbstan',
     msg='Use tmbstan::tmbstan instead, this function is no longer under development/maintenance')
  control <- .update_control(control)
  ## Argument checking.
  if(is.null(init)){
    if(chains>1) warning('Using same starting values for each chain -- strongly recommended to use dispersed inits')
    init <- lapply(1:chains, function(i) as.numeric(unlist(obj$par)))
  } else if(is.function(init)){
    init <- lapply(1:chains, function(i) unlist(init()))
  } else if(length(init) != chains){
    stop("Length of init does not equal number of chains.")
  } else if(any(unlist(lapply(init, function(x) length(unlist(x)) != length(obj$par))))){
    stop("Initial parameter vector is wrong length")
  }
  algorithm <- match.arg(algorithm, choices=c("NUTS", "RWM", "HMC"))
  stopifnot(thin >=1)
  stopifnot(chains >= 1)
  if(iter < 10 | !is.numeric(iter)) stop("iter must be > 10")
  obj$env$beSilent()                  # silence console output
  ## if(control$adapt_mass)
  ##   warning("Mass matrix adaptation is experimental -- use with caution")

  ## Ignore parameters declared as random? Borrowed from tmbstan.
  if(laplace){
    par <- obj$env$last.par.best[-obj$env$random]
    fn0 <- obj$fn
    gr0 <- obj$gr
  } else {
    par <- obj$env$last.par.best
    fn0 <- obj$env$f
    gr0 <- function(x) obj$env$f(x, order=1)
  }

  ## Parameter constraints, if provided, require the fn and gr functions to
  ## be modified to account for differents in volume. There are four cases:
  ## no constraints, bounded below, bounded above, or both (box
  ## constraint).
  bounded <- !(is.null(lower) & is.null(upper))
  if(bounded){
    if(is.null(lower)) lower <- rep(-Inf, len=length(upper))
    if(is.null(upper)) upper <- rep(Inf, len=length(lower))
    cases <- .transform.cases(lower, upper)
    fn <- function(y){
      x <- .transform(y, lower, upper, cases)
      scales <- .transform.grad(y, lower, upper, cases)
      -fn0(x) + sum(log(scales))
    }
    gr <- function(y){
      x <- .transform(y, lower, upper, cases)
      scales <- .transform.grad(y, lower, upper, cases)
      scales2 <- .transform.grad2(y, lower, upper, cases)
      -as.vector(gr0(x))*scales + scales2
    }
    init <- lapply(init, function(x) .transform.inv(x=unlist(x), a=lower, b=upper, cases=cases))
  } else {
    fn <- function(y) -fn0(y)
    gr <- function(y) -as.vector(gr0(y))
  }

  ## Make parameter names unique if vectors exist
  par.names <- names(par)
  par.names <- as.vector((unlist(sapply(unique(par.names), function(x){
    temp <- par.names[par.names==x]
    if(length(temp)>1) paste0(temp,'[',1:length(temp),']') else temp
  }))))

  ## Select and run the chain.
  if(!parallel){
    if(algorithm=="HMC"){
      mcmc.out <- lapply(1:chains, function(i)
        sample_tmb_hmc(iter=iter, fn=fn, gr=gr, init=init[[i]], warmup=warmup,
                      chain=i, thin=thin, seed=seeds[i], control=control, ...))
    } else if(algorithm=="NUTS"){
      mcmc.out <- lapply(1:chains, function(i)
        sample_tmb_nuts(iter=iter, fn=fn, gr=gr, init=init[[i]], warmup=warmup,
                      chain=i, thin=thin, seed=seeds[i], control=control, ...))
    } else if(algorithm=="RWM")
      mcmc.out <- lapply(1:chains, function(i)
        sample_tmb_rwm(iter=iter, fn=fn, init=init[[i]], warmup=warmup,
                       chain=i, thin=thin, seed=seeds[i], control=control, ...))
  } else {
    if(is.null(path))
      stop("path argument must be supplied to use parallel chains")
    if(!requireNamespace("snowfall", quietly=TRUE)) stop("snowfall package not found")
    stopifnot(is.character(path))
    if(file.exists('mcmc_progress.txt')) trash <- file.remove('mcmc_progress.txt')
    snowfall::sfInit(parallel=TRUE, cpus=cores, slaveOutfile='mcmc_progress.txt')
    ## snowfall::sfLibrary("TMB")
    snowfall::sfExportAll()
    on.exit(snowfall::sfStop())
    message("Starting parallel chains... ")
    ##mcmc.out <- lapply(1:chains, function(i)
    mcmc.out <- snowfall::sfLapply(1:chains, function(i)
      sample_tmb_parallel(parallel_number=i, iter=iter, obj=obj, path=path,
                          init=init[[i]], algorithm=algorithm, warmup=warmup,
                          lower=lower, upper=upper, seed=seeds[i],
                          laplace=laplace,
                          control=control, ...))
    message("... Finished parallel chains")
  }
  warmup <- mcmc.out[[1]]$warmup
  ## Clean up returned output
  samples <-  array(NA, dim=c(nrow(mcmc.out[[1]]$par), chains, 1+length(par.names)),
                    dimnames=list(NULL, NULL, c(par.names,'lp__')))
  ## Before transforming, get estimated covariance to be used as metrix
  ## later.
  covar.est <- cov(do.call(rbind, lapply(1:chains, function(i) mcmc.out[[i]]$par[-(1:warmup),1:length(par.names)])))
  dimnames(covar.est) <- NULL
  for(i in 1:chains){
    if(bounded){
      temp <- mcmc.out[[i]]$par
      temp[,-ncol(temp)] <-
        t(apply(temp[,-ncol(temp)], 1, function(x)
          .transform(x, lower, upper, cases)))
      samples[,i,] <- temp
    } else {
      samples[,i,] <- mcmc.out[[i]]$par
    }
  }
  ## message("... Calculating ESS and Rhat")
  ## temp <- (rstan::monitor(samples, warmup=warmup, probs=.5, print=FALSE))
  ## Rhat <- temp[,6]; ess <- temp[,5]
  sampler_params <- lapply(mcmc.out, function(x) x$sampler_params)
  time.warmup <- unlist(lapply(mcmc.out, function(x) as.numeric(x$time.warmup)))
  time.total <- unlist(lapply(mcmc.out, function(x) as.numeric(x$time.total)))
  result <- list(samples=samples, sampler_params=sampler_params,
                 time.warmup=time.warmup, time.total=time.total,
                 algorithm=algorithm, warmup=warmup,
                 model=obj$env$DLL, covar.est=covar.est)#, Rhat=Rhat, ess=ess)
  if(algorithm=="NUTS") result$max_treedepth <- mcmc.out[[1]]$max_treedepth
  return(invisible(result))
}
