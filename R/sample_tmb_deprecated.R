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
#' @param ... Further arguments to be passed to samplers
#' @return A list containing the samples, and properties of the sampler
#'   useful for diagnosing behavior and efficiency.
#' @seealso \code{\link{extract_samples}} to extract samples and
#'   \code{\link{launch_shinytmb}} to explore the results graphically which
#'   is a wrapper for the \code{\link[shinystan]{launch_shinystan}} function.
#' @inheritParams sample_admb
#' @inheritSection sample_admb Warning
#' @export
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
     msg='Use tmbstan::tmbstan instead, sample_tmb is deprecated and no longer under development/maintenance')
  control <- .update_control_tmb(control)
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
  message('Calculating ESS and Rhat...')
  mon <- rstan::monitor(samples, warmup, print=FALSE)
  sampler_params <- lapply(mcmc.out, function(x) x$sampler_params)
  time.warmup <- unlist(lapply(mcmc.out, function(x) as.numeric(x$time.warmup)))
  time.total <- unlist(lapply(mcmc.out, function(x) as.numeric(x$time.total)))
  result <- list(samples=samples, sampler_params=sampler_params,
                 time.warmup=time.warmup, time.total=time.total,
                 algorithm=algorithm, warmup=warmup,
                 model=obj$env$DLL, covar.est=covar.est, monitor=mon)
  if(algorithm=="NUTS") result$max_treedepth <- mcmc.out[[1]]$max_treedepth
  result <- adfit(result)
  return(invisible(result))
}

#' Draw MCMC samples from a model posterior using the No-U-Turn (NUTS)
#' sampler with dual averaging.
#'
#' @details
#' This function implements algorithm 6 of Hoffman and Gelman
#'   (2014), which includes adaptive step sizes (\code{eps}) via an
#'   algorithm called dual averaging. It also includes an adaptation scheme
#'   to tune a diagonal mass matrix (metric) during warmup.
#'
#' These \code{fn} and \code{gr} functions must have Jacobians already
#'   applied if there are transformations used.
#'
#' @param fn A function that returns the log of the posterior density.
#' @param gr A function that returns a vector of gradients of the log of
#'   the posterior density (same as \code{fn}).
#' @param chain The chain number, for printing only.
#' @param seed The random seed to use.
#' @references
#' Hoffman and Gelman (2014). The No-U-Turn sampler: Adaptively setting
#'   path lengths in Hamiltonian Monte Carlo. J. Mach. Learn. Res.
#'   15:1593-1623.
#' @inheritParams sample_tmb
#' @seealso \code{sample_tmb}
sample_tmb_nuts <- function(iter, fn, gr, init, warmup=floor(iter/2),
                          chain=1, thin=1, seed=NULL, control=NULL){
  ## Now contains all required NUTS arguments
  if(!is.null(seed)) set.seed(seed)
  control <- .update_control_tmb(control)
  eps <- control$stepsize
  init <- as.vector(unlist(init))
  npar <- length(init)
  M <- control$metric
  if(is.null(M)) M <- rep(1, len=npar)
  if( !(is.vector(M) | is.matrix(M) | is(M,"Matrix")) )
    stop("Metric must be vector or matrix or Matrix")
  max_td <- control$max_treedepth
  adapt_delta <- control$adapt_delta
  adapt_mass <- control$adapt_mass
  ## Mass matrix adapatation algorithm arguments. Same as Stan defaults.
  w1 <- control$w1; w2 <- control$w2; w3 <- control$w3
  aws <- w2 # adapt window size
  anw <- w1+w2 # adapt next window
  if(warmup < (w1+w2+w3) & adapt_mass){
    warning("Too few warmup iterations to do mass matrix adaptation.. disabled")
    adapt_mass <- FALSE
  }

  ## Using a mass matrix means redefining what fn and gr do and
  ## backtransforming the initial value.
  rotation <- .rotate_space(fn=fn, gr=gr, M=M, y.cur=init)
  fn2 <- rotation$fn2; gr2 <- rotation$gr2
  theta.cur <- rotation$x.cur
  chd <- rotation$chd
  J <- rotation$J

  sampler_params <-
    matrix(numeric(0), nrow=iter, ncol=6, dimnames=list(NULL,
      c("accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__",
        "divergent__", "energy__")))
  ## This holds the rotated but untransformed variables ("y" space)
  theta.out <- matrix(NA, nrow=iter, ncol=npar)
  ## how many steps were taken at each iteration, useful for tuning
  j.results <- lp <- rep(NA, len=iter)
  useDA <- is.null(eps)               # whether to use DA algorithm
  if(useDA){
    epsvec <- Hbar <- epsbar <- rep(NA, length=warmup+1)
    eps <- epsvec[1] <- epsbar[1] <-
      .find.epsilon(theta=theta.cur, fn=fn2, gr=gr2, eps=.1, verbose=FALSE)
    mu <- log(10*eps)
    Hbar[1] <- 0; gamma <- 0.05; t0 <- 10; kappa <- 0.75
  } else {
    ## dummy values to return
    epsvec <- epsbar <- Hbar <- NULL
  }
  ## Start of MCMC chain
  time.start <- Sys.time()
  message('')
  message(paste('Starting NUTS at', time.start))
  for(m in 1:iter){
    ## Initialize this iteration from previous in case divergence at first
    ## treebuilding. If successful trajectory they are overwritten
    theta.minus <- theta.plus <- theta.cur
    if(is.vector(M)){
      theta.out[m,] <- chd*theta.cur
    }else if(is.matrix(M)){
      theta.out[m,] <- t(chd %*% theta.cur)
    }else if(is(M,"Matrix")){
      theta.out[m,] <- t(as.numeric(J%*%solve(chd, solve(chd, J%*%theta.cur, system="Lt"), system="Pt")))
    }else{
      stop("M not viable")
    }
    lp[m] <- if(m==1) fn2(theta.cur) else lp[m-1]
    r.cur <- r.plus <- r.minus <-  rnorm(npar,0,1)
    H0 <- .calculate.H(theta=theta.cur, r=r.cur, fn=fn2)

    ## Draw a slice variable u in log space
    logu <-
      log(runif(1)) + .calculate.H(theta=theta.cur,r=r.cur, fn=fn2)
    j <- 0; n <- 1; s <- 1; divergent <- 0
    ## Track steps and divergences; updated inside .buildtree
    info <- as.environment(list(n.calls=0, divergent=0))
    while(s==1) {
      v <- sample(x=c(1,-1), size=1)
      if(v==1){
        ## move in right direction
        res <- .buildtree(theta=theta.plus, r=r.plus, logu=logu, v=v,
                          j=j, eps=eps, H0=H0,
                          fn=fn2, gr=gr2, info=info)
        theta.plus <- res$theta.plus
        r.plus <- res$r.plus
      } else {
        ## move in left direction
        res <- .buildtree(theta=theta.minus, r=r.minus, logu=logu, v=v,
                          j=j, eps=eps, H0=H0,
                          fn=fn2, gr=gr2, info=info)
        theta.minus <- res$theta.minus
        r.minus <- res$r.minus
      }
      ## test whether to accept this state
      if(!is.finite(res$s)) res$s <- 0
      if(res$s==1) {
        if(runif(n=1, min=0,max=1) <= res$n/n){
          theta.cur <- res$theta.prime
          lp[m] <- fn2(theta.cur)
          ## Rotate parameters
          #theta.out[m,] <-
          #  if(is.vector(M)) chd*theta.cur else t(chd %*% theta.cur)
          if(is.vector(M)){
            theta.out[m,] <- chd*theta.cur
          }else if(is.matrix(M)){
            theta.out[m,] <- t(chd %*% theta.cur)
          }else if(is(M,"Matrix")){
            theta.out[m,] <- t(as.numeric(J%*%solve(chd, solve(chd, J%*%theta.cur, system="Lt"), system="Pt")))
          }else{
            stop("M not viable")
          }
        }
      }
      n <- n+res$n
      s <- as.vector(res$s*.test.nuts(theta.plus, theta.minus, r.plus, r.minus))
      if(!is.finite(s)) s <- 0
      j <- j+1
      ## Stop doubling if too many or it's diverged enough
      if(j>=max_td) {
        ## warning("j larger than max_treedepth, skipping to next m")
        break
      }
    }
    j.results[m] <- j-1

    alpha2 <- res$alpha/res$nalpha
    if(!is.finite(alpha2)) alpha2 <- 0
    ## ---------------
    ## Step size adapation with the
    ## Do the adapting of eps.
    if(useDA){
      if(m <= warmup){
        ## Adaptation during warmup:
        Hbar[m+1] <- (1-1/(m+t0))*Hbar[m] +
          (adapt_delta-alpha2)/(m+t0)
        ## If logalpha not defined, skip this updating step and use
        ## the last one.
        ## if(is.nan(Hbar[m+1])) Hbar[m+1] <- abs(Hbar[m])
        logeps <- mu-sqrt(m)*Hbar[m+1]/gamma
        epsvec[m+1] <- exp(logeps)
        logepsbar <- m^(-kappa)*logeps + (1-m^(-kappa))*log(epsbar[m])
        epsbar[m+1] <- exp(logepsbar)
        eps <- epsvec[m+1]
      } else {
        ## Fix eps for sampling period
        eps <- epsbar[warmup]
      }
    }
    ## ---------------
    ## Do the adaptation of mass matrix. The algorithm is working in X
    ## space but I need to calculate the mass matrix in Y space. So need to
    ## do this coversion in the calcs below.
    if(adapt_mass & .slow_phase(m, warmup, w1, w3)){
      ## If in slow phase, update running estimate of variances
      ## The Welford running variance calculation, see
      ## https://www.johndcook.com/blog/standard_deviation/
      if(m== w1){
        ## Initialize algorithm from end of first fast window
        m1 <- theta.out[m,]; s1 <- rep(0, len=npar); k <- 1
      } else if(m==anw){
        ## If at end of adaptation window, update the mass matrix to the estimated
        ## variances
        M <- as.numeric(s1/(k-1)) # estimated variance
        ## Update density and gradient functions for new mass matrix
        if(any(!is.finite(M))){
          warning("Non-finite estimates in mass matrix adaptation -- reverting to unit")
          M <- rep(1, length(M))
        }
        rotation <- .rotate_space(fn=fn, gr=gr, M=M,  y.cur=theta.out[m,])
        fn2 <- rotation$fn2; gr2 <- rotation$gr2; chd <- rotation$chd;
        theta.cur <- rotation$x.cur
        ## Reset the running variance calculation
        k <- 1; m1 <- theta.out[m,]; s1 <- rep(0, len=npar)
        ## Calculate the next end window. If this overlaps into the final fast
        ## period, it will be stretched to that point (warmup-w3)
        aws <- 2*aws
        anw <- .compute_next_window(m, anw, warmup, w1, aws, w3)
        ## Find new reasonable eps since it can change dramatically when M
        ## updates
        eps <- .find.epsilon(theta=theta.cur, fn=fn2, gr=gr2, eps=.1, verbose=FALSE)
        if(!is.null(control$verbose))
        print(paste(m, ": new range(M) is:",
                    round(min(M),5), round(max(M),5), ", pars",
                    which.min(M), which.max(M), ", eps=", eps))
      } else {
        k <- k+1; m0 <- m1; s0 <- s1
        ## Update M and S
        m1 <- m0+(theta.out[m,]-m0)/k
        s1 <- s0+(theta.out[m,]-m0)*(theta.out[m,]-m1)
      }
    }
    ## End of mass matrix adaptation
    ##---------------

    ## Save adaptation info.
    sampler_params[m,] <-
      c(alpha2, eps, j, info$n.calls, info$divergent, fn2(theta.cur))
    if(m==warmup) time.warmup <- difftime(Sys.time(), time.start, units='secs')
    .print.mcmc.progress(m, iter, warmup, chain)
  } ## end of MCMC loop

  ## Process the output for returning
  theta.out <- cbind(theta.out, lp)
  theta.out <- theta.out[seq(1, nrow(theta.out), by=thin),]
  warm <- warmup/thin
  sampler_params <- sampler_params[seq(1, nrow(sampler_params), by=thin),]
  ndiv <- sum(sampler_params[-(1:warm),5])
  if(ndiv>0)
    message(paste0("There were ", ndiv, " divergent transitions after warmup"))
  msg <- paste0("Final acceptance ratio=", sprintf("%.2f", mean(sampler_params[-(1:warm),1])))
  if(useDA) msg <- paste0(msg,", and target=", adapt_delta)
  message(msg)
  if(useDA) message(paste0("Final step size=", round(eps, 3),
                           "; after ", warmup, " warmup iterations"))
  time.total <- difftime(Sys.time(), time.start, units='secs')
  .print.mcmc.timing(time.warmup=time.warmup, time.total=time.total)
  return(list(par=theta.out, sampler_params=sampler_params,
              time.total=time.total, time.warmup=time.warmup,
              warmup=warm, max_treedepth=max_td))
}

## Calculate the log joint density (Hamiltonian) value for given position and
## momentum variables.
## @details This function currently assumes iid standard normal momentum
## variables.
.calculate.H <- function(theta, r, fn) fn(theta)-(1/2)*sum(r^2)
## Test whether a "U-turn" has occured in a branch of the binary tree
## created by \ref\code{.buildtree} function. Returns TRUE if no U-turn,
## FALSE if one occurred
.test.nuts <- function(theta.plus, theta.minus, r.plus, r.minus){
  theta.temp <- theta.plus-theta.minus
  res <- (crossprod(theta.temp,r.minus) >= 0) *
    (crossprod(theta.temp, r.plus) >= 0)
  return(res)
}

## A recursive function that builds a leapfrog trajectory using a balanced
## binary tree.
##
## @references This is from the No-U-Turn sampler with dual averaging
## (algorithm 6) of Hoffman and Gelman (2014).
##
## @details The function repeatedly doubles (in a random direction) until
## either a U-turn occurs or the trajectory becomes unstable. This is the
## 'efficient' version that samples uniformly from the path without storing
## it. Thus the function returns a single proposed value and not the whole
## trajectory.
##
.buildtree <- function(theta, r, logu, v, j, eps, H0, fn, gr,
                       delta.max=1000, info = environment() ){
  if(j==0){
    ## ## Useful code for debugging. Returns entire path to global env.
    ## if(!exists('theta.trajectory'))
    ##   theta.trajectory <<- data.frame(step=0, t(theta))
    ## base case, take one step in direction v
    r <- r+(v*eps/2)*gr(theta)
    theta <- theta+(v*eps)*r
    r <- r+(v*eps/2)*gr(theta)
    ## verify valid trajectory. Divergences occur if H is NaN, or drifts
    ## too from from true H.
    H <- .calculate.H(theta=theta, r=r, fn=fn)
    n <- logu <= H
    s <- logu < delta.max + H
    if(!is.finite(H) | s == 0){
     info$divergent <- 1; s <- 0
    }
    ## Acceptance ratio in log space: (Hnew-Hold)
    logalpha <- H-H0
    alpha <- min(exp(logalpha),1)
    info$n.calls <- info$n.calls + 1
    ## theta.trajectory <<-
    ##   rbind(theta.trajectory, data.frame(step=tail(theta.trajectory$step,1),t(theta)))
    return(list(theta.minus=theta, theta.plus=theta, theta.prime=theta, r.minus=r,
                r.plus=r, s=s, n=n, alpha=alpha, nalpha=1))
  } else {
    ## recursion - build left and right subtrees
    xx <- .buildtree(theta=theta, r=r, logu=logu, v=v, j=j-1, eps=eps,
                       H0=H0, fn=fn, gr=gr, info=info)
    theta.minus <- xx$theta.minus
    theta.plus <- xx$theta.plus
    theta.prime <- xx$theta.prime
    r.minus <- xx$r.minus
    r.plus <- xx$r.plus
    alpha <- xx$alpha
    nalpha <- xx$nalpha
    s <- xx$s
    if(!is.finite(s)) s <- 0
    nprime <- xx$n
    ## If it didn't fail, update the above quantities
    if(s==1){
      if(v== -1){
        yy <- .buildtree(theta=theta.minus, r=r.minus, logu=logu, v=v,
                         j=j-1, eps=eps, H0=H0,
                         fn=fn, gr=gr, info=info)
        theta.minus <- yy$theta.minus
        r.minus <- yy$r.minus
      } else {
        yy <- .buildtree(theta=theta.plus, r=r.plus, logu=logu, v=v,
                         j=j-1, eps=eps, H0=H0,
                         fn=fn, gr=gr, info=info)
        theta.plus <- yy$theta.plus
        r.plus <- yy$r.plus
      }
      ### Update elements:
      ## If both slice variables fail you get 0/0.
      nprime <- yy$n+ xx$n
      if(!is.finite(nprime)) {nprime <- 0}
      ## choose whether to keep this theta
      if(nprime>0)
        if(runif(1) <= yy$n/nprime)
          theta.prime <- yy$theta.prime
      alpha <- xx$alpha+yy$alpha
      nalpha <- xx$nalpha+yy$nalpha
      ## check for valid proposal
      b <- .test.nuts(theta.plus=theta.plus,
                      theta.minus=theta.minus, r.plus=r.plus,

                      r.minus=r.minus)
      s <- yy$s*b
    }
    return(list(theta.minus=theta.minus, theta.plus=theta.plus,
                theta.prime=theta.prime,
                r.minus=r.minus, r.plus=r.plus, s=s, n=nprime,
                alpha=alpha, nalpha=nalpha))
  }
}

## Estimate a reasonable starting value for epsilon (step size) for a given
## model, for use with Hamiltonian MCMC algorithms.
##
## This is Algorithm 4 from Hoffman and Gelman (2010) and is used in the
## dual-averaging algorithms for both HMC and NUTS to find a reasonable
## starting value.
## @title Estimate step size for Hamiltonian MCMC algorithms
## @param theta An initial parameter vector.
## @param fn A function returning the log-likelihood (not the negative of
## it) for a given parameter vector.
## @param gr A function returning the gradient of the log-likelihood of a
## model.
## @param eps A value for espilon to initiate the algorithm. Defaults to
## 1. If this is far too big the algorithm won't work well and an
## alternative value can be used.
## @return Returns the "reasonable" espilon invisible, while printing how
## many steps to reach it.
## @details The algorithm uses a while loop and will break after 50
## iterations.
##
.find.epsilon <- function(theta,  fn, gr, eps=1, verbose=TRUE){
  r <- rnorm(n=length(theta), mean=0, sd=1)
  ## Do one leapfrog step
  r.new <- r+(eps/2)*gr(theta)
  theta.new <- theta+eps*r.new
  r.new <- r.new+(eps/2)*gr(theta.new)
  H1 <- .calculate.H(theta=theta, r=r, fn=fn)
  H2 <- .calculate.H(theta=theta.new, r=r.new, fn=fn)
  a <- 2*(exp(H2)/exp(H1)>.5)-1
  ## If jumped into bad region, a can be NaN so setup algorithm to keep
  ## halving eps instead of throwing error
  if(!is.finite(a)) a <- -1
  k <- 1
  ## Similarly, keep going if there are infinite values
  while (!is.finite(H1) | !is.finite(H2) | a*H2-a*H1 > -a*log(2)) {
    eps <- (2^a)*eps
    ## Do one leapfrog step
    r.new <- r+(eps/2)*gr(theta)
    theta.new <- theta+eps*r.new
    r.new <- r.new+(eps/2)*gr(theta.new)
    H2 <- .calculate.H(theta=theta.new, r=r.new, fn=fn)
    k <- k+1
    if(k>50) {
      stop("More than 50 iterations to find reasonable eps. Model is likely misspecified or some other issue.")
    }
  }
  if(verbose) message(paste("Reasonable epsilon=", eps, "found after", k, "steps"))
  return(invisible(eps))
}


#' [Deprecated] Draw MCMC samples from a model posterior using a
#' Random Walk Metropolis (RWM) sampler.
#'
#' @details This algorithm does not yet contain adaptation of \code{alpha}
#' so some trial and error may be required for efficient sampling.
#' @param alpha The amount to scale the proposal, i.e,
#' Xnew=Xcur+alpha*Xproposed where Xproposed is generated from a mean-zero
#' multivariate normal. Varying \code{alpha} varies the acceptance rate.
#' @return A list containing samples and other metadata.
#' @inheritParams sample_tmb_nuts
#' @seealso \code{\link{sample_tmb}}
#' @references
#' Metropolis, N., Rosenbluth, A.W., Rosenbluth, M.N., Teller, A.H.,
#'   Teller, E., 1953. Equation of state calculations by fast computing
#'   machines.  J Chem Phys. 21:1087-1092.
sample_tmb_rwm <- function(iter, fn, init, alpha=1, chain=1,
                         warmup=floor(iter/2), thin=1,
                         seed=NULL, control=NULL){
  if(!is.null(seed)) set.seed(seed)
  control <- .update_control_tmb(control)
  lp <- accepted <- rep(0, length=iter)
  init <- as.vector(unlist(init))
  n.params <- length(init)
  theta.out <- matrix(NA, nrow=iter, ncol=n.params)
  ## If using covariance matrix and Cholesky decomposition, redefine
  ## these functions to include this transformation. The algorithm will
  ## work in the transformed space.
  metric <- control$metric
  if(!is.null(metric)){
    fn2 <- function(theta) fn(chd %*% theta)
    chd <- t(chol(metric))               # lower triangular Cholesky decomp.
    chd.inv <- solve(chd)               # inverse
    theta.cur <- chd.inv %*% init
  } else {
    fn2 <- fn
    theta.cur <- init
  }
  fn.cur <- fn2(theta.cur)
  time.start <- Sys.time()
  message('')
  message(paste('Starting RWM at', time.start))
  for(m in 1:iter){
    ## generate proposal
    theta.new <- theta.cur + alpha*rnorm(n=n.params, mean=0, sd=1)
    fn.new <- fn2(theta.new)
    if(log(runif(1))< fn.new-fn.cur){
      ## accept
      accepted[m] <- 1
      theta.cur <- theta.out[m,] <- theta.new
      fn.cur <- fn.new
    } else {
      ## do not accept
      theta.out[m,] <- theta.cur
    }
    lp[m] <- fn.cur
    if(m==warmup) time.warmup <- difftime(Sys.time(), time.start, units='secs')
    .print.mcmc.progress(m, iter, warmup, chain)
  } ## end of MCMC loop

  ## Back transform parameters if metric is used
  if(!is.null(metric)) {
    theta.out <- t(apply(theta.out, 1, function(x) chd %*% x))
  }
  theta.out <- cbind(theta.out, lp)
  theta.out <- theta.out[seq(1, nrow(theta.out), by=thin),]
  message(paste0("Final acceptance ratio=", sprintf("%.2f", mean(accepted[-(1:warmup)]))))
  time.total <- difftime(Sys.time(), time.start, units='secs')
  .print.mcmc.timing(time.warmup=time.warmup, time.total=time.total)
  return(list(par=theta.out, sampler_params=NULL,
              time.total=time.total, time.warmup=time.warmup, warmup=warmup/thin))
}




## These fucntions were modeled after those in Stan at:
## https://github.com/stan-dev/stan/blob/967f5717edc8a101c8712838d47e6d393ac4a3db/src/stan/mcmc/windowed_adaptation.hpp

## Helper functions for mass matrix adapatation.

## Compute the next window size
##
## @param i MCMC iteration number
## @param warmup Number of warmup iterations
## @param w1 The first adapation window (usually 75)
## @param aws The last adapt window size
## @param w3 The last adaptation window (usually 50)
## @return The next adaptation window (anw)
## @details This function calculates the size of the next window for
##   adapation. If the next window size would be too long then this is
##   extended to the end of that window.
.compute_next_window <- function(i, anw, warmup, w1, aws, w3){
  anw <- i+aws
  if(anw== (warmup-w3) ) return(anw)
  ## Check that the next anw is not too long. This will be the anw for the
  ## next time this is computed. If the next one is too long, extend this
  ## one to the very end.
  nwb <- anw+2*aws
  if(nwb >= warmup-w3){
    ## if(i != warmup-w3)
    ##   message(paste("Extending last slow window from", anw, "to", warmup-w3))
    anw <- warmup-w3
  }
  return(anw)
}

## Check whether adaptation is in the slow phase
##
## @param i MCMC iteration number
## @param warmup Number of warmup iterations
## @param w1 The first adapation window (usually 75)
## @param w3 The last adaptation window (usually 50)
## @return Bool whether in slow phase
## @details During the slow phase the mass matrix is updated in a series of
##   expanding windows. See Stan manual on adaptation.
.slow_phase <- function(i, warmup, w1, w3){
  ## After w1, before start of w3
  x1 <- i>= w1 # after initial fast window
  x2 <- i<= (warmup-w3) # but before last fast window
  x3 <- i < warmup # definitely not during sampling
  return(x1 & x2 & x3)
}

## Determine which transformation is to be used for each parameter.
##
## @details The 4 cases are (0) none, (1) lower only [a,Inf], (2) upper
##   only [-Inf, b], and (3) both lower and upper [a,b]. Each case requires
##   a different set of transformation functions.
## @param lower Vector of lower bounds, potentially -infinity for some
## @param upper Vector of upper bounds, potentially infinity for some
## @return Vector of cases, in 0:3, to be used in all other transformation
##   functions. Error checking is only done here, not in other functions.
## @seealso \code{\link{.transform}}, \code{\link{.transform.inv}},
##   \code{\link{.transform.grad}}, \code{\link{.transform.grad2}}
## @export
##
.transform.cases <- function(lower, upper){
  if(length(lower) != length(upper))
    stop("Lengths of lower and upper do not match")
  if(any(is.na(c(lower, upper))) | any(is.nan(c(lower, upper))))
    stop("Bounds must be finite or -Inf/Inf -- NA and NaN not allowed")
  if(any(lower >= upper))
    stop("Lower bound >= upper bound")
  cases <- rep(NA, length(lower))
  cases[!is.finite(lower) & !is.finite(upper)] <- 0
  cases[is.finite(lower) & !is.finite(upper)] <- 1
  cases[!is.finite(lower) & is.finite(upper)] <- 2
  cases[is.finite(lower) & is.finite(upper)] <- 3
  if(any(is.na(cases)))
    stop("Something unexpected went wrong determining the bounding functions.
 Check lower and upper.")
  return(cases)
}

## This function returns the transformed variable, x=f(y).
##
## @export
.transform <- function(y, a, b, cases){
  x <- y
  ind <- cases==1
  if(length(ind)>0)
    x[ind] <- exp(y[ind])+a[ind]
  ind <- cases==2
  if(length(ind)>0)
    x[ind] <- b[ind]-exp(y[ind])
  ind <- cases==3
  if(length(ind)>0)
   x[ind] <- a[ind]+(b[ind]-a[ind])/(1+exp(-y[ind]))
  return(x)
}

## The inverse of the transformation, y=f-1(x).
## @export
.transform.inv <- function(x, a, b, cases){
  if(any(x<a) | any(x>b)) stop("x outside limits provided -- not meaningful")
  y <- sapply(1:length(x), function(i) {
    if(cases[i]==0) return(x[i])
    else if(cases[i]==1) return(log(x[i]-a[i]))
    else if(cases[i]==2) return(log(b[i]-x[i]))
    else if(cases[i]==3) return(-log( (b[i]-x[i])/(x[i]-a[i]) ))
  })
  return(y)
}

## The absolute value of the derivative of transformation.
## @export
.transform.grad <- function(y, a, b, cases){
  x <- rep(1, length(y))
  ind <- cases %in% 1:2
  if(length(ind)>0)
    x[ind] <- exp(y[ind])
  ind <- cases==3
  if(length(ind)>0){
    tmp <- exp(-y[ind])
    x[ind] <- (b[ind]-a[ind])*tmp/(1+tmp)^2
  }
  return(x)
}

## The derivative of the log of the derivative of the transformation. I.e.,
## d/dy[log(.transform.grad(y,a,b))].
## @export
.transform.grad2 <- function(y, a, b, cases){
  x <- rep(0, len=length(y))
  ind <- cases %in% 1:2
  if(length(ind)>0)
    x[ind] <- 1
  ind <- cases==3
  if(length(ind)>0)
    x[ind] <- -1+2/(1+exp(y[ind]))
  return(x)
}


#' Draw MCMC samples from a model posterior using a static HMC sampler.
#'
#' @details This function implements algorithm 5 of Hoffman and Gelman
#'   (2014), which includes adaptive step sizes (\code{eps}) via an
#'   algorithm called dual averaging.
#' @inheritParams sample_tmb_nuts
#' @param L The number of leapfrog steps to take. The NUTS algorithm does
#'   not require this as an input. If \code{L=1} this function will perform
#'   Langevin sampling. In some contexts \code{L} can roughly be thought of
#'   as a thinning rate.
#' @param eps The step size. If a numeric value is passed, it will be used
#'   throughout the entire chain. A \code{NULL} value will initiate
#'   sampler_params of \code{eps} using the dual averaging algorithm during
#'   the first \code{warmup} steps.
#' @references \itemize{ \item{Neal, R. M. (2011). MCMC using Hamiltonian
#'   dynamics. Handbook of Markov Chain Monte Carlo.}  \item{Hoffman and
#'   Gelman (2014). The No-U-Turn sampler: Adaptively setting path lengths
#'   in Hamiltonian Monte Carlo. J. Mach. Learn. Res.  15:1593-1623.}  }
#' @seealso \code{\link{sample_tmb}}
#' @return A list containing samples ('par') and algorithm details such as
#'   step size adaptation and acceptance probabilities per iteration
#'   ('sampler_params').
#' @references
#' Hoffman and Gelman (2014). The No-U-Turn sampler: Adaptively setting
#'   path lengths in Hamiltonian Monte Carlo. J. Mach. Learn. Res.
#'   15:1593-1623.
#' @seealso \code{\link{sample_tmb}}
sample_tmb_hmc <-
  function(iter, fn, gr, init, L, eps, warmup=floor(iter/2), seed=NULL,
           chain=1,thin=1, control=NULL){
  warning("NUTS should be prefered to sHMC except in rare, specific cases")
  if(!is.null(seed)) set.seed(seed)
  control <- .update_control_tmb(control)
  metric <- control$metric
  adapt_delta <- control$adapt_delta
  ## If using covariance matrix and Cholesky decomposition, redefine
  ## these functions to include this transformation. The algorithm will
  ## work in the transformed space
  if(!is.null(metric)){
    fn2 <- function(theta) fn(chd %*% theta)
    gr2 <- function(theta) as.vector( t( gr(chd %*% theta) ) %*% chd )
    chd <- t(chol(metric))               # lower triangular Cholesky decomp.
    chd.inv <- solve(chd)               # inverse
    theta.cur <- chd.inv %*% init
  } else {
    fn2 <- fn; gr2 <- gr
    theta.cur <- init
  }
  accepted <- divergence <- lp <- rep(NA, iter)
  theta.out <- matrix(NA, nrow=iter, ncol=length(init))
  sampler_params <- matrix(numeric(0), nrow=iter, ncol=4, # holds DA info by iteration
                       dimnames=list(NULL, c("accept_stat__",
                                                "stepsize__", "int_time__", "energy__")))
  ## A NULL value for eps signifies to use the dual averaging algorithm
  useDA <- is.null(eps)
  if(useDA){
    ## Initialize the dual-averaging algorithm.
    epsvec <- Hbar <- epsbar <- rep(NA, length=warmup+1)
    eps <- epsvec[1] <- epsbar[1] <-
      .find.epsilon(theta=theta.cur, fn=fn2, gr=gr2, eps=.1, verbose=FALSE)
    mu <- log(10*eps)
    Hbar[1] <- 0; gamma <- 0.05; t0 <- 10; kappa <- 0.75
  } else {
    ## dummy values to return
    epsvec <- epsbar <- Hbar <- NULL
  }
  ## Start of MCMC chain
  time.start <- Sys.time()
  message('')
  message(paste('Starting static HMC at', time.start))
  for(m in 1:iter){
    ## Jitter step size to mitigate potential negative autocorrelations,
    ## only once fixed though
    if(useDA & m > warmup) eps <- eps*runif(1,.9,1.1)
    r.cur <- r.new <- rnorm(length(init),0,1)
    theta.new <- theta.cur
    theta.leapfrog <- matrix(NA, nrow=L, ncol=length(theta.cur))
    r.leapfrog <- matrix(NA, nrow=L, ncol=length(theta.cur))
    ## Make a half step for first iteration
    r.new <- r.new+eps*gr2(theta.new)/2
    for(i in 1:L){
      theta.leapfrog[i,] <- theta.new
      r.leapfrog[i,] <- r.new
      theta.new <- theta.new+eps*r.new
      ## Full step except at end
      if(i!=L) r.new <- r.new+eps*gr2(theta.new)
      ## If divergence, stop trajectory earlier to save computation
      if(any(!is.finite(r.new)) | any(!is.finite(theta.new))) break
    }
    ## half step for momentum at the end
    r.new <- r.new+eps*gr2(theta.new)/2
    logalpha <- -fn2(theta.cur)+fn2(theta.new)+ sum(r.cur^2)/2-sum(r.new^2)/2
    ## Numerical divergence is registered as a NaN above. In this case we
    ## want to reject the proposal, mark the divergence, and adjust the
    ## step size down if still adapting (see below).
    if(!is.finite(logalpha)){
      divergence[m] <- 1
      logalpha <- -Inf
    } else {
      divergence[m] <- 0
    }
    ## Test whether to accept or reject proposed state
    if(log(runif(1)) < logalpha){
      ## accept the proposed state
      theta.cur <- theta.new
      accepted[m] <- TRUE
    } else {
      ## otherwise reject it and stay there
      accepted[m] <- FALSE
    }
    theta.out[m,] <- theta.cur
    lp[m] <- fn(theta.cur)
    if(useDA){
      ## Do the adapting of eps.
      if(m <= warmup){
        Hbar[m+1] <-
          (1-1/(m+t0))*Hbar[m] + (adapt_delta-min(1,exp(logalpha)))/(m+t0)
        logeps <- mu-sqrt(m)*Hbar[m+1]/gamma
        epsvec[m+1] <- exp(logeps)
        logepsbar <- m^(-kappa)*logeps + (1-m^(-kappa))*log(epsbar[m])
        epsbar[m+1] <- exp(logepsbar)
        eps <- epsvec[m+1]
      } else {
        eps <- epsbar[warmup]
      }
    }
    ## Save adaptation info.
    sampler_params[m,] <- c(min(1,exp(logalpha)), eps, eps*L, fn2(theta.cur))
    if(m==warmup) time.warmup <- difftime(Sys.time(), time.start, units='secs')
    .print.mcmc.progress(m, iter, warmup, chain)
  } ## end of MCMC loop
  ## Back transform parameters if metric is used
  if(!is.null(metric)) {
    theta.out <- t(apply(theta.out, 1, function(x) chd %*% x))
  }
  theta.out <- cbind(theta.out, lp)
  theta.out <- theta.out[seq(1, nrow(theta.out), by=thin),]
  sampler_params <- sampler_params[seq(1, nrow(sampler_params), by=thin),]
  if(sum(divergence[-(1:warmup)])>0)
    message(paste0("There were ", sum(divergence[-(1:warmup)]),
                   " divergent transitions after warmup"))
  message(paste0("Final acceptance ratio=", sprintf("%.2f", mean(accepted[-(1:warmup)])),
                 " and target is ", adapt_delta))
  if(useDA) message(paste0("Final step size=", round(epsbar[warmup], 3),
                           "; after ", warmup, " warmup iterations"))
  time.total <- difftime(Sys.time(), time.start, units='secs')
  .print.mcmc.timing(time.warmup=time.warmup, time.total=time.total)
  return(list(par=theta.out, sampler_params=sampler_params,
              time.total=time.total, time.warmup=time.warmup, warmup=warmup/thin))
}


## Update the control list.
##
## @param control A list passed from \code{sample_tmb}.
## @return A list with default control elements updated by those supplied
##   in \code{control}
.update_control_tmb <- function(control){
  default <- list(adapt_delta=0.8, metric=NULL, stepsize=NULL,
                  adapt_mass=TRUE, max_treedepth=12, w1=75, w2=50, w3=25)
  if(is.matrix(control$metric) & !is.null(control$adapt_mass)){
    if(control$adapt_mass==TRUE){
      warning("Mass matrix adaptation disabled if metric is a matrix")
    }
    control$adapt_mass <- FALSE
  }
  new <- default
  if(!is.null(control))
    for(i in names(control))  new[[i]] <- control[[i]]
  if(is.matrix(new$metric)) new$adapt_mass <- FALSE
  return(new)
}
