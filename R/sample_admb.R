#' Bayesian inference of an ADMB model using the no-U-turn sampler.
#'
#' Draw Bayesian posterior samples from an AD Model Builder (ADMB) model
#' using an MCMC algorithm. This function generates posterior samples from
#' which inference can be made. Adaptation schemes are used so
#' specifying tuning parameters is not necessary, and parallel
#' execution reduces overall run time.
#'
#' @details This function implements algorithm 6 of Hoffman and Gelman (2014),
#' and loosely follows package \code{rstan}. The step size can be
#'   adapted or specified manually. The metric (i.e., mass matrix) can be
#'   unit diagonal, adapted diagonal (default and recommended), or a dense
#'   matrix specified by the user. Further control of algorithms can be
#'   specified with the \code{control} argument.  Elements are:
#' \describe{
#' \item{adapt_delta}{The target acceptance rate. D}
#' \item{metric}{The mass metric to use. Options are: "unit" for a unit diagonal
#'   matrix; \code{NULL} to estimate a diagonal matrix during warmup; a matrix
#'   to be used directly (in untransformed space).}
#' \item{adapt_delta}{Whether adaptation of step size is turned on.}
#' \item{adapt_mass}{Whether adaptation of mass matrix is turned
#'   on. Currently only allowed for diagonal metric.}
#' \item{max_treedepth}{Maximum treedepth for the NUTS algorithm.}
#' \item{stepsize}{The stepsize for the NUTS algorithm. If \code{NULL} it
#'   will be adapted during warmup.}
#' }
#'
#' @author Cole Monnahan
#' @param model Name of model (i.e., model.tpl)
#' @param path Path to model executable. Defaults to working
#'   directory. Often best to have model files in a separate subdirectory,
#'   particularly for parallel.
#' @param mceval Whether to run the model with \code{-mceval} on samples
#'   from merged chains.
#' @param duration The number of minutes after which the model will quit
#'   running.
#' @param algorithm Which algorithm to use, either "NUTS" or "RWM".
#' @inheritParams sample_tmb
#' @inheritSection sample_tmb Warning
#' @export
#' @examples
#' \dontrun{
#' ## This is the packaged simple regression model
#' path.simple <- system.file('examples', 'simple', package='adnuts')
#' ## It is best to have your ADMB files in a separate folder and provide that
#' ## path, so make a copy of the model folder locally.
#' path <- 'simple'
#' dir.create(path)
#' trash <- file.copy(from=list.files(path.simple, full.names=TRUE), to=path)
#' ## Compile and run model
#' oldwd <- getwd()
#' setwd(path)
#' system('admb simple.tpl')
#' system('simple')
#' setwd('..')
#' init <- function() rnorm(2)
#' ## Run NUTS with defaults
#' fit <- sample_admb(model='simple', init=init, path=path)
#' unlink(path, TRUE) # cleanup folder
#' setwd(oldwd)
#' }
#'
sample_admb <-
  function(model, path=getwd(), iter=2000, init=NULL, chains=3, warmup=NULL,
           seeds=NULL, thin=1, mceval=FALSE, duration=NULL,
           parallel=FALSE, cores=NULL, control=NULL, algorithm="NUTS", ...){
  ## Argument checking and processing
  stopifnot(thin >=1); stopifnot(chains >= 1)
  if(is.null(seeds)) seeds <- sample.int(1e7, size=chains)
  if(iter < 10 | !is.numeric(iter)) stop("iter must be > 10")
  ## Update control with defaults
  control <- .update_control(control)
  if(is.null(warmup)) warmup <- floor(iter/2)
  if(!(algorithm %in% c('NUTS', 'RWM')))
    stop("Invalid algorithm specified")
  if(is.null(init)){
    warning('Using MLE inits for each chain -- strongly recommended to use dispersed inits')
  }  else if(is.function(init)){
    init <- lapply(1:chains, function(x) init())
  } else if(!is.list(init)){
    stop("init must be NULL, a list, or a function")
  }
  if(!is.null(init) & length(init) != chains){
    stop("Length of init does not equal number of chains.")
  }
  ## Delete any psv files in case something goes wrong we dont use old
  ## values by accident
  trash <- file.remove(list.files()[grep('.psv', x=list.files())])
  ## Run in serial
  if(!parallel){
    if(algorithm=="NUTS"){
      mcmc.out <- lapply(1:chains, function(i)
        sample_admb_nuts(path=path, model=model, warmup=warmup, duration=duration,
                         iter=iter, init=init[[i]], chain=i,
                         seed=seeds[i], thin=thin, control=control, ...))
    } else {
      mcmc.out <- lapply(1:chains, function(i)
        sample_admb_rwm(path=path, model=model, warmup=warmup, duration=duration,
                        iter=iter, init=init[[i]], chain=i,
                        seed=seeds[i], thin=thin, control=control, ...))
    }
    ## Parallel execution
  } else {
    if(!requireNamespace("snowfall", quietly=TRUE)) stop("snowfall package not found")
    snowfall::sfInit(parallel=TRUE, cpus=cores)
    snowfall::sfExportAll()
    on.exit(snowfall::sfStop())
    mcmc.out <- snowfall::sfLapply(1:chains, function(i)
      sample_admb_parallel(parallel_number=i, path=path, model=model,
                           duration=duration,
                           algorithm=algorithm,
                           iter=iter, init=init[[i]], warmup=warmup,
                           seed=seeds[i], thin=thin, control=control, ...))
  }
    warmup <- mcmc.out[[1]]$warmup
    mle <- .read_mle_fit(model=model, path=path)
    if(is.null(mle)){
      par.names <- dimnames(mcmc.out[[1]]$samples)[[2]]
      par.names <- par.names[-length(par.names)]
    } else {
      par.names <- mle$par.names
    }
    iters <- unlist(lapply(mcmc.out, function(x) dim(x$samples)[1]))
    if(any(iters!=iter/thin)){
    N <- min(iters)
    warning(paste("Variable chain lengths, truncating to minimum=", N))
  } else {
    N <- iter/thin
  }
  samples <- array(NA, dim=c(N, chains, 1+length(par.names)),
                   dimnames=list(NULL, NULL, c(par.names,'lp__')))
  for(i in 1:chains)
    samples[,i,] <- mcmc.out[[i]]$samples[1:N,]
  if(algorithm=="NUTS")
    sampler_params <-
      lapply(mcmc.out, function(x) x$sampler_params[1:N,])
  else sampler_params <- NULL
  time.warmup <- unlist(lapply(mcmc.out, function(x) as.numeric(x$time.warmup)))
  time.total <- unlist(lapply(mcmc.out, function(x) as.numeric(x$time.total)))
  cmd <- unlist(lapply(mcmc.out, function(x) x$cmd))
  if(N < warmup) warning("Duration too short to finish warmup period")
  ## When running multiple chains the psv files will either be overwritten
  ## or in different folders (if parallel is used). Thus mceval needs to be
  ## done posthoc by recombining chains AFTER thinning and warmup and
  ## discarded into a single chain, written to file, then call -mceval.
  ## Merge all chains together and run mceval
  message(paste("... Merging post-warmup chains into main folder:", path))
  samples2 <- do.call(rbind, lapply(1:chains, function(i)
    samples[-(1:warmup), i, -dim(samples)[3]]))
  .write_psv(fn=model, samples=samples2, model.path=path)
  ## These already exclude warmup
  unbounded <- do.call(rbind, lapply(mcmc.out, function(x) x$unbounded))
  oldwd <- getwd(); on.exit(setwd(oldwd))
  setwd(path)
  write.table(unbounded, file='unbounded.csv', sep=",", col.names=FALSE, row.names=FALSE)
  if(mceval){
    message("... Running -mceval on merged chains")
    system(paste(model, "-mceval -noest -nohess"), ignore.stdout=FALSE)
  }
  ## message("... Calculating ESS and Rhat")
  ## temp <- (rstan::monitor(samples, warmup=warmup, probs=.5, print=FALSE))
  ## Rhat <- temp[,6]; ess <- temp[,5]
  covar.est <- cov(unbounded)
  result <- list(samples=samples, sampler_params=sampler_params,
                 ## ess=ess, Rhat=Rhat,
                 time.warmup=time.warmup, time.total=time.total,
                 algorithm=algorithm, warmup=warmup,
                 model=model, max_treedepth=mcmc.out[[1]]$max_treedepth,
                 cmd=cmd, covar.est=covar.est, mle=mle)
  return(invisible(result))
}

#' Run a single NUTS chain for an ADMB model
#'
#' A low level function to run a single chain. Unlikely to be used by a
#' user, instead prefer \code{\link{sample_admb}}
#' @inheritParams sample_admb
#' @param seed Random seed to use.
#' @param chain Chain number, for printing purposes only.
#' @param extra.args Character string of extra command line argument to
#' pass to ADMB.
#' @param verbose Boolean for whether to print ADMB output to console.
#' @seealso \code{\link{sample_admb}}
sample_admb_nuts <- function(path, model, iter=2000,
                             init=NULL, chain=1,
                             thin=1, warmup=NULL,
                             seed=NULL, duration=NULL, control=NULL,
                             verbose=TRUE, extra.args=NULL){
  wd.old <- getwd(); on.exit(setwd(wd.old))
  setwd(path)
  ## Now contains all required NUTS arguments
  control <- .update_control(control)
  eps <- control$stepsize
  metric <- control$metric
  stopifnot(iter >= 1)
  stopifnot(warmup <= iter)
  stopifnot(duration > 0)
  stopifnot(thin >=1)
  if(is.null(warmup)) stop("Must provide warmup")
  if(thin < 1 | thin > iter) stop("Thin must be >1 and < iter")
  max_td <- control$max_treedepth
  adapt_delta <- control$adapt_delta
  adapt_mass <- control$adapt_mass

  ## Build the command to run the model
  cmd <- paste(model,"-nox -nohess -maxfn 0 -phase 1000 -nuts -mcmc ",iter)
  cmd <- paste(cmd, "-warmup", warmup, "-chain", chain)
  if(!is.null(seed)) cmd <- paste(cmd, "-mcseed", seed)
  if(!is.null(duration)) cmd <- paste(cmd, "-duration", duration)
  cmd <- paste(cmd, "-max_treedepth", max_td, "-adapt_delta", adapt_delta)
  if(!is.null(eps)) cmd <- paste(cmd, "-hyeps", eps)

  ## Three options for metric. (1) 'mle' is to use the MLE estimates in
  ## admodel.cov without mass adaptation. (2) If a matrix is passed, this
  ## is written to file admodel.cov and no adaptation is done. (3) (default)
  ## Adaptation starting with diagonal. (4) Diagonal without mass adaptation.
  if(is.null(metric)) {
    ## The default: Use mass matrix adaptating starting from unit
    ## diag. Currently the only option where mass adaptation is used.
    cmd <- paste(cmd, '-adapt_mass')
  } else if(is.matrix(metric)){
    ## User defined one will be writen to admodel.cov
    if(!require(matrixcalc)) stop("Package matrixcalc required to pass a matrix")
    cor.user <- metric/ sqrt(diag(metric) %o% diag(metric))
    if(!matrixcalc::is.positive.definite(x=cor.user))
      stop("Invalid mass matrix, not positive definite")
    .write.admb.cov(metric, hbf=1)
    warning("admodel.cov overwritten, revert admodel_original.cov if needed")
    if(adapt_mass){
      warning("Mass matrix adaptation not allowed with user-specified matrix")
      adapt_mass <- FALSE
    }
  } else if(metric=='unit') {
    ## Identity in unbounded space, without mass adapataion
    cmd <- paste(cmd, "-mcdiag")
  } else if(metric=='mle') {
    ## ADMB default so do nothing special. No adaptation, will use
    ## estimated MLE covariance matrix in unbounded space (read from
    ## admodel.cov)
  } else {
    stop("Invalid metric option")
  }
  ## Write the starting values to file. A NULL value means to use the MLE,
  ## so need to run model
  if(!is.null(init)){
    cmd <- paste(cmd, "-mcpin init.pin")
    write.table(file="init.pin", x=unlist(init), row.names=F, col.names=F)
  } else {
    ## Use MLE values which are read in from the admodel.hes file
    ## which is the default behavior
  }
  if(!is.null(extra.args)) cmd <- paste(cmd, extra.args)

  ## Run it and get results
  time <- system.time(system(cmd, ignore.stdout=!verbose))[3]
  if(!file.exists('adaptation.csv'))
    stop("NUTS output files missing. Check that ADMB version >= 12.0.")
  sampler_params<- as.matrix(read.csv("adaptation.csv"))
  unbounded <- as.matrix(read.csv("unbounded.csv", header=FALSE))
  dimnames(unbounded) <- NULL
  pars <- .get_psv(model)
  par.names <- names(pars)
  pars[,'log-posterior'] <- sampler_params[,'energy__']
  pars <- as.matrix(pars)
  ## Thin samples and adaptation post hoc for NUTS
  pars <- pars[seq(1, nrow(pars), by=thin),]
  unbounded <- unbounded[seq(1, nrow(unbounded), by=thin),]
  sampler_params <- sampler_params[seq(1, nrow(sampler_params), by=thin),]
  time.total <- time; time.warmup <- NA
  warmup <- warmup/thin
  return(list(samples=pars, sampler_params=sampler_params,
              time.total=time.total, time.warmup=time.warmup,
              warmup=warmup, max_treedepth=max_td,
              model=model, par.names=par.names, cmd=cmd,
              unbounded=unbounded))
}


#' Run a single random walk Metropolis chain for an ADMB model
#'
#' A low level function to run a single chain. Unlikely to be used by a
#' user, instead prefer \code{\link{sample_admb}}
#' @inheritParams sample_admb_nuts
#' @seealso \code{\link{sample_admb}}
sample_admb_rwm <-
  function(path, model, iter=2000, thin=1, warmup=ceiling(iter/2),
           init=NULL,  chain=1, seed=NULL, control=NULL,
           verbose=TRUE, extra.args=NULL, duration=NULL){
    wd.old <- getwd(); on.exit(setwd(wd.old))
    setwd(path)
    ## Now contains all required NUTS arguments
    control <- .update_control(control)
    metric <- control$metric
    stopifnot(iter >= 1)
    stopifnot(warmup <= iter)
    stopifnot(duration > 0)
    stopifnot(thin >=1)
    if(is.null(warmup)) stop("Must provide warmup")
    if(thin < 1 | thin > iter) stop("Thin must be >1 and < iter")

    ## Build the command to run the model
    cmd <- paste(model,"-nox -maxfn 0 -phase 1000 -nohess -rwm -mcmc",iter)
    cmd <- paste(cmd, "-mcscale", warmup, "-chain", chain)
    if(!is.null(seed)) cmd <- paste(cmd, "-mcseed", seed)
    if(!is.null(duration)) cmd <- paste(cmd, "-duration", duration)
    cmd <- paste(cmd, "-mcsave", thin)

    ## Three options for metric. NULL (default) is to use the MLE estimates
    ## in admodel.cov.  If a matrix is passed, this is written to file and
    ## no scaling is done. Option 'unit' means identity. Note: these are
    ## all in unbounded space.
    if(is.matrix(metric)){
      ## User defined one will be writen to admodel.cov
      cor.user <- metric/ sqrt(diag(metric) %o% diag(metric))
      if(!matrixcalc::is.positive.definite(x=cor.user))
        stop("Invalid mass matrix, not positive definite")
      .write.admb.cov(metric)
    } else if(is.null(metric)){
      ## NULL means default of MLE
    } else if(metric=='mle'){
      ## also use mle (i.e., do nothing)
    } else if(metric=='unit') {
      ## Identity in unbounded space
      cmd <- paste(cmd, "-mcdiag")
    } else {
      stop("Invalid metric option")
    }
    ## Write the starting values to file. A NULL value means to use the MLE,
    ## so need to run model
    if(!is.null(init)){
      cmd <- paste(cmd, "-mcpin init.pin")
      write.table(file="init.pin", x=unlist(init), row.names=F, col.names=F)
    }
    if(!is.null(extra.args)) cmd <- paste(cmd, extra.args)

    ## Run it and get results
    time <- system.time(system(cmd, ignore.stdout=!verbose))[3]
    unbounded <- as.matrix(read.csv("unbounded.csv", header=FALSE))
    dimnames(unbounded) <- NULL
    pars <- .get_psv(model)
    par.names <- names(pars)
    lp <- as.vector(read.table('rwm_lp.txt', header=TRUE)[,1])
    pars[,'log-posterior'] <- lp
    pars <- as.matrix(pars)
    ## Thinning is done interally for RWM (via -mcsave) so don't need to do
    ## it here
    time.total <- time; time.warmup <- NA
    warmup <- warmup/thin
    return(list(samples=pars, sampler_params=NULL, time.total=time.total,
                time.warmup=time.warmup, warmup=warmup,  model=model,
                par.names=par.names, cmd=cmd, unbounded=unbounded))
  }



