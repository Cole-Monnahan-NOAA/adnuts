

#' @rdname wrappers
#' @export
sample_nuts <- function(model, path=getwd(), iter=2000, init=NULL, chains=3, warmup=NULL,
                        seeds=NULL, thin=1, mceval=FALSE, duration=NULL,
                        parallel=FALSE, cores=NULL, control=NULL,
                        skip_optimization=TRUE,
                        skip_monitor=FALSE, skip_unbounded=TRUE,
                        admb_args=NULL, ...){
  ## Argument checking and processing
  if (!missing(parallel)) {
    warning("Argument parallel is deprecated, set cores=1 for serial, and cores>1 for parallel.",
            call. = FALSE)
  }
  .sample_admb(model=model, path=getwd(), iter=2000, init=init, chains=chains, warmup=warmup,
              seeds=seeds, thin=thin, mceval=mceval, duration=duration,
              cores=cores, control=control,
              algorithm="NUTS", skip_optimization=skip_optimization,
              skip_monitor=skip_monitor, skip_unbounded=skip_unbounded,
              admb_args=admb_args, ...)
}

#' @rdname wrappers
#' @export
sample_rwm <- function(model, path=getwd(), iter=2000, init=NULL, chains=3, warmup=NULL,
                        seeds=NULL, thin=1, mceval=FALSE, duration=NULL,
                        parallel=FALSE, cores=NULL, control=NULL,
                        skip_optimization=TRUE,
                        skip_monitor=FALSE, skip_unbounded=TRUE,
                        admb_args=NULL, ...){
  ## Argument checking and processing
  if (!missing(parallel)) {
    warning("Argument parallel is deprecated, set cores=1 for serial, and cores>1 for parallel.",
            call. = FALSE)
  }
  .sample_admb(model=model, path=getwd(), iter=2000, init=init, chains=chains, warmup=warmup,
              seeds=seeds, thin=thin, mceval=mceval, duration=duration,
              cores=cores, control=control,
              algorithm="RWM", skip_optimization=skip_optimization,
              skip_monitor=skip_monitor, skip_unbounded=skip_unbounded,
              admb_args=admb_args, ...)
}



#' Bayesian inference of an ADMB model using the no-U-turn
#' sampler (NUTS) or random walk Metropolis (RWM) algorithms.
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
#'   unit diagonal, adapted diagonal (default and recommended), a dense
#'   matrix specified by the user, or an adapted dense matrix.
#'  Further control of algorithms can be
#'   specified with the \code{control} argument.  Elements are:
#' \describe{
#' \item{adapt_delta}{The target acceptance rate. D}
#' \item{metric}{The mass metric to use. Options are: "unit" for a unit diagonal
#'   matrix; \code{NULL} to estimate a diagonal matrix during warmup; a matrix
#'   to be used directly (in untransformed space).}
#' \item{adapt_delta}{Whether adaptation of step size is turned on.}
#' \item{adapt_mass}{Whether adaptation of mass matrix is turned
#'   on. Currently only allowed for diagonal metric.}
#' \item{adapt_mass_dense}{Whether dense adaptation of mass
#'   matrix is turned on.}
#' \item{max_treedepth}{Maximum treedepth for the NUTS algorithm.}
#' \item{stepsize}{The stepsize for the NUTS algorithm. If \code{NULL} it
#'   will be adapted during warmup.}
#' \item{adapt_init_buffer}{The initial buffer size during mass matrix
#'   adaptation where sample information is not used (default
#'   50)}
#' \item{adapt_term_buffer}{The terminal buffer size (default 75)
#'   during mass
#'   matrix adpatation (final fast phase)}
#' \item{adapt_window}{The initial size of the mass matrix
#'   adaptation window, which gets doubled each time thereafter.}
#' \item{refresh}{The rate at which to refresh progress to the
#'   console. Defaults to even 10%. A value of 0 turns off
#'   progress updates.}
#' }
#' The adaptation scheme (step size and mass matrix) is based heavily on those by the
#'   software Stan, and more details can be found in that
#'   documentation and this vignette.
#'
#' @author Cole Monnahan
#' @name wrappers
#' @param model Name of model (i.e., model.tpl)
#' @param path Path to model executable. Defaults to working
#'   directory. Often best to have model files in a separate
#'   subdirectory, particularly for parallel.
#' @param iter The number of samples to draw.
#' @param init A list of lists containing the initial parameter
#'   vectors, one for each chain or a function. It is strongly
#'   recommended to initialize multiple chains from dispersed
#'   points. A of NULL signifies to use the starting values
#'   present in the model (i.e., \code{obj$par}) for all chains.
#' @param chains The number of chains to run.
#' @param warmup The number of warmup iterations.
#' @param seeds A vector of seeds, one for each chain.
#' @param thin The thinning rate to apply to samples. Typically
#'   not used with NUTS.
#' @param mceval Whether to run the model with \code{-mceval} on
#'   samples from merged chains.
#' @param duration The number of minutes after which the model
#'   will quit running.
#' @param algorithm Which algorithm to use, either "NUTS" or
#'   "RWM".
#' @param parallel A boolean for whether to use parallel
#'   cores. The package snowfall is used if TRUE.
#' @param cores The number of cores to use for parallel
#'   execution.
#' @param control A list to control the sampler. See details for
#'   further use.
#' @param algorithm The algorithm to use."NUTS" is the default
#'   and recommended one, but "RWM" for the random walk
#'   Metropolis sampler and "HMC" for the static HMC sampler are
#'   available. HMC is deprecated but may be of use in special
#'   situations. These algorithms require different arguments;
#'   see their help files for more information.
#' @param skip_optimization Whether to run the optimizer before
#'   running MCMC. This is rarely need as it is better to run it
#'   once before to get the covariance matrix, or the estimates
#'   are not needed with adaptive NUTS.
#' @param skip_monitor Whether to skip calculating diagnostics
#'   (effective sample size, Rhat) via the \code{rstan::monitor}
#'   function. This can be slow for models with high dimension or
#'   many iterations. The result is used in plots and summaries
#'   so it is recommended to turn on. If model run with
#'   \code{skip_monitor=FALSE} you can recreate it post-hoc by
#'   setting \code{fit$monitor=rstan::monitor(fit$samples,
#'   fit$warmup, print=FALSE)}.
#' @param skip_unbounded Whether to skip returning the unbounded
#'   version of the posterior samples in addition to the bounded
#'   ones. It may be advisable to set to FALSE for very large
#'   models to save space.
#' @param admb_args A character string which gets passed to the
#'   command line, allowing finer control
#' @param ... Further arguments to be passed to the
#'   algorithm. See help files for the samplers for further
#'   arguments.
#' @section Warning: The user is responsible for specifying the
#'   model properly (priors, starting values, desired parameters
#'   fixed, etc.), as well as assessing the convergence and
#'   validity of the resulting samples (e.g., through the
#'   \code{coda} package), or with function
#'   \code{\link{launch_shinytmb}} before making
#'   inference. Specifically, priors must be specified in the
#'   template file for each parameter. Unspecified priors will be
#'   implicitly uniform.
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
NULL

#' Deprecated version of wrapper function. Use sample_nuts or
#' sample_rwm instead.
#'
#' @inheritParams wrappers
#' @section Warning: This is deprecated and will cease to exist
#'   in future releases
sample_admb <- function(model, path=getwd(), iter=2000, init=NULL, chains=3, warmup=NULL,
                        seeds=NULL, thin=1, mceval=FALSE, duration=NULL,
                        parallel=FALSE, cores=NULL, control=NULL,
                        skip_optimization=TRUE,
                        skip_monitor=FALSE, skip_unbounded=TRUE,
                        admb_args=NULL, ...){
  ## Argument checking and processing
  if (!missing(parallel)) {
    warning("Argument parallel is deprecated, set cores=1 for serial, and cores>1 for parallel.",
            call. = FALSE)
  }
  warning("Function sample_admb is deprecated, use sample_nuts or sample_rwm instead",
          call. = FALSE)
  .sample_admb(model=model, path=getwd(), iter=2000, init=init, chains=chains, warmup=warmup,
              seeds=seeds, thin=thin, mceval=mceval, duration=duration,
              cores=cores, control=control,
              algorithm="NUTS", skip_optimization=skip_optimization,
              skip_monitor=skip_monitor, skip_unbounded=skip_unbounded,
              admb_args=admb_args, ...)
}


#' Hidden wrapper function for sampling from ADMB models
#'
#' @inheritParams wrappers
#'
.sample_admb <- function(model, path=getwd(), iter=2000, init=NULL, chains=3, warmup=NULL,
                        seeds=NULL, thin=1, mceval=FALSE, duration=NULL,
                        parallel=FALSE, cores=NULL, control=NULL,
                        algorithm="NUTS", skip_optimization=TRUE,
                        skip_monitor=FALSE, skip_unbounded=TRUE,
                        admb_args=NULL,
                        ...){
  if(is.null(cores)) cores <- parallel::detectCores()-1
  cores.max  <- parallel::detectCores()
  if(cores > cores.max) {
    cores <- cores.max-1
    warning(paste('Specified cores larger than available, using total-1=', cores))
  }
  stopifnot(is.numeric(cores))
  if(cores<1) stop(paste("Cores must be >=1, but is", cores))
  parallel <- ifelse(cores==1 | chains ==1, FALSE, TRUE)
  if(parallel & cores < chains)
    message(paste("Recommend using chains < cores=", cores))
  stopifnot(thin >=1); stopifnot(chains >= 1)
  if(is.null(seeds)) seeds <- sample.int(1e7, size=chains)
  if(iter < 1 | !is.numeric(iter)) stop("iter must be > 1")
  ## Catch path and model name errors early
  stopifnot(is.character(path)); stopifnot(is.character(model))
  if(!dir.exists(path)) stop(paste('Folder', path, 'does not exist. Check argument \'path\''))
  if (.Platform$OS.type=="windows") {
    ff <- file.path(path, paste(model,".exe",sep=""))
  } else {
    ff <- file.path(path, paste("./",model,sep=""))
  }
  if(!file.exists(ff)) stop(paste('File', ff, 'not found. Check \'path\' and \'model\' arguments'))
  v <- .check_ADMB_version(model=model, path=path, warn= algorithm=='NUTS')
  if(v<=12.0 & !skip_unbounded) {
    warning(paste('Version', v, 'of ADMB is incompatible with skip_unbounded=FALSE, ignoring'))
    skip_unbounded <- TRUE
  }
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
  trash <- suppressWarnings(file.remove(list.files(path)[grep('.psv', x=list.files())]))
  trash <- suppressWarnings(file.remove(file.path(path, 'adaptation.csv'),
                       file.path(path, 'unbounded.csv')))
  ## Run in serial
  if(!parallel){
    if(algorithm=="NUTS"){
      mcmc.out <- lapply(1:chains, function(i)
        sample_admb_nuts(path=path, model=model, warmup=warmup, duration=duration,
                         iter=iter, init=init[[i]], chain=i,
                         seed=seeds[i], thin=thin,
                         control=control, admb_args=admb_args,
                         skip_optimization=skip_optimization,
                         ...))
    } else {
      mcmc.out <- lapply(1:chains, function(i)
        sample_admb_rwm(path=path, model=model, warmup=warmup, duration=duration,
                        iter=iter, init=init[[i]], chain=i,
                        seed=seeds[i], thin=thin,
                        control=control,
                        skip_optimization=skip_optimization,
                        admb_args=admb_args, ...))
    }
    ## Parallel execution
  } else {
    snowfall::sfStop()
    snowfall::sfInit(parallel=TRUE, cpus=cores)
    ## errors out with empty workspace
    if(length(ls(envir=globalenv()))>0)
      snowfall::sfExportAll()
    on.exit(snowfall::sfStop())
    mcmc.out <- snowfall::sfLapply(1:chains, function(i)
      sample_admb_parallel(parallel_number=i, path=path, model=model,
                           duration=duration,
                           algorithm=algorithm,
                           iter=iter, init=init[[i]], warmup=warmup,
                           seed=seeds[i], thin=thin,
                           control=control,
                           skip_optimization=skip_optimization,
                           admb_args=admb_args, ...))
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
  if(!skip_unbounded){
    samples.unbounded <- samples
  } else {
    samples.unbounded= NULL
  }
  for(i in 1:chains){
    samples[,i,] <- mcmc.out[[i]]$samples[1:N,]
    if(!skip_unbounded)
      samples.unbounded[,i,] <- cbind(mcmc.out[[i]]$unbounded[1:N,],
                                      mcmc.out[[i]]$samples[,1+length(par.names)])
  }
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
  covar.est <- cov(unbounded)
  if(!skip_monitor){
    if(!requireNamespace("rstan", quietly = TRUE))
      stop("Package 'rstan' is required to calculate diagnostics.\n Install it and try again, or set skip_monitor=FALSE.")
    message('Calculating ESS and Rhat (skip_monitor=TRUE will skip)...')
    mon <- rstan::monitor(samples, warmup, print=FALSE)
  } else {
    message('Skipping ESS and Rhat statistics..')
    mon <- NULL
  }
  result <- list(samples=samples, sampler_params=sampler_params,
                 samples_unbounded=samples.unbounded,
                 time.warmup=time.warmup, time.total=time.total,
                 algorithm=algorithm, warmup=warmup,
                 model=model, max_treedepth=mcmc.out[[1]]$max_treedepth,
                 cmd=cmd, init=init, covar.est=covar.est, mle=mle,
                 monitor=mon)
  result <- adfit(result)
  return(invisible(result))
}

#' Run a single NUTS chain for an ADMB model
#'
#' A low level function to run a single chain. Unlikely to be used by a
#' user, instead prefer \code{\link{sample_admb}}
#' @inheritParams wrappers
#' @param seed Random seed to use.
#' @param chain Chain number, for printing purposes only.
#' @param admb_args Character string of extra command line argument to
#' pass to ADMB.
#' @param extra.args Deprecated, use \code{admb_args} instead
#' @param verbose Boolean for whether to print ADMB output to console.
#' @seealso \code{\link{sample_admb}}
sample_admb_nuts <- function(path, model, iter=2000,
                             init=NULL, chain=1,
                             thin=1, warmup=NULL,
                             seed=NULL, duration=NULL,
                             control=NULL,
                             skip_optimization=TRUE,
                             verbose=TRUE, admb_args=admb_args, extra.args=NULL){

  if (!missing(extra.args)) {
    warning("Argument extra.args is deprecated, use admb_args instead",
            call. = FALSE)
    admb_args <- extra.args
  }

  wd.old <- getwd(); on.exit(setwd(wd.old))
  setwd(path)
  ## Now contains all required NUTS arguments
  control <- .update_control(control)
  eps <- control$stepsize
  stopifnot(iter >= 1)
  stopifnot(warmup <= iter)
  stopifnot(duration > 0)
  stopifnot(thin >=1)
  if(is.null(warmup)) stop("Must provide warmup")
  if(thin < 1 | thin > iter) stop("Thin must be >1 and < iter")
  max_td <- control$max_treedepth
  adapt_delta <- control$adapt_delta

  ## Build the command to run the model
  if(skip_optimization){
    cmd <- paste(model,"-nox -nohess -maxfn 0 -phase 1000 -nuts -mcmc ",iter)
  } else {
    cmd <- paste(model,"-hbf -nuts -mcmc ",iter)
  }
  cmd <- paste(cmd, "-warmup", warmup, "-chain", chain)
  if(!is.null(seed)) cmd <- paste(cmd, "-mcseed", seed)
  if(!is.null(duration)) cmd <- paste(cmd, "-duration", duration)
  cmd <- paste(cmd, "-max_treedepth", max_td, "-adapt_delta", adapt_delta)
  if(!is.null(eps)) cmd <- paste(cmd, "-hyeps", eps)
  if(!is.null(control$adapt_init_buffer))
    cmd <- paste(cmd, "-adapt_init_buffer", control$adapt_init_buffer)
  if(!is.null(control$adapt_term_buffer))
    cmd <- paste(cmd, "-adapt_term_buffer", control$adapt_term_buffer)
  if(!is.null(control$adapt_window))
    cmd <- paste(cmd, "-adapt_window", control$adapt_window)
  if(!is.null(control$refresh))
    cmd <- paste(cmd, "-refresh", control$refresh)
  if(control$adapt_mass)
    cmd <- paste(cmd, "-adapt_mass")
  if(control$adapt_mass_dense)
    cmd <- paste(cmd, "-adapt_mass_dense")

  ## Three options for metric. (1) 'mle' is to use the MLE estimates in
  ## admodel.cov without mass adaptation. (2) If a matrix is passed, this
  ## is written to file admodel.cov and no adaptation is done. (3) (default)
  ## Adaptation starting with diagonal. (4) Diagonal without mass adaptation.
  metric <- control$metric
  stopifnot(!is.null(metric))
  if(is.matrix(metric)){
    ## User defined one will be writen to admodel.cov
    if(!requireNamespace("matrixcalc", quietly = TRUE))
      stop("Package 'matrixcalc' is required to pass a matrix.\n Install it and try again.")
    cor.user <- metric/ sqrt(diag(metric) %o% diag(metric))
    if(!matrixcalc::is.positive.definite(x=cor.user))
      stop("Invalid mass matrix passed: it is not positive definite.\n Check 'metric' argument or use different option.")
    .write.admb.cov(metric, hbf=1)
    warning("admodel.cov overwritten, revert admodel_original.cov if needed")
  } else if(is.character(metric) && metric == 'unit') {
    ## The default: Start from unit diag.
    cmd <- paste(cmd, '-mcdiag')
  } else if(is.character(metric) && metric=='mle') {
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
  if(!is.null(admb_args)) cmd <- paste(cmd, admb_args)

  ## Run it and get results
  time <- system.time(system(cmd, ignore.stdout=!verbose))[3]
  if(!file.exists('adaptation.csv') | !file.exists('unbounded.csv'))
    stop(paste0("NUTS failed to run in chain ", chain, ". Check inputs."))
  sampler_params <- as.matrix(read.csv("adaptation.csv"))
  unbounded <- as.matrix(read.csv("unbounded.csv", header=FALSE))
  dimnames(unbounded) <- NULL
  pars <- .get_psv(model)
  par.names <- names(pars)
  if(!"lp__" %in% dimnames(sampler_params)[[2]]){
    ## Previous version had a bug where energy__ was stored as
    ## the log-posterior. So energy is wrong, but log-posterior
    ## is right here.
    ## warning("ADMB version <= 12.0 has a bug where the energy statistic is wrong. Please consider updating")
    pars[,'log-posterior'] <- sampler_params[,'energy__']
  } else {
    ## Later versions has a 7th column containing the LP and 6 is
    ## the energy. Both enegy and lp are correct
    pars[,'log-posterior'] <- sampler_params[,'lp__']
    ## Drop the lp__ here since not used and may cause issues
    ## downstream.
    sampler_params <- sampler_params[,-7]
  }
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
#' @inheritParams wrappers
#' @seealso \code{\link{sample_admb}}
sample_admb_rwm <-
  function(path, model, iter=2000, thin=1, warmup=ceiling(iter/2),
           init=NULL,  chain=1, seed=NULL, control=NULL,
           verbose=TRUE, duration=NULL, extra.args=NULL,
           admb_args=NULL, skip_optimization=TRUE){

    if (!missing(extra.args)) {
      warning("Argument extra.args is deprecated, use admb_args instead",
              call. = FALSE)
      admb_args <- extra.args
    }

    wd.old <- getwd(); on.exit(setwd(wd.old))
    setwd(path)
    ## Now contains all required NUTS arguments
    control <- .update_control(control)
    metric <- 'mle' ## only one allowed
    stopifnot(iter >= 1)
    stopifnot(warmup <= iter)
    stopifnot(duration > 0)
    stopifnot(thin >=1)
    if(is.null(warmup)) stop("Must provide warmup")
    if(thin < 1 | thin > iter) stop("Thin must be >1 and < iter")

    ## Build the command to run the model
    if(skip_optimization){
      cmd <- paste(model,"-nox -nohess -maxfn 0 -phase 1000 -rwm -mcmc ",iter)
    } else {
      cmd <- paste(model,"-rwm -mcmc ",iter)
    }

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
    if(!is.null(admb_args)) cmd <- paste(cmd, admb_args)

    ## Run it and get results
    time <- system.time(system(cmd, ignore.stdout=!verbose))[3]
    if(!file.exists('unbounded.csv'))
      stop(paste0("RWM failed to run in chain ", chain, ". Check inputs."))
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



