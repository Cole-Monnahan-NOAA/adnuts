

#' @rdname wrappers
#' @export
sample_nuts <- function(model, path=getwd(), iter=2000, init=NULL, chains=3, warmup=NULL,
                        seeds=NULL, thin=1, mceval=FALSE, duration=NULL,
                        parallel=FALSE, cores=NULL, control=NULL,
                        skip_optimization=TRUE, verbose=TRUE,
                        skip_monitor=FALSE, skip_unbounded=TRUE,
                        admb_args=NULL, extra.args=NULL){
  ## Argument checking and processing
  if (!missing(parallel)) {
    warning("Argument parallel is deprecated, set cores=1 for serial, and cores>1 for parallel.",
            call. = FALSE)
  }
  if (!missing(extra.args)) {
    warning("Argument extra.args is deprecated, use admb_args instead",
            call. = FALSE)
    admb_args <- extra.args
  }
  if(is.null(init) & verbose)
    warning('Default init of MLE used for each chain. Consider using dispersed inits')
  .sample_admb(model=model, path=path, iter=iter, init=init,
               chains=chains, warmup=warmup, seeds=seeds,
               thin=thin, mceval=mceval, duration=duration,
               cores=cores, control=control, algorithm="NUTS",
               skip_optimization=skip_optimization,
               skip_monitor=skip_monitor,
               skip_unbounded=skip_unbounded,
               admb_args=admb_args, verbose=verbose)
}

#' @rdname wrappers
#' @export
sample_rwm <- function(model, path=getwd(), iter=2000, init=NULL, chains=3, warmup=NULL,
                        seeds=NULL, thin=1, mceval=FALSE, duration=NULL,
                        parallel=FALSE, cores=NULL, control=NULL,
                        skip_optimization=TRUE, verbose=TRUE,
                        skip_monitor=FALSE, skip_unbounded=TRUE,
                        admb_args=NULL, extra.args=NULL){
  ## Argument checking and processing
  if (!missing(parallel)) {
    warning("Argument parallel is deprecated, set cores=1 for serial, and cores>1 for parallel.",
            call. = FALSE)
  }
    if (!missing(extra.args)) {
    warning("Argument extra.args is deprecated, use admb_args instead",
            call. = FALSE)
    admb_args <- extra.args
    }
  if(is.null(init) & verbose)
    warning('Default init of MLE used for each chain. Consider using dispersed inits')
  .sample_admb(model=model, path=path, iter=iter, init=init,
               chains=chains, warmup=warmup, seeds=seeds,
               thin=thin, mceval=mceval, duration=duration,
               cores=cores, control=control, algorithm="RWM",
               skip_optimization=skip_optimization,
               skip_monitor=skip_monitor, verbose=verbose,
               skip_unbounded=skip_unbounded,
               admb_args=admb_args)
}



#' Bayesian inference of an ADMB model using the no-U-turn
#' sampler (NUTS) or random walk Metropolis (RWM) algorithms.
#'
#' Draw Bayesian posterior samples from an AD Model Builder
#' (ADMB) model using an MCMC algorithm. `sample_nuts` and
#' `sample_rwm` generates posterior samples from which inference
#' can be made.
#'
#' Adaptation schemes are used with NUTS so specifying tuning
#' parameters is not necessary. See vignette for options for
#' adaptation of step size and mass matrix.  The RWM algorithm
#' provides no new functionality not available from previous
#' versions of ADMB. However, `sample_rwm` has an improved
#' console output, is setup for parallel execution, and a smooth
#' workflow for diagnostics.
#'
#' Parallel chains will be run if argument `cores` is greater
#' than one. This entails copying the folder, and starting a new
#' R session to run that chain, which are then merged back
#' together. Note that console output is inconsistent when using
#' parallel, and may not show. On Windows the R terminal shows
#' output live, but the GUI does not. RStudio is a special case
#' and will not show live, and instead is captured and returned
#' at the end. It is strongly recommended to start with serial
#' execution as debugging parallel chains is very difficult.
#'
#' Note that the algorithm code is in the ADMB source code, and
#' 'adnuts' provides a wrapper for it. The command line arguments
#' are returned and can be examined by the user. See vignette for
#' more information.
#'
#' @details This function implements algorithm 6 of Hoffman and Gelman (2014),
#' and loosely follows package \code{rstan}. The step size can be
#'   adapted or specified manually. The metric (i.e., mass matrix) can be
#'   unit diagonal, adapted diagonal (default and recommended), a dense
#'   matrix specified by the user, or an adapted dense matrix.
#'  Further control of algorithms can be
#'   specified with the \code{control} argument.  Elements are:
#' \describe{
#' \item{adapt_delta}{The target acceptance rate, with values
#'   closer to 1 forcing smaller step sizes. Defaults to 0.8.  }
#' \item{metric}{The mass metric to use. Options are: "unit" for a unit diagonal
#'   matrix; \code{NULL} to estimate a diagonal matrix during warmup; a matrix
#'   to be used directly (in untransformed space).}
#' \item{adapt_mass}{Whether adaptation of diagonal mass matrix is turned
#'   on.}
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
#'   matrix adaptation (final fast phase)}
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
#' @param model Name of model (i.e., 'model' for model.tpl). For
#'   non-Windows systems this will automatically be converted to
#'   './model' internally. For Windows, long file names are
#'   sometimes shortened from e.g., 'long_model_filename' to
#'   'LONG_~1'. This should work, but will throw warnings. Please
#'   shorten the model name. See
#'   https://en.wikipedia.org/wiki/8.3_filename.
#' @param path Path to model executable. Defaults to working
#'   directory. Often best to have model files in a separate
#'   subdirectory, particularly for parallel.
#' @param iter The number of samples to draw.
#' @param init Can be either a list containing a vector for each
#'   chain, a function which returns a vector of parameters, or
#'   NULL which specifies to use the MLE as stored in the
#'   admodel.hes file. It is generally recommended to use
#'   dispersed initial values to improve diagnostic checks
#'   (starting from the same point makes it less likely to find
#'   multiple modes).
#' @param chains The number of chains to run.
#' @param warmup The number of warmup iterations.
#' @param seeds A vector of seeds, one for each chain.
#' @param thin The thinning rate to apply to samples. Typically
#'   not used with NUTS.
#' @param mceval Whether to run the model with \code{-mceval} on
#'   samples from merged chains.
#' @param duration The number of minutes after which the model
#'   will quit running. It is recommended to set the warmup
#'   carefully and iter higher than expected so it runs through
#'   duration. This usually results in chains with different
#'   lengths, so the minimum is taken across them all.
#' @param parallel A deprecated argument, use cores=1 for serial
#'   execution or cores>1 for parallel (default is to parallel
#'   with cores equal to the available-1)
#' @param cores The number of cores to use for parallel
#'   execution. Default is number available in the system minus
#'   1. If \code{cores=1}, serial execution occurs (even if
#'   \code{chains>1}), otherwise parallel execution via package
#'   snowfall is used. For slow analyses it is recommended to set
#'   \code{chains}<=\code{cores} so each core needs to run only a
#'   single chain.
#' @param control A list to control the sampler. See details for
#'   further use.
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
#' @param verbose Flag whether to show console output (default)
#'   or suppress it completely except for warnings and
#'   errors. Works for serial or parallel execution.
#' @param admb_args A character string which gets passed to the
#'   command line, allowing finer control
#' @param extra.args Deprecated, use a \code{admb_args} instead.
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
#' fit <- sample_nuts(model='simple', init=init, path=path)
#' unlink(path, TRUE) # cleanup folder
#' setwd(oldwd)
#' }
#'
NULL

#' Deprecated version of wrapper function. Use sample_nuts or
#' sample_rwm instead.
#'
#' @inheritParams wrappers
#' @param algorithm The algorithm to use, one of "NUTS" or "RWM"
#' @section Warning: This is deprecated and will cease to exist
#'   in future releases
#' @export
sample_admb <- function(model, path=getwd(), iter=2000, init=NULL, chains=3, warmup=NULL,
                        seeds=NULL, thin=1, mceval=FALSE, duration=NULL,
                        parallel=FALSE, cores=NULL, control=NULL,
                        skip_optimization=TRUE, algorithm='NUTS',
                        skip_monitor=FALSE, skip_unbounded=TRUE,
                        admb_args=NULL){
  ## Argument checking and processing
  if (!missing(parallel)) {
    warning("Argument parallel is deprecated, set cores=1 for serial, and cores>1 for parallel.",
            call. = FALSE)
  }
  warning("Function sample_admb is deprecated, use sample_nuts or sample_rwm instead",
          call. = FALSE)
  .sample_admb(model=model, path=path, iter=iter, init=init,
               chains=chains, warmup=warmup, seeds=seeds,
               thin=thin, mceval=mceval, duration=duration,
               cores=cores, control=control, algorithm=algorithm,
               skip_optimization=skip_optimization,
               skip_monitor=skip_monitor,
               skip_unbounded=skip_unbounded,
               admb_args=admb_args)
}


#' Hidden wrapper function for sampling from ADMB models
#'
#' @inheritParams wrappers
#' @param algorithm The algorithm to use, one of "NUTS" or "RWM"
#'
.sample_admb <- function(model, path=getwd(), iter=2000, init=NULL, chains=3, warmup=NULL,
                         seeds=NULL, thin=1, mceval=FALSE, duration=NULL,
                         cores=NULL, control=NULL, verbose=TRUE,
                         algorithm="NUTS", skip_optimization=TRUE,
                         skip_monitor=FALSE, skip_unbounded=TRUE,
                         admb_args=NULL){
  if(is.null(cores)) cores <- parallel::detectCores()-1
  cores.max  <- parallel::detectCores()
  if(cores > cores.max) {
    cores <- cores.max-1
    warning(paste('Specified cores larger than available, using total-1=', cores))
  }
  stopifnot(is.numeric(cores))
  if(!is.null(control) & !is.list(control))
    stop("control argument invalid, must be a list")
  if(cores<1) stop(paste("Cores must be >=1, but is", cores))
  parallel <- ifelse(cores==1 | chains ==1, FALSE, TRUE)
  if(parallel & cores < chains)
    if(verbose) message(paste("Recommend using chains < cores=", cores))
  stopifnot(thin >=1); stopifnot(chains >= 1)
  if(is.null(seeds)) seeds <- sample.int(1e7, size=chains)
  if(length(seeds) != chains) stop("Length of seeds must match chains")
  if(iter < 1 | !is.numeric(iter)) stop("iter must be > 1")

  ## Catch path and model name errors early
  .check_model_path(model=model, path=path)
  ## Check verison; warnings only meaningful for NUTS at the moment.
  v <- .check_ADMB_version(model=model, path=path, warn= (algorithm=='NUTS'))
  if(v<=12.0 & !skip_unbounded) {
    warning(paste('Version', v, 'of ADMB is incompatible with skip_unbounded=FALSE, ignoring'))
    skip_unbounded <- TRUE
  }

  ## Update control with defaults
  if(is.null(warmup)) warmup <- floor(iter/2)
  if(!(algorithm %in% c('NUTS', 'RWM')))
    stop("Invalid algorithm specified")
  if(algorithm=='NUTS')
    control <- .update_control(control)
  if(is.null(init)){
    ## warning moved to higher functions
  }  else
    if(is.function(init)){
    init <- lapply(1:chains, function(x) init())
  } else if(!is.list(init)){
    stop("init must be NULL, a list, or a function")
  }
  if(!is.null(init) & length(init) != chains){
    stop("Length of init does not equal number of chains.")
  }
  ## Delete any psv files in case something goes wrong we dont use old
  ## values by accident. Also windows short names might cause
  ## there to be two
  ff <- list.files(path)[grep('.psv', x=list.files(path))]
  if(length(ff)>1){
      if(.Platform$OS.type == "windows" & length(grep("~", ff))>0){
        warning("It appears a shortened Windows filename exists,",
                "which occurs with long\nmodel names. Try shortening it.",
                " See help for argument 'model'")
      } else {
        warning("Found more than one .psv file. Deleting: ", paste(ff, collapse=' '))
      }
  }
  trash <- suppressWarnings(file.remove(file.path(path, ff)))
  trash <- suppressWarnings(file.remove(file.path(path, 'adaptation.csv'),
                                        file.path(path, 'unbounded.csv')))
  ## Run in serial
  if(!parallel){
    if(algorithm=="NUTS"){
      mcmc.out <- lapply(1:chains, function(i)
        sample_admb_nuts(path=path, model=model, warmup=warmup, duration=duration,
                         iter=iter, init=init[[i]], chain=i,
                         seed=seeds[i], thin=thin, verbose=verbose,
                         control=control, admb_args=admb_args,
                         skip_optimization=skip_optimization,
                         parallel=parallel))
    } else {
      mcmc.out <- lapply(1:chains, function(i)
        sample_admb_rwm(path=path, model=model, warmup=warmup, duration=duration,
                        iter=iter, init=init[[i]], chain=i,
                        seed=seeds[i], thin=thin,
                        control=control, verbose=verbose,
                        skip_optimization=skip_optimization,
                        admb_args=admb_args,
                        parallel=parallel))
    }
    ## Parallel execution
  } else {
    console <- .check_console_printing(parallel)
    snowfall::sfStop()
    snowfall::sfInit(parallel=TRUE, cpus=cores)
    if(verbose){
      if(console)
        message("Parallel output to console is inconsistent between consoles.\n",
                "For live updates try using Rterm. See help for info on console output")
      else
        message("RStudio detected so output will display at conclusion. \n",
                "For live updates try using Rterm. See help for info on console output")
    }
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
                           control=control, verbose=verbose,
                           skip_optimization=skip_optimization,
                           admb_args=admb_args,
                           parallel=TRUE))
    if(!console & !is.null(mcmc.out[[1]]$progress)){
      trash <- lapply(mcmc.out, function(x) writeLines(x$progress))
    }
  }

  ## Build output list
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
    ## This can happen if 'duration' arg used, or if chain errors
    ## out.
    N <- floor(min(iters)/thin)
    warning(paste0("Incomplete chain lengths, iter=(",
                   paste0(iters, collapse=','),
                   "), truncating to minimum after thinning=", N))
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
  if(verbose) message(paste("Merging post-warmup chains into main folder:", path))
  samples2 <- do.call(rbind, lapply(1:chains, function(i)
    samples[-(1:warmup), i, -dim(samples)[3]]))
  .write_psv(fn=model, samples=samples2, model.path=path)
  ## These already exclude warmup
  unbounded <- do.call(rbind, lapply(mcmc.out, function(x) x$unbounded))
  oldwd <- getwd(); on.exit(setwd(oldwd))
  setwd(path)
  write.table(unbounded, file='unbounded.csv', sep=",", col.names=FALSE, row.names=FALSE)
  if(mceval){
    if(verbose) message("Running -mceval on merged chains...")
    system(paste(.update_model(model), "-mceval", admb_args), ignore.stdout=FALSE)
  }
  covar.est <- cov(unbounded)
  if(!skip_monitor){
    if(!requireNamespace("rstan", quietly = TRUE))
      stop("Package 'rstan' is required to calculate diagnostics.\n Install it and try again, or set skip_monitor=FALSE.")
    if(verbose) message('Calculating ESS and Rhat (skip_monitor=TRUE will skip)...')
    mon <- rstan::monitor(samples, warmup, print=FALSE)
  } else {
    if(verbose) message('Skipping ESS and Rhat statistics...')
    mon <- NULL
  }
  par_names <- dimnames(samples)[[3]]
  result <- list(samples=samples, sampler_params=sampler_params,
                 par_names=par_names,
                 samples_unbounded=samples.unbounded,
                 time.warmup=time.warmup, time.total=time.total,
                 algorithm=algorithm, warmup=warmup,
                 model=model, max_treedepth=mcmc.out[[1]]$max_treedepth,
                 cmd=cmd, init=init, covar.est=covar.est, mle=mle,
                 monitor=mon)
  result <- adfit(result)
  return(result)
}

