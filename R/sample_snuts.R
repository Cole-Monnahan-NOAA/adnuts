
#' NUTS sampling for TMB models using a sparse metric (BETA).
#'
#' @param obj The TMB object with random effects turned on and
#'   optimized.
#' @param num_samples The number of post-warmup iterations to
#' run per chain.
#' @param num_warmup The number of warmup iterations to run per
#'   chain. The default of NULL indicates to automatically
#'   determine it based on other settings (recommended).
#' @param metric A character specifying which metric to use.
#'   Defaults to "auto" which uses an algorithm to select the
#'   best metric (see details), otherwise one of "sparse",
#'   "dense", "diag", "unit", "Stan", or "sparse-naive" can be
#'   specified.
#' @param init Either 'last.par.best' (default), 'random',
#'   'random-t', or 'unif'. The former starts from the joint
#'   mode, while 'random' and 'random-t' draw from multivariate
#'   normal or multivariate t with 2 degrees of freedom
#'   distributions using the inverse joint precision matrix as a
#'   covariance matrix. 'random-t' is provided to allow for more
#'   dispersed initial values. 'unif' will draw U(-2,2) samples
#'   for all parameters, similar ot Stan's default behavior. If
#'   the joint NLL is undefined at the initial values then the
#'   model will exit and return the initial vector for further
#'   investigation by the user, if desired. Note that
#'   \code{\link{StanEstimators::stan_sample}} only allows for
#'   the same init vector for all chains currently. If a seed is
#'   specified it will be set and thus the inits used will be
#'   reproducible. The inits are also returned in the 'inits'
#'   slot of the fitted object.
#' @param chains Number of chains
#' @param cores Number of parallel cores to use, defaults to
#'   \code{chains} so set this to 1 to execute serial chains.
#' @param thin The thinning rate (defaults to 1).
#' @param adapt_stan_metric A boolean whether Stan should engage
#'   diagonal mass matrix adaptation. The default of NULL
#'   indicates to automatically select depending on other
#'   settings. See details.
#' @param control NUTS control list, currently available options
#'   are 'adapt_delta', 'max_treedepth', and 'metric' which is
#'   the type of metric adaptation for Stan to do with options
#'   ('unit_e', 'diag_e', or 'dense_e'). For dense and sparse
#'   metrics this usually can be 'unit_e' to skip adaptation.
#'   NULL values (default) revert to \code{stan_sample} defaults.
#' @param seed Random number seed, used for generating initial
#'   values (if 'random") and for NUTS.
#' @param laplace Whether to leave the Laplace approximation on
#'   and only use NUTS to sample the fixed effects, or turn it
#'   off and sample from the joint parameter space (default).
#' @param skip_optimization Whether to skip optimization or not
#'   (default).
#' @param Q The sparse precision matrix. It is calculated
#'   internally if not specified (default).
#' @param Qinv The dense matrix (M). It is calculated internally
#'   if not specified (default).
#' @param globals A named list of objects to pass to new R
#'   sessions when running in parallel and using RTMB. Typically
#'   this is the `data` object for now.
#' @param model_name An optional character giving the model name.
#'   If NULL it will use the DLL name which for RTMB models is
#'   just 'RTMB'. The name is used only for printing.
#' @param refresh How often to print updates to console
#'   (integer). 0 will turn off printing. The default is 100.
#' @param print Whether to print summary of run (default) or not
#' @param rotation_only Whether to return only the rotation object
#'  (for debugging purposes)
#' @param iter (Deprecated) Total iterations to run (warmup + sampling)
#' @param warmup (Deprecated) Total warmup iterations. Defaults to
#'   \code{iter}/2 based on Stan defaults, but when using dense,
#'   sparse, or diag metrics a much shorter warmup can be used
#'   (e.g., 150), especially if paired with a 'unit_e' Stan
#'   metric. Use \code{\link{plot_sampler_params}} to investigate
#'   warmup performance and adjust as necessary for subsequent
#'   runs.
#' @param ... Additional arguments to pass to
#'   \code{\link{StanEstimators::stan_sample}}.
#' @return A fitted MCMC object of class 'adfit'
#'
#' @details
#' The \strong{TMB metric} is used to decorrelate and descale the posterior
#' distribution before sampling with Stan's algorithms. The chosen metric
#' and the reasoning for its selection are printed to the console before
#' sampling begins.
#'
#' \strong{Metric Options}
#' \itemize{
#'  \item{\code{'auto'}: This is the default setting. It uses an internal
#'    algorithm to determine the optimal metric for the model. The choice
#'    depends on the availability of the precision matrix (\eqn{Q}) and/or
#'    the covariance matrix (\eqn{M=Q^{-1}}), the extent of parameter
#'    correlations, and the speed of gradient calculations.}
#'  \item{\code{'dense'} and \code{'sparse'}: Both of these options
#'    decorrelate and descale the posterior. However, the
#'    \code{'sparse'} metric is more computationally efficient
#'    for models with high-dimensional, sparse precision matrices.
#'    For models without random effects the \code{'sparse'} option
#'    is not available}
#'  \item{\code{'diag'}: This option only descales the posterior, using the
#'    marginal standard deviations derived from the covariance matrix
#'    \eqn{M}. It does not account for correlations.}
#'  \item{\code{'unit'}: This option uses an identity matrix, which is the
#'    default behavior in Stan. Unlike the \code{'Stan'} option below, the
#'    model is still optimized to find the mode, and the \eqn{Q} matrix is
#'    calculated. This ensures that a full \code{mle} object (containing the
#'    mode, standard errors, and correlations) is returned.}
#'  \item{\code{'sparse-naive'}: This metric is constructed to be
#'    mathematically equivalent to \code{'dense'} but is often
#'    computationally faster. It is generally recommended only for testing
#'    and exploration.}
#'  \item{\code{'stan'}: This is a special flag that reverts the sampler to
#'    the standard Stan behavior. It skips the optimization and all \eqn{Q}
#'    matrix calculations and ensures that Stan's mass matrix adaptation is
#'    engaged during warmup.}
#' }
#'
#' \strong{Important Distinction}
#'
#' Note that the \code{metric} parameter described here is
#' specific to \strong{TMB} and is distinct from the Stan metric,
#' which is controlled via the \code{control} list argument in
#' the sampling function. Stan by default adapts a diagonal mass
#' matrix (metric_e) using a series of expanding windows. If Q is
#' a good estimate of the global covariance then this is not
#' needed and disabling Stan's metric adaptation is recommended.
#' This can be done by setting `adapt_stan_metric=FALSE`. If left
#' at NULL Stan's adaptation will only be done for metrics 'stan'
#' and 'unit' because those two options do not descale the
#' posterior. In this case, it is recommended to use a longer
#' warmup period to account for this adaptive procedure.
#' @export
sample_snuts <-
  function(obj, num_samples=1000, num_warmup=NULL,
           chains=4, cores=chains, thin=1,
           adapt_stan_metric=NULL,
           control=NULL, seed=NULL, laplace=FALSE,
           init=c('last.par.best', 'random', 'random-t', 'unif'),
           metric=c('auto', 'unit', 'diag', 'dense',
                    'sparse', 'stan', 'sparse-naive'),
           skip_optimization=FALSE, Q=NULL, Qinv=NULL,
           globals=NULL, model_name=NULL, refresh=NULL,
           print=TRUE, rotation_only=FALSE,
           iter=2000, warmup=floor(iter/2),
           ...){
    if (!requireNamespace("StanEstimators", quietly = TRUE))
      stop("StanEstimators package must be installed manually to use this function. \nSee https://github.com/andrjohns/StanEstimators")
    if(!is.list(obj)) stop("obj appears not to be a valid TMB object")
    if(is.null(obj$env$DLL)) stop("obj appears not to be a valid TMB object")
    if(!missing(iter)){
      warning("The arguments iter and warmup are deprecated in favor of num_samples and num_warmup")
      if(is.null(warmup)) warmup <- floor(iter/2)
      num_samples <- iter-warmup
      num_warmup <- warmup
    }
    metric <- match.arg(metric)
    init <- match.arg(init)
    obj$env$beSilent()
    if(!is.null(model_name)){
      stopifnot(is.character(model_name))
    } else {
      model_name <- obj$env$DLL
    }
    if(metric=='stan'){
      skip_optimization <- TRUE
      laplace <- FALSE
    }
    # This function optimizes and gets Q and Qinv, depending on
    # metric, except 'stan' which skips that all.
    inputs <- .get_inputs(obj=obj, skip_optimization=skip_optimization,
                          laplace=laplace, metric=metric, Q=Q, Qinv=Qinv)
    # build new joint object here so can test initial values next
    mydll <- unclass(getLoadedDLLs()[[obj$env$DLL]])$path
    isRTMB <- ifelse(obj$env$DLL=='RTMB', TRUE, FALSE)
    if(!isRTMB){
      packages <- c("TMB", "Matrix")
      obj2 <- obj
      if(!laplace){
        message("Rebuilding TMB obj without random effects...")
        obj2 <- TMB::MakeADFun(data=obj$env$data, parameters=obj$env$parList(),
                               map=obj$env$map,
                               random=NULL, silent=TRUE,
                               DLL=obj$env$DLL)
      }
    } else {
      packages <- c("RTMB", "Matrix")
      obj2 <- obj
      if(!laplace){
        message("Rebuilding RTMB obj without random effects...")
        obj2 <- RTMB::MakeADFun(func=obj$env$data, parameters=obj$env$parList(),
                                map=obj$env$map,
                                random=NULL, silent=TRUE,
                                DLL=obj$env$DLL)
      }
    }
    # reset back to the mode
    if(!is.null(inputs$mle$est)) dummy <- obj2$fn(inputs$mle$est)
    # this is the initial value in the untransformed (original)
    # parameter space
    yinits <- .get_inits(init=init, obj2=obj2,
                        seed=seed, inputs=inputs)
    # now can transform the parameter space via the metric selected
    rotation <- .rotate_posterior(metric=metric, fn=obj2$fn,
                                  gr=obj2$gr, Q=inputs$Q,
                                  Qinv=inputs$Qinv,
                                  y.cur=yinits)
    if(rotation_only) return(rotation)
    fsparse <- function(v) {dyn.load(mydll); -rotation$fn2(v)}
    gsparse <- function(v) -as.numeric(rotation$gr2(v))
    inits <- rotation$x.cur
    if(!is.null(inputs$mle$est)){
      nll0=round(-obj2$fn(inputs$mle$est),3)
      message("log-posterior at inits=", round(fsparse(inits),3),
              "; at conditional mode=",nll0)
    } else {
      message("log-posterior at inits=", round(fsparse(inits),3))
    }
    globals2 <- list(obj2 = obj2, mydll=mydll, rotation=rotation)
    ## the user must pass data objects
    if(isRTMB) globals2 <- c(globals2,globals)
    if(is.null(adapt_stan_metric)){
      adapt_stan_metric <- ifelse(metric %in% c('stan', 'unit'), TRUE, FALSE)
    }

    if(!adapt_stan_metric){
      # this turns off diagonal mass matrix adaptation and
      # disables windowed step size adaptation which is wasteful
      # if not doing M
      control$metric <- 'unit_e' # Stan metric, not TMB's
      control$adapt_window <- 0
      if(is.null(num_warmup)) num_warmup <- 150
    }
    if(is.null(num_warmup)) num_warmup <- num_samples
    message("Starting MCMC sampling...")
    if(cores>1) message("Preparing parallel workspace...")
    fit <- StanEstimators::stan_sample(fn=fsparse, par_inits=inits,
                       grad_fun=gsparse, num_samples=num_samples,
                       num_warmup=num_warmup, thin=thin,
                       globals = globals2, packages=packages,
                       adapt_delta=control$adapt_delta,
                       adapt_window=control$adapt_window,
                       adapt_init_buffer=control$adapt_init_buffer,
                       adapt_term_buffer=control$adapt_term_buffer,
                       metric=control$metric, #Stan metric!
                       max_treedepth=control$max_treedepth,
                       parallel_chains=cores, save_warmup=TRUE,
                       num_chains = chains, seed = seed,
                       refresh=refresh, ...)

    fit2 <- as.tmbfit(fit, parnames=inputs$parnames, mle=inputs$mle,
                      invf=rotation$finv, metric=rotation$metric,
                      model=model_name)
    fit2$time.Q <- inputs$time.Q; fit2$time.Qinv <- inputs$time.Qinv;
    fit2$time.opt <- inputs$time.opt
    fit2$inits <- yinits
    ## gradient timings to check for added overhead
    if(require(microbenchmark)){
      bench <- microbenchmark(obj2$gr(inits),
                              gsparse(inits),
                              times=500, unit='s')
      fit2$time.gr <- summary(bench)$median[1]
      fit2$time.gr2 <- summary(bench)$median[2]
    } else {
      warning("Package microbenchmark required to do accurate gradient timings, using system.time() instead")
      fit2$time.gr <-
        as.numeric(system.time(trash <- replicate(1000, obj2$gr(inits)))[3])
      fit2$time.gr2 <-
        as.numeric(system.time(trash <- replicate(1000, gsparse(inits)))[3])
    }
    cat('\n\n')
    if(print) print(fit2)
    return(invisible(fit2))
  }
