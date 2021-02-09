
## #' Run a single NUTS chain for an ADMB model
## #'
## #' A low level function to run a single chain. Unlikely to be used by a
## #' user, instead prefer \code{\link{sample_nuts}}
## #' @inheritParams wrappers
## #' @param chain Chain number, for printing purposes only.
## #' @param admb_args Character string of extra command line argument to
## #' pass to ADMB.
## #' @param extra.args Deprecated, use \code{admb_args} instead
## #' @param verbose Boolean for whether to print ADMB output to console.
## #' @seealso \code{\link{sample_nuts}}
sample_admb_nuts <- function(path, model, iter=2000,
                             init=NULL, chain=1,
                             thin=1, warmup=NULL,
                             seed=NULL, duration=NULL,
                             control=NULL,
                             skip_optimization=TRUE,
                             verbose=TRUE, admb_args=admb_args){

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
  model2 <- .update_model(model)
  if(skip_optimization){
    cmd <- paste(model2,"-nox -nohess -maxfn 0 -phase 1000 -nuts -mcmc ",iter)
  } else {
    cmd <- paste(model2,"-hbf -nuts -mcmc ",iter)
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
    stop(paste0("NUTS failed to run. Command attempted was:\n", cmd))
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


## #' Run a single random walk Metropolis chain for an ADMB model
## #'
## #' A low level function to run a single chain. Unlikely to be used by a
## #' user, instead prefer \code{\link{sample_rwm}}
## #' @inheritParams wrappers
## #' @seealso \code{\link{sample_rwm}}
sample_admb_rwm <- function(path, model, iter=2000, thin=1, warmup=ceiling(iter/2),
                            init=NULL,  chain=1, seed=NULL, control=NULL,
                            verbose=TRUE, duration=NULL,
                            admb_args=NULL, skip_optimization=TRUE){

  wd.old <- getwd(); on.exit(setwd(wd.old))
  setwd(path)
  ## Only refresh is used by RWM
  if(any(names(control) !='refresh'))
    warning("Only refresh control argument is used with RWM, ignoring: ",
            paste(names(control)[names(control)!='refresh'],
                  collapse=', '), call.=FALSE)
  refresh <- control$refresh
  if(!is.null(refresh) & !is.numeric(refresh))
    stop("Invalid refresh value ", refresh)
  metric <- 'mle' ## only one allowed
  stopifnot(iter >= 1)
  stopifnot(warmup <= iter)
  stopifnot(duration > 0)
  stopifnot(thin >=1)
  if(is.null(warmup)) stop("Must provide warmup")
  if(thin < 1 | thin > iter) stop("Thin must be >1 and < iter")


  ## Build the command to run the model
  if(skip_optimization){
    cmd <- paste("-nox -nohess -maxfn 0 -phase 1000 -rwm -mcmc ",iter)
  } else {
    cmd <- paste("-rwm -mcmc ",iter)
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
  if(!is.null(refresh)) cmd <- paste(cmd, "-refresh", refresh)
  if(!is.null(admb_args)) cmd <- paste(cmd, admb_args)


  ## Run it and get results
  model2 <- .update_model(model)
  console <- adnuts:::.check_console_printing()
  progress <- NULL
  if(console){
    ## Normal case
    time <- system.time(system2(model2, cmd, stdout=''))[3]
  } else {
    ## RStudio won't print output so capture it and print at
    ## end. Better than nothing
    fn <- 'mcmc_progress.txt'
    if(file.exists(fn)) file.remove(fn)
    time <- system.time(system2(model2, cmd, stdout=fn))[3]
    if(file.exists(fn)){
      progress <- readLines('mcmc_progress.txt')
      cat(progress, sep='\n')
    } else {
      warning("Progress output file not found. Try troubleshooting in serial model")
    }
  }
  if(!file.exists('unbounded.csv'))
    stop(paste0("RWM failed to run. Command attempted was:\n", cmd))
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
              par.names=par.names, cmd=cmd, unbounded=unbounded,
              progress=progress))
}



