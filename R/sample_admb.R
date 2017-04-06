#' @export
#'
#'
sample_admb <- function(model, iter, init, chains=1, warmup=NULL, seeds=NULL,
                        thin=1, dir=getwd(), mceval=FALSE, duration=NULL,
                        parallel=FALSE, cores=NULL, control=NULL,
                        algorithm="NUTS", ...){
  ## Argument checking and processing
  control <- adnuts:::update_control(control)
  if(is.null(warmup)) warmup <- floor(iter/2)
  if(!(algorithm %in% c('NUTS', 'RWM')))
    stop("Invalid algorithm specified")
  if(is.null(init)){
    warning('Using MLE inits for each chain -- strongly recommended to use dispersed inits')
  } else if(length(init) != chains){
    stop("Length of init does not equal number of chains.")
  }
  if(is.null(seeds)) seeds <- sample.int(1e7, size=chains)
  stopifnot(thin >=1) ;stopifnot(chains >= 1)
  if(iter < 10 | !is.numeric(iter)) stop("iter must be > 10")
  ## Delete any psv files in case something goes wrong we dont use old
  ## values by accident
  trash <- file.remove(list.files()[grep('.psv', x=list.files())])
  ## Run in serial
  if(!parallel){
    if(algorithm=="NUTS"){
      mcmc.out <- lapply(1:chains, function(i)
        sample_admb_nuts(dir=dir, model=model, warmup=warmup, duration=duration,
                         iter=iter, init=init[[i]], chain=i,
                         seed=seeds[i], thin=thin, control=control, ...))
    } else {
      mcmc.out <- lapply(1:chains, function(i)
        sample_admb_rwm(dir=dir, model=model, warmup=warmup, duration=duration,
                        iter=iter, init=init[[i]], chain=i,
                        seed=seeds[i], thin=thin, control=control, ...))
    }
    ## Parallel execution
  } else {
    mcmc.out <- sfLapply(1:chains, function(i)
      sample_admb_parallel(parallel_number=i, dir=dir, model=model,
                           duration=duration,
                           algorithm=algorithm,
                           iter=iter, init=init[[i]], warmup=warmup,
                           seed=seeds[i], thin=thin, control=control, ...))
  }
  warmup <- mcmc.out[[1]]$warmup
  par.names <- mcmc.out[[1]]$par.names
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
  message("... Writing samples from all chains to psv file")
  samples2 <- do.call(rbind, lapply(1:chains, function(i)
    samples[, i, -dim(samples)[3]]))
  write_psv(fn=model, samples=samples2, model.path=dir)
  ## These already exclude warmup
  rotated <- do.call(rbind, lapply(mcmc.out, function(x) x$rotated))
  bounded <- do.call(rbind, lapply(mcmc.out, function(x) x$bounded))
  unbounded <- do.call(rbind, lapply(mcmc.out, function(x) x$unbounded))
  oldwd <- getwd(); on.exit(setwd(oldwd))
  setwd(dir)
  write.table(unbounded, file='unbounded.csv', sep=",", col.names=FALSE, row.names=FALSE)
  write.table(rotated, file='rotated.csv', sep=",", col.names=FALSE, row.names=FALSE)
  write.table(bounded, file='bounded.csv', sep=",", col.names=FALSE, row.names=FALSE)
  if(mceval)
    system(paste(model, "-mceval -noest -nohess"), ignore.stdout=FALSE)
  covar.est <- cov(unbounded)
  result <- list(samples=samples, sampler_params=sampler_params,
                 time.warmup=time.warmup, time.total=time.total,
                 algorithm=algorithm, warmup=warmup,
                 model=model, max_treedepth=mcmc.out[[1]]$max_treedepth,
                 cmd=cmd, covar.est=covar.est)
  return(invisible(result))
}


#' Run an MCMC using an ADMB model, return (1) the posterior
#' draws, MLE fits and covariance/correlation matrices, and some
#' MCMC convergence diagnostics using CODA.
#'
#' @param dir (Character) A path to the folder containing
#'   the model. NULL indicates the current folder.
#' @param mode.name (Character) The name of the model
#'   executable. A character string, without '.exe'.
#' @param iter (Integer) The number of draws after thinning and
#'   burn in.
#' @param thin (Integer) Controls thinning of samples. Save every
#'   thin value, such that 1 corresponds to keeping all draws,
#'   and 100 saving every 100th draw.
#' @param warmup (Integer) How many samples to discard from the
#'   beginning of the chain, *after* thining. The burn in period
#'   (i.e., the first warmup*thin draws) should be at least large
#'   enough to cover dynamic scaling.
#' @param covar (Numeric matrix) A manually defined covariance
#'   matrix (in bounded space) to use in the Metropolis-Hastings
#'   algorithm.
#' @param init (Numeric vector) A vector of initial values, which
#'   are written to file and used in the model via the -mcpin
#'   option.
#' @param seed (Integer) Which seed (integer value) to pass
#'   ADMB. Used for reproducibility.
#' @param mcdiag (Logical) Whether to use the \code{mcdiag}
#'   feature. This uses an identity matrix for the covariance
#'   matrix.
#' @param eps (Numeric) The size of the leapfrog jump in the
#'   hybrid method, with smaller values leading to smaller but
#'   more accurate jumps. Must be a positive value.
#' @param verbose (Logical) Whether to print ADMB warnings and
#'   other information. Useful for testing and troubleshooting.
#' @param extra.args (Character) A string which is passed to ADMB
#'   at runtime. Useful for passing additional arguments to the
#'   model executable.
#' @export
#' @return Returns a list containing (1) the posterior draws, (2)
#'   and object of class 'admb', read in using the results read
#'   in using \code{read_admb}, and (3) some MCMC convergence
#'   diagnostics using CODA.
sample_admb_nuts <- function(dir, model, iter, thin, warmup, duration=NULL,
                             init=NULL,  chain=1, par.names=NULL, seed=NULL, control=NULL,
                             verbose=TRUE, extra.args=NULL){
  wd.old <- getwd(); on.exit(setwd(wd.old))
  setwd(dir)
  ## Now contains all required NUTS arguments
  control <- adnuts:::update_control(control)
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

  ## Build the command to run the model
  cmd <- paste(model," -nohess -nuts -mcmc ",iter)
  cmd <- paste(cmd, "-warmup", warmup, "-chain", chain)
  if(!is.null(seed)) cmd <- paste(cmd, "-mcseed", seed)
  if(!is.null(duration)) cmd <- paste(cmd, "-duration", duration)
  cmd <- paste(cmd, "-max_treedepth", max_td, "-adapt_delta", adapt_delta)
  if(!is.null(eps)) cmd <- paste(cmd, "-hyeps", eps)

  ## Three options for metric. NULL (default) is to use the MLE estimates
  ## in admodel.cov. These need to be rescaled (see below), which means
  ## the model needs to be re-estimated. If a matrix is passed, this is
  ## written to file and no scaling is done. Option 'unit' means
  ## identity. Note: these are all in unbounded space.
  est <- FALSE
  if(is.matrix(metric)){
    ## User defined one will be writen to admodel.cov
    cor.user <- metric/ sqrt(diag(metric) %o% diag(metric))
    if(!matrixcalc:::is.positive.definite(x=cor.user))
      stop("Invalid mass matrix, not positive definite")
    write.admb.cov(metric, hbf=1)
  } else if(is.null(metric)) {
    ## MLE one. Need to re-estimate model to rescale covar
    est <- TRUE
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
  if(!est) cmd <- paste(cmd, "-noest ")
  if(!is.null(extra.args)) cmd <- paste(cmd, extra.args)

  ## Run it and get results
  time <- system.time(system(cmd, ignore.stdout=!verbose))[3]
  sampler_params<- as.matrix(read.csv("adaptation.csv"))
  unbounded <- as.matrix(read.csv("unbounded.csv", header=FALSE))
  rotated <- as.matrix(read.csv("rotated.csv", header=FALSE))
  bounded <- as.matrix(read.csv("bounded.csv", header=FALSE))
  dimnames(unbounded) <- dimnames(rotated) <- dimnames(bounded) <- NULL
  pars <- get_psv(model)
  if(is.null(par.names)) par.names <- names(pars)
  pars[,'log-posterior'] <- sampler_params[,'energy__']
  pars <- as.matrix(pars)
  ## Thin samples and adaptation post hoc for NUTS
  pars <- pars[seq(1, nrow(pars), by=thin),]
  bounded <- bounded[seq(1, nrow(bounded), by=thin),]
  unbounded <- unbounded[seq(1, nrow(unbounded), by=thin),]
  rotated <- rotated[seq(1, nrow(rotated), by=thin),]
  sampler_params <- sampler_params[seq(1, nrow(sampler_params), by=thin),]
  time.total <- time; time.warmup <- NA
  warmup <- warmup/thin
  return(list(samples=pars, sampler_params=sampler_params,
              time.total=time.total, time.warmup=time.warmup,
              warmup=warmup, max_treedepth=max_td,
              model=model, par.names=par.names, cmd=cmd,
              unbounded=unbounded, rotated=rotated, bounded=bounded))
}


get_psv <- function(model){
      if(!file.exists(paste0(model, '.psv'))){
      ## Sometimes ADMB will shorten the name of the psv file for some
      ## reason, so need to catch that here.
      ff <- list.files()[grep(x=list.files(), pattern='psv')]
      if(length(ff)==1){
        warning(paste("No .psv file found, using", ff))
        pars <- R2admb::read_psv(sub('.psv', '', x=ff))
      } else {
        stop(paste("No .psv file found -- did something go wrong??"))
      }
    } else {
      ## If model file exists
      pars <- R2admb::read_psv(model)
    }
  return(pars)
}


sample_admb_rwm <-
  function(dir, model, iter=2000, thin=1, warmup=ceiling(iter/2),
           init=NULL,  chain=1, seed=NULL, control=NULL, par.names=NULL,
           verbose=TRUE, extra.args=NULL, duration=NULL,
           mceval=TRUE){
    wd.old <- getwd(); on.exit(setwd(wd.old))
    setwd(dir)
    ## Now contains all required NUTS arguments
    control <- update_control(control)
    metric <- control$metric
    stopifnot(iter >= 1)
    stopifnot(warmup <= iter)
    stopifnot(duration > 0)
    stopifnot(thin >=1)
    if(is.null(warmup)) stop("Must provide warmup")
    if(thin < 1 | thin > iter) stop("Thin must be >1 and < iter")

    ## Build the command to run the model
    cmd <- paste(model," -nohess -mcmc ",iter)
    cmd <- paste(cmd, "-mcscale", warmup, "-chain", chain)
    if(!is.null(seed)) cmd <- paste(cmd, "-mcseed", seed)
    if(!is.null(duration)) cmd <- paste(cmd, "-duration", duration)
    cmd <- paste(cmd, "-nosdmcmc -mcsave", thin)

    ## Three options for metric. NULL (default) is to use the MLE estimates
    ## in admodel.cov. These need to be rescaled (see below), which means
    ## the model needs to be re-estimated. If a matrix is passed, this is
    ## written to file and no scaling is done. Option 'unit' means
    ## identity. Note: these are all in unbounded space.
    est <- FALSE
    if(is.matrix(metric)){
      ## User defined one will be writen to admodel.cov
      cor.user <- metric/ sqrt(diag(metric) %o% diag(metric))
      if(!matrixcalc:::is.positive.definite(x=cor.user))
        stop("Invalid mass matrix, not positive definite")
      write.admb.cov(metric)
    } else if(is.null(metric)) {
      ## MLE one. Should not need to re-estimate model to rescale covar
      est <- FALSE
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
    if(!est) cmd <- paste(cmd, "-noest ")
    if(!is.null(extra.args)) cmd <- paste(cmd, extra.args)

    ## Run it and get results
    time <- system.time(system(cmd, ignore.stdout=!verbose))[3]
    unbounded <- as.matrix(read.csv("unbounded.csv", header=FALSE))
    rotated <- as.matrix(read.csv("rotated.csv", header=FALSE))
    bounded <- as.matrix(read.csv("bounded.csv", header=FALSE))
    dimnames(unbounded) <- dimnames(rotated) <- dimnames(bounded) <- NULL
    pars <- get_psv(model)
    if(is.null(par.names))  par.names <- names(pars)
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
                rotated=rotated, bounded=bounded))
  }



