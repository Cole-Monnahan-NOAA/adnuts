#' @export
#'
#'
sample_admb <- function(dir, model, iter, init, chains=1, seeds=NULL,
                        control=NULL, ...){
  control <- update_control(control)
  algorithm <- control$algorithm
  ## Argument checking
  if(is.null(init)){
    warning('Using MLE inits for each chain -- strongly recommended to use dispersed inits')
  } else if(length(init) != chains){
    stop("Length of init does not equal number of chains.")
  }
  if(is.null(seeds)) seeds <- sample.int(1e7, size=chains)
  thin <- floor(control$thin)
  stopifnot(thin >=1)
  thin.ind <- seq(1, iter, by=thin)
  stopifnot(chains >= 1)
  if(iter < 10 | !is.numeric(iter)) stop("iter must be > 10")
  if(algorithm=="NUTS"){
    mcmc.out <- lapply(1:chains, function(i)
      sample_admb_nuts(dir=dir, model=model,
                       iter=iter, init=init[[i]], chain=i,
                       mcseed=seeds[i], control=control, ...))
  } else if(algorithm=="RWM"){
    mcmc.out <- lapply(1:chains, function(i)
      sample_admb_rwm(dir=dir, model=model,
                       iter=iter, init=init[[i]], chain=i,
                       mcseed=seeds[i], control=control, ...))
  } else { stop("Invalid algorithm specified")}

  par.names=mcmc.out[[1]]$par.names
  samples <-  array(NA, dim=c(nrow(mcmc.out[[1]]$samples), chains, 1+length(par.names)),
                    dimnames=list(NULL, NULL, c(par.names,'lp__')))
  for(i in 1:chains){samples[,i,] <- mcmc.out[[i]]$samples}
  ## Drop=FALSE prevents it from dropping 2nd dimension when chains=1
  samples <- samples[thin.ind,,, drop=FALSE]
  sampler_params <- lapply(mcmc.out,
                           function(x) x$sampler_params[thin.ind,])
  time.warmup <- unlist(lapply(mcmc.out, function(x) as.numeric(x$time.warmup)))
  time.total <- unlist(lapply(mcmc.out, function(x) as.numeric(x$time.total)))
  result <- list(samples=samples, sampler_params=sampler_params,
                 time.warmup=time.warmup, time.total=time.total,
                 algorithm=algorithm, warmup=mcmc.out[[1]]$warmup/thin,
                 model=mcmc.out[[1]]$model,
                 max_treedepth=mcmc.out[[1]]$max_treedepth)
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
#' @param mcseed (Integer) Which seed (integer value) to pass
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
sample_admb_nuts <-
  function(dir, model, iter=2000, thin=1, warmup=ceiling(iter/2),
           init=NULL,  chain=1, mcseed=NULL, control=NULL,
           verbose=TRUE, extra.args=NULL,
           mceval=TRUE){
    wd.old <- getwd(); on.exit(setwd(wd.old))
    setwd(dir)
    ## Now contains all required NUTS arguments
    control <- update_control(control)
    eps <- control$stepsize
    metric <- control$metric
    if(metric=='unit') covar <- NULL
    if(metric=='diag') covar <- NULL
    max_td <- control$max_treedepth
    adapt_delta <- control$adapt_delta

    ## Grab original admb fit and metrics
    if(iter <1)
      stop("Iterations must be >1")
    if(!file.exists(paste0(model,'.par'))) {
      system(paste("admb", model))
      system(paste(model))
    }
    ## If user provided covar matrix, write it to file and save to
    ## results
    mle <- R2admb::read_admb(model, verbose=TRUE)
    if(!is.null(covar)){
      cor.user <- covar/ sqrt(diag(covar) %o% diag(covar))
      if(!matrixcalc:::is.positive.definite(x=cor.user))
        stop("Invalid covar matrix, not positive definite")
      write.admb.cov(covar)
      mle$covar <- covar
    }
    ## Write the starting values to file. Always using a
    ## init file b/c need to use -nohess -noest so
    ## that the covar can be specified and not
    ## overwritten. HOwever, this feature then starts the
    ## mcmc chain from the initial values instead of the
    ## MLEs. So let the user specify the init values, or
    ## specify the MLEs manually
    est <- FALSE
    if(is.null(init)){
      init <- mle$coefficients[1:mle$npar]
    } else if(init[1]=='mle') {
      est <- TRUE
    }
    write.table(file="init.pin", x=init, row.names=F, col.names=F)
    ## Separate the options by algorithm, first doing the shared
    ## arguments
    cmd <- model
    if(!est)
      cmd <- paste(cmd, " -noest -mcpin init.pin")
    cmd <- paste(cmd," -nohess -nuts -mcmc ",iter)
    cmd <- paste(cmd, "-chain", chain)
    cmd <- paste(cmd, "-max_treedepth", max_td)
    if(!is.null(extra.args))
      cmd <- paste(cmd, extra.args)
    if(!is.null(mcseed))
      cmd <- paste(cmd, "-mcseed", mcseed)
    if(!is.null(eps)){
      cmd <- paste(cmd, "-hyeps", eps)
    }
    ## Run it and get results
    time <- system.time(system(cmd, ignore.stdout=!verbose))[3]
    if(mceval) system(paste(model, "-mceval -noest -nohess"),
                      ignore.stdout=!verbose)
    sampler_params<- as.matrix(read.csv("adaptation.csv"))
    pars <- read_psv(model)
    pars[,'log-posterior'] <- sampler_params[,'energy__']
    pars <- as.matrix(pars)
    time.total <- time; time.warmup <- NA
    ## Thin
    pars <- pars[seq(1, nrow(pars), by=thin),]
    sampler_params <- sampler_params[seq(1, nrow(sampler_params), by=thin),]
    ndiv <- sum(sampler_params[-(1:warmup),5])
    par.names <- names(mle$coefficients[1:mle$npar])
    return(list(samples=pars, sampler_params=sampler_params,
                time.total=time.total, time.warmup=time.warmup,
                warmup=warmup/thin, max_treedepth=max_td,
                model=model, mle=mle, par.names=par.names))
  }




sample_admb_rwm <-
  function(dir, model, iter=2000, thin=1, warmup=ceiling(iter/2),
           init=NULL,  chain=1, mcseed=NULL, control=NULL,
           verbose=TRUE, extra.args=NULL,
           mceval=TRUE){
    wd.old <- getwd(); on.exit(setwd(wd.old))
    setwd(dir) ## Now contains all required NUTS arguments
    control <- update_control(control)
    metric <- control$metric
    if(metric=='unit') covar <- NULL
    if(metric=='diag') covar <- NULL
    ## Grab original admb fit and metrics
    if(iter <1)
      stop("Iterations must be >1")
    if(!file.exists(paste0(model,'.par'))) {
      system(paste("admb", model))
      system(paste(model))
    }
    ## If user provided covar matrix, write it to file and save to
    ## results
    mle <- R2admb::read_admb(model, verbose=TRUE)
    if(!is.null(covar)){
      cor.user <- covar/ sqrt(diag(covar) %o% diag(covar))
      if(!matrixcalc:::is.positive.definite(x=cor.user))
        stop("Invalid covar matrix, not positive definite")
      write.admb.cov(covar)
      mle$covar <- covar
    }
    ## Write the starting values to file. Always using a init file b/c need
    ## to use -nohess -noest so that the covar can be specified and not
    ## overwritten. HOwever, this feature then starts the mcmc chain from
    ## the initial values instead of the MLEs. So let the user specify the
    ## init values, or specify the MLEs manually
    est <- FALSE
    if(is.null(init)){
      init <- mle$coefficients[1:mle$npar]
    } else if(init[1]=='mle') {
      est <- TRUE
    }
    write.table(file="init.pin", x=init, row.names=F, col.names=F)
    ## Separate the options by algorithm, first doing the shared
    ## arguments
    cmd <- model
    if(!est)
      cmd <- paste(cmd, " -noest -mcpin init.pin")
    cmd <- paste(cmd," -nohess -mcmc ",iter,  "-mcball", warmup)
    cmd <- paste(cmd, "-nosdmcmc -mcsave", thin)
    if(!is.null(extra.args))
      cmd <- paste(cmd, extra.args)
    if(!is.null(mcseed))
      cmd <- paste(cmd, "-mcseed", mcseed)
    ## Run it and get results
    time <- system.time(system(cmd, ignore.stdout=!verbose))[3]
    if(mceval) system(paste(model, "-mceval -noest -nohess"),
                      ignore.stdout=!verbose)
    pars <- read_psv(model)
    pars[,'log-posterior'] <- pars[,1]
    warning("log posterior column in RWM is still broken!!")
    pars <- as.matrix(pars)
    time.total <- time; time.warmup <- NA
    par.names <- names(mle$coefficients[1:mle$npar])
    return(list(samples=pars,  time.total=time.total,
                time.warmup=time.warmup,
                warmup=warmup,  model=model, mle=mle,
                par.names=par.names))
  }
