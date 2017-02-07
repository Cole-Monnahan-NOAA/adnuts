
#' Run an MCMC using an ADMB model, return (1) the posterior
#' draws, MLE fits and covariance/correlation matrices, and some
#' MCMC convergence diagnostics using CODA.
#'
#' @param model.path (Character) A path to the folder containing
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
run_admb_nuts <-
  function(model.path, model.name, iter=2000, thin=1, warmup=ceiling(iter/2),
           init=NULL, eps=NULL, covar=NULL,  mcseed=NULL,
           mcdiag=FALSE, verbose=TRUE, extra.args=NULL, max_treedepth=12,
           mceval=TRUE, chain=1){
    wd.old <- getwd(); on.exit(setwd(wd.old))
    setwd(model.path)
    ## Grab original admb fit and metrics
    if(iter <1)
      stop("Iterations must be >1")
    if(!file.exists(paste0(model.name,'.par'))) {
      system(paste("admb", model.name))
      system(paste(model.name))
    }
    ## If user provided covar matrix, write it to file and save to
    ## results
    mle <- R2admb::read_admb(model.name, verbose=TRUE)
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
    cmd <- model.name
    if(!est)
      cmd <- paste(cmd, " -noest -mcpin init.pin")
    cmd <- paste(cmd," -nohess -nuts -mcmc ",iter)
    cmd <- paste(cmd, "-chain", chain)
    if(!is.null(extra.args))
      cmd <- paste(cmd, extra.args)
    if(!is.null(mcseed))
      cmd <- paste(cmd, "-mcseed", mcseed)
    if(!is.null(eps)){
      cmd <- paste(cmd, "-hyeps", eps)
    }
    ## Run it and get results
    system(cmd, ignore.stdout=!verbose)
    if(mceval) system(paste(model.name, "-mceval -noest -nohess"),
                      ignore.stdout=!verbose)
    sampler_params<- as.matrix(read.csv("adaptation.csv"))
    pars <- read_psv(model.name)
    pars[,'log-posterior'] <- sampler_params[,'energy__']
    pars <- as.matrix(pars)
    time.total <- 5; time.warmup <- 3
    ## Thin
    pars <- pars[seq(1, nrow(pars), by=thin),]
    sampler_params <- sampler_params[seq(1, nrow(sampler_params), by=thin),]
    ndiv <- sum(sampler_params[-(1:warmup),5])
    par.names <- names(mle$coefficients[1:mle$npar])
    return(list(par=pars, sampler_params=sampler_params,
                time.total=time.total, time.warmup=time.warmup,
                warmup=warmup/thin, max_treedepth=max_treedepth,
                model=model.name, mle=mle, par.names=par.names))
  }

#' @export
run_admb_mcmc <- function(model.path, model.name, iter, chains=1, init=NULL,
                     seeds=NULL, covar=NULL, thin=1, ...){
  ## Argument checking
  if(is.null(init)){
    warning('Using MLE inits for each chain -- strongly recommended to use dispersed inits')
  } else if(length(init) != chains){
    stop("Length of init does not equal number of chains.")
  }
  if(is.null(seeds)) seeds <- sample.int(1e7, size=chains)
  thin <- floor(thin)
  stopifnot(thin >=1)
  thin.ind <- seq(1, iter, by=thin)
  stopifnot(chains >= 1)
  if(iter < 10 | !is.numeric(iter)) stop("iter must be > 10")
  mcmc.out <- lapply(1:chains, function(i)
    run_admb_nuts(model.path=model.path, model.name=model.name,
                  iter=iter, init=init[[i]], covar=covar,
                  chain=i,  mcseed=seeds[i], ...))
  par.names=mcmc.out[[1]]$par.names
  samples <-  array(NA, dim=c(nrow(mcmc.out[[1]]$par), chains, 1+length(par.names)),
                    dimnames=list(NULL, NULL, c(par.names,'lp__')))
  for(i in 1:chains){samples[,i,] <- mcmc.out[[i]]$par}
  ## Drop=FALSE prevents it from dropping 2nd dimension when chains=1
  samples <- samples[thin.ind,,, drop=FALSE]
  sampler_params <- lapply(mcmc.out,
       function(x) x$sampler_params[thin.ind,])
  time.warmup <- unlist(lapply(mcmc.out, function(x) as.numeric(x$time.warmup)))
  time.total <- unlist(lapply(mcmc.out, function(x) as.numeric(x$time.total)))
  result <- list(samples=samples, sampler_params=sampler_params,
                 time.warmup=time.warmup, time.total=time.total,
                 algorithm="NUTS", warmup=mcmc.out[[1]]$warmup/thin,
                 model=mcmc.out[[1]]$model,
                 max_treedepth=mcmc.out[[1]]$max_treedepth)
  return(invisible(result))
}

#' @export
as.shinystan.admb <- function(admb.fit){
  sso <- with(admb.fit,
    as.shinystan(samples, warmup=warmup, max_treedepth=max_treedepth,
                 sampler_params=sampler_params, algorithm='NUTS',
                 model_name=model))
  ## pars2 <- array(0, dim=c(nrow(pars), 1, ncol(pars)))
  ## pars2[,1,] <- as.matrix(pars)
  ## dimnames(pars2) <-
  ##   list(iter=1:nrow(pars), chains="chain:1",
  ##        parameters=dimnames(pars)[[2]])
  ## ss <- monitor(sims=pars2)
  ## y <- vector("list", length=length(dimnames(pars2)[[3]]))
  ## names(y) <- dimnames(pars2)[[3]]
  ## z <- lapply(y, function(x) x=numeric(0))
  ## sso <-
  ##   shinystan:::shinystan(
  ##   model_name=model.name, param_names=names(pars), param_dims=z,
  ##   posterior_sample=pars2, sampler_params=list(adapt),
  ##   summary=ss, n_chain=1, n_iter=nrow(pars),
  ##   n_warmup=nrow(pars)/2, model_code='NA',
  ##   misc=list(max_td=12, stan_method='sampling',
  ##             stan_algorithm='NUTS',
  ##             sso_version=utils::packageVersion('shinystan')))

    return(sso)
}

#' @export
launch_shinystan_admb <- function(admb.fit){
  launch_shinystan(as.shinystan.admb(admb.fit))
}

