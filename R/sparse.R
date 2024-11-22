#' Fit a TMB model using a sparse inverse mass matrix
#' @param obj The TMB object with random effects turned on and
#'   optimized
#' @param iter Total iterations to run (warmup + sampling)
#' @param init Either 'last.par.best' (default) or 'random'. The
#'   former starts from the joint mode and the latter draws from
#'   a multivariate t distribution with df=1 centered at the mode
#'   using the inverse joint precision matrix as a covariance
#'   matrix. Note that StanEstimators only allows for the same
#'   init vector for all chains currently. If a seed is specified
#'   it will be set and thus the inits used will be reproducible.
#' @param warmup Total warmup iterations
#' @param chains Number of chains
#' @param cores Number of parallel cores to use
#' @param control NUTS control list, currently available options
#'   are 'adapt_delta', 'max_treedepth', and 'metric' which is
#'   the type of metric adaptation for Stan to do with options
#'   ('unit_e', 'diag_e', or 'dense_e'). For dense and sparse
#'   metrics this usually can be 'unit_e' to skip adaptation.
#'   NULL values (default) revert to \code{stan_sample} defaults.
#' @param seed Random number seed
#' @param laplace Whether to leave the Laplace approximation on
#'   and only use NUTS to sample the fixed effects, or turn it
#'   off and sample from the joint parameter space (default).
#' @param skip_optimization Whether to skip optimization or not
#'   (default).
#' @param globals A named list of objects to pass to new R
#'   sessions when running in parallel and using RTMB. Typically
#'   this is the `data` object for now.
#' @param ... Additional arguments to pass to
#'   \code{\link{StanEstimators::stan_sample}}.
#' @return A fitted MCMC object of class 'adfit'
#' @export
sample_sparse_tmb <- function(obj, iter, warmup, cores, chains,
                              control=NULL, seed=NULL, laplace=FALSE,
                              init=c('last.par.best', 'random'),
                              metric=c('sparse','dense','diag', 'unit'),
                              skip_optimization=FALSE, Q=NULL,
                              Qinv=NULL,
                              globals=NULL, ...){



  iter <- iter-warmup
  metric <- match.arg(metric)
  init <- match.arg(init)
  obj$env$beSilent()
  time.opt <- time.Q <- time.Qinv <- 0
  if(!skip_optimization){
    message("Optimizing...")
    time.opt <-
      as.numeric(system.time(opt <- with(obj, nlminb(par, fn, gr)))[3])
  }
  hasRE <-  !is.null(obj$env$random)
  if(laplace & !hasRE)
    stop("No random effects found so laplace=TRUE fails, set to FALSE")
  if( (laplace | !hasRE) & metric=='sparse')
    stop("sparse metric only allowed with random effects
         and laplace=FALSE")
  if(!laplace){
    if(is.null(Q) & hasRE){
      message("Getting Q...")
      time.Q <- as.numeric(system.time(
        sdr <- sdreport(obj, getJointPrecision=TRUE))[3])
      Q <- sdr$jointPrecision
    }
    if(is.null(Qinv)){
      if(!is.null(Q)){
        ## Q found above
        message("Inverting Q...")
        time.Qinv <- as.numeric(system.time(Qinv <- solve(Q))[3])
      } else if(!hasRE){
        ## fixed effect only model
        time.Qinv <- as.numeric(system.time(Qinv <- sdreport(obj)$cov.fixed)[3])
      } else {
       stop("something wrong here")
     }
    }
    .print.mat.stats(Q)
    .print.mat.stats(Qinv)
    mle <- obj$env$last.par.best
  ## Make parameter names unique if vectors exist
  parnames <- names(mle)
  parnames <- as.vector((unlist(sapply(unique(parnames), function(x){
    temp <- parnames[parnames==x]
    if(length(temp)>1) paste0(temp,'[',1:length(temp),']') else temp
  }))))
  stopifnot(all.equal(length(mle), nrow(Qinv)))
  } else {
      message("Getting M for fixed effects...")
    time.Qinv <- as.numeric(system.time(sdr <- sdreport(obj)))
    Qinv <- sdr$cov.fixed
    .print.mat.stats(Qinv)
    if(is.null(opt))
      stop("No opt object found, rerun with 'skip_optimization=FALSE'")
    mle <- opt$par
    ## Make parameter names unique if vectors exist
    parnames <- names(mle)
    parnames <- as.vector((unlist(sapply(unique(parnames), function(x){
      temp <- parnames[parnames==x]
      if(length(temp)>1) paste0(temp,'[',1:length(temp),']') else temp
    }))))
    stopifnot(all.equal(length(mle), nrow(Qinv)))
  }
  ses <- suppressWarnings(sqrt(diag(Qinv)))
  mycor <- suppressWarnings(cov2cor(Qinv))
  if(!all(is.finite(ses))){
    if(metric=='unit'){
    warning("Some standard errors estimated to be NaN, filling with dummy values so unit metric works. The 'mle' slot will be wrong so do not use it")
    cor <- diag(length(mle))
    ses <- rep(1,length(mle))
    } else {
      stop("Some standard errors estimated to be NaN, use 'unit' metric for models without a mode or positive definite Hessian")
    }
  }
  mle <- list(nopar=length(mle), est=mle, se=ses,
              cor=mycor, parnames=parnames, Q=Q,
              Qinv=Qinv)

  ## prepare to run via StanEstimators
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

  if(metric=='dense') Q <- as.matrix(Qinv)
  if(metric=='diag') Q <- as.numeric(diag(Qinv))
  if(metric=='unit') Q <- rep(1, length(mle$nopar))
  if(metric=='dense'){
    if(!matrixcalc::is.symmetric.matrix(Q) ||
       !matrixcalc::is.positive.definite(Q))
    warning("Estimated dense matrix was not positive definite so may be unreliable. Try different metric or turn on the laplace if there are random effects if it fails.")
  }
  rotation <- .rotate_space(fn=obj2$fn, gr=obj2$gr, M=Q, y.cur=mle$est)
  fsparse <- function(v) {dyn.load(mydll); -rotation$fn2(v)}
  gsparse <- function(v) -as.numeric(rotation$gr2(v))
  inits <- rotation$x.cur
  if(init=='random'){
    if(!is.null(seed)) set.seed(seed)
    inits <- as.numeric(rotation$x.cur + mvtnorm::rmvt(n=1, sigma=diag(length(inits)), df=1))
    if(!is.finite(obj2$fn(inits))) {
      warning("random inits resulted in NaN NLL, trying parameter mode instead")
      inits <- rotation$x.cur
    }
  }
  finv <- rotation$finv
  globals2 <- list(obj2 = obj2, mydll=mydll, rotation=rotation)
  ## the user must pass data objects
  if(isRTMB) globals2 <- c(globals2,globals)
  message("Starting MCMC sampling...")
  fit <- stan_sample(fn=fsparse, par_inits=inits,
                     grad_fun=gsparse, num_samples=iter,
                     num_warmup=warmup,
                     globals = globals2, packages=packages,
                     adapt_delta=control$adapt_delta,
                     adapt_window=control$adapt_window,
                     adapt_init_buffer=control$adapt_init_buffer,
                     adapt_term_buffer=control$adapt_term_buffer,
                     metric=control$metric,
                     max_treedepth=control$max_treedepth,
                     parallel_chains=cores, save_warmup=TRUE,
                     num_chains = chains, seed = seed, ...)
  fit2 <- as.tmbfit(fit, mle=mle, invf=finv)
  fit2$time.Q <- time.Q; fit2$time.Qinv <- time.Qinv; fit2$time.opt <- time.opt
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
  fit2$metric <- metric
  print(fit2)
  fit2
}


#' Extract posterior samples from a tmbfit object
#' @param x A fitted tmbfit object
#' @param invf The inverse function to decorrelate the parameters
#' @param parnames A vector of parameter names, excluding lp__
#' @param array Whether to return a data.frame (default) or array
#'   which is used in constructing other objects downstream
#' @export
get_post <- function(x, invf, parnames, array=FALSE) {
  p <- unconstrain_draws(x) |> as.data.frame()
  q <- subset(p, select=-c(lp__, .iteration, .draw, .chain))
  names(q) <- parnames
  q <- as.data.frame(t(apply(q, 1, invf))) |> cbind(p$lp__)
  colnames(q) <- c(parnames, 'lp__')
  ## build array
  if(array){
    samples <- array(NA, dim=c(max(p$.iter), max(p$.chain), 1 + length(parnames)),
                     dimnames = list(NULL, NULL, c(parnames, "lp__")))
    for(chain in 1:max(p$.chain))
      samples[,chain,] <- as.matrix(q[p$.chain==chain,])
    return(samples)
  }
  return(q)
}


#' Construtor for tmbfit objects
#' @param x A fitted MCMC object
#' @param mle A list of MLE parameters
#' @param invf The inverse function for the parameters
#' @export
as.tmbfit <- function(x, mle, invf){
  parnames <- mle$parnames
  ## move lp__ to end to match order of draws
  mon <- StanEstimators::summary(x)
  mon$variable <- c('lp__', parnames)
  mon <- rbind(mon[-1,], mon[1,])
  mon$n_eff <- mon$ess_tail
  mon$Rhat <- mon$rhat
  ## prepare objects to use the pairs_admb function
  post <- get_post(x, invf, parnames=parnames, TRUE)
  sp <- as.data.frame(x@diagnostics)
  spl <- list()
  for(chain in 1:max(sp$.chain)){
    spl[[chain]] <- as.matrix(sp[sp$.chain==chain,1:6])
  }
  timing <- sapply(x@timing, function(x) unlist(x))
  x <- list(samples=post, sampler_params=spl, mle=mle,
            monitor=mon, model='test',
            par_names=mle$parnames,
            max_treedepth=x@metadata$max_depth,
            warmup=as.numeric(x@metadata$num_warmup),
            time.warmup=timing[1,],
            time.sampling=timing[2,],
            time.total=timing[1,]+timing[2,],
            ## iter=as.numeric(x@metadata$num_samples)+as.numeric(x@metadata$num_warmup),
            algorithm='NUTS')
  adfit(x)
}

#' Print matrix stats
#'
#' @param x matrix object
#' @param name
#'
.print.mat.stats <- function(x){
  if(is.null(x)) return(NULL)
  nm <- deparse(substitute(x))
  e <- eigen(x,TRUE)
  mine <- min(e$value); maxe <- max(e$value); ratio <- maxe/mine
  ## why does this return NA without the as.numeric?
  pct.sparsity <- round(100*mean(as.numeric(x==0)),2)
  message(nm, " is ", pct.sparsity,
          " % zeroes, with condition factor=",round(ratio,0),
          ' (min=',round(mine,3), ', max=', round(maxe,1),")")
}
