#' Fit a TMB model using a sparse inverse mass matrix
#' @param obj The TMB object with random effects turned on and
#'   optimized
#' @param iter Total iterations to run (warmup + sampling)
#' @param warmup Total warmup iterations
#' @param chains Number of chains
#' @param cores Number of parallel cores to use
#' @param control NUTS control list
#' @param seed Random number seed
#' @return A fitted MCMC object of class 'adfit'
#' @export
sample_sparse_tmb <- function(obj, iter, warmup, cores, chains,
                              control=NULL, seed=NULL,
                              init=c('last.par.best', 'random'),
                              metric=c('sparse','dense','diag', 'unit')){
  iter <- iter-warmup
  metric <- match.arg(metric)
  init <- match.arg(init)
  obj$env$beSilent()
  message("Getting Q...")
  time.Q <- as.numeric(system.time(sdr <- sdreport(obj, getJointPrecision=TRUE))[3])
  Q <- sdr$jointPrecision
  message("Inverting Q...")
  time.Qinv <- as.numeric(system.time(Qinv <- solve(Q))[3])
  init.mle <- obj$env$last.par.best
  if(metric=='dense') Q <- as.matrix(Qinv)
  if(metric=='diag') Q <- as.numeric(diag(Qinv))
  if(metric=='unit') Q <- rep(1, nrow(Q))
  ## Make parameter names unique if vectors exist
  parnames <- names(init.mle)
  parnames <- as.vector((unlist(sapply(unique(parnames), function(x){
    temp <- parnames[parnames==x]
    if(length(temp)>1) paste0(temp,'[',1:length(temp),']') else temp
  }))))

  par <- obj$env$last.par.best
  mle <- list(nopar=length(par), est=par, se=sqrt(diag(Qinv)),
              cor=cov2cor(solve(Qinv)), parnames=parnames, Q=Q)
  ## rebuild without random effects
  mydll <- unclass(getLoadedDLLs()[[obj$env$DLL]])$path
  isRTMB <- ifelse(obj$env$DLL=='RTMB', TRUE, FALSE)
  message("Rebuilding obj without random effects...")
  if(!isRTMB){
    obj2 <- TMB::MakeADFun(data=obj$env$data, parameters=obj$env$parList(),
                    map=obj$env$map,
                    random=NULL, silent=TRUE, DLL=obj$env$DLL)
  } else {
    if(is.null(obj$myfun))
      stop("Slot 'myfun' not found in RTMB obj. Please add it manually and retry")
    obj2 <- RTMB::MakeADFun(obj$myfun, parameters=obj$env$parList(),
                            map=obj$env$map,
                            random=NULL, silent=TRUE,
                            DLL=obj$env$DLL)
  }
  rotation <- .rotate_space(fn=obj2$fn, gr=obj2$gr, M=Q, y.cur=init.mle)
  fsparse <- function(v) {dyn.load(mydll); -rotation$fn2(v)}
  gsparse <- function(v) -as.numeric(rotation$gr2(v))
  inits <- rotation$x.cur
  if(init=='random'){
    if(!is.null(seed)) set.seed(seed)
    inits <- as.numeric(rotation$x.cur + mvtnorm::rmvt(n=1, sigma=diag(length(inits)), df=1))
  }
  finv <- rotation$finv
  globals <- list(obj2 = obj2, mydll=mydll, rotation=rotation)
  if(!isRTMB){
    packages = c("TMB", "Matrix")
  } else {
    packages = c("RTMB", "Matrix")
  }
  message("Starting MCMC sampling...")
  fit <- stan_sample(fn=fsparse, par_inits=inits,
                     grad_fun=gsparse, num_samples=iter,
                     num_warmup=warmup,
                     globals = globals, packages=packages,
                     adapt_delta=control$adapt_delta,
                     parallel_chains=cores, save_warmup=TRUE,
                     num_chains = chains, seed = seed)
  fit2 <- as.tmbfit(fit, mle=mle, invf=finv)
  fit2$time.Q <- time.Q; fit2$time.Qinv <- time.Qinv
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
            time.warmup=timing[1,], time.total=timing[1,]+timing[2,],
            ## iter=as.numeric(x@metadata$num_samples)+as.numeric(x@metadata$num_warmup),
            algorithm='NUTS')
  adfit(x)
}


