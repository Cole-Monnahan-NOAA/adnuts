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
                              control=NULL, seed=NULL){
  iter <- iter-warmup
  obj$env$beSilent()
  sdr <- sdreport(obj, getJointPrecision=TRUE)
  Q <- sdr$jointPrecision
  Qinv <- solve(Q)
  init <- obj$env$last.par.best
  ## Make parameter names unique if vectors exist
  parnames <- names(init)
  parnames <- as.vector((unlist(sapply(unique(parnames), function(x){
    temp <- parnames[parnames==x]
    if(length(temp)>1) paste0(temp,'[',1:length(temp),']') else temp
  }))))

  par <- obj$env$last.par.best
  mle <- list(nopar=length(par), est=par, se=sqrt(diag(Qinv)),
              cor=cov2cor(solve(Qinv)), parnames=parnames)
  ## rebuild without random effects
  mydll <- unclass(getLoadedDLLs()[[obj$env$DLL]])$path
  obj2 <- MakeADFun(data=obj$env$data, parameters=obj$env$parList(),
                    map=obj$env$map,
                    random=NULL, silent=TRUE, DLL=obj$env$DLL)
  rotation <- .rotate_space(fn=obj2$fn, gr=obj2$gr, M=Q, y.cur=init)
  fsparse <- function(v) {dyn.load(mydll); -rotation$fn2(v)}
  gsparse <- function(v) -as.numeric(rotation$gr2(v))
  initssparse <- rotation$x.cur
  invsparse <- function(theta.cur){
    J <- rotation$J
    chd <- rotation$chd
    t(as.numeric(J%*%Matrix::solve(chd, Matrix::solve(chd, J%*%theta.cur, system="Lt"), system="Pt")))
  }
  globals <- list(obj = obj2, mydll=mydll, rotation=rotation)
  fit <- stan_sample(fn=fsparse, par_inits=initssparse,
                     grad_fun=gsparse, num_samples=iter,
                     num_warmup=warmup, globals = globals,
                     packages = c("TMB", "Matrix"),
                     adapt_delta=control$adapt_delta,
                     parallel_chains=cores, save_warmup=TRUE,
                     num_chains = chains, seed = seed)

  fit2 <- as.tmbfit(fit, mle=mle, invf=invsparse)
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
  x <- list(samples=post, sampler_params=spl, mle=mle,
            monitor=mon, model='test',
            max_treedepth=x@metadata$max_depth,
            warmup=as.numeric(x@metadata$num_warmup),
            ## iter=as.numeric(x@metadata$num_samples)+as.numeric(x@metadata$num_warmup),
            algorithm='NUTS')
  adfit(x)
}
