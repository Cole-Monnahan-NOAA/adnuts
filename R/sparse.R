#' NUTS sampling for TMB models using a sparse inverse mass matrix (beta)
#'
#' @param obj The TMB object with random effects turned on and
#'   optimized
#' @param iter Total iterations to run (warmup + sampling)
#' @param warmup Total warmup iterations. Defaults to
#'   \code{iter}/2 based on Stan defaults, but when using dense,
#'   sparse, or diag metrics a much shorter warmup can be used
#'   (e.g., 150), especially if paired with a 'unit_e' Stan
#'   metric. Use \code{\link{plot_sampler_params}} to investigate
#'   warmup performance and adjust as necessary for subsequent
#'   runs.
#' @param metric A character specifying which metric to use.
#'   Defaults to "auto" which uses an algorithm to select the
#'   best metric (see details), otherwise one of "sparse",
#'   "dense", "diag", or "unit" can be specified.
#' @param init Either 'last.par.best' (default) or 'random'. The
#'   former starts from the joint mode and the latter draws from
#'   a multivariate t distribution with df=1 centered at the mode
#'   using the inverse joint precision matrix as a covariance
#'   matrix. Note that StanEstimators only allows for the same
#'   init vector for all chains currently. If a seed is specified
#'   it will be set and thus the inits used will be reproducible.
#' @param chains Number of chains
#' @param cores Number of parallel cores to use, defaults to
#'   \code{chains} so set this to 1 to execute serial chains.
#' @param control NUTS control list, currently available options
#'   are 'adapt_delta', 'max_treedepth', and 'metric' which is
#'   the type of metric adaptation for Stan to do with options
#'   ('unit_e', 'diag_e', or 'dense_e'). For dense and sparse
#'   metrics this usually can be 'unit_e' to skip adaptation.
#'   NULL values (default) revert to \code{stan_sample} defaults.
#' @param seed Random number seed, used for generating inital
#'   values (if 'random") and for NUTS.
#' @param laplace Whether to leave the Laplace approximation on
#'   and only use NUTS to sample the fixed effects, or turn it
#'   off and sample from the joint parameter space (default).
#' @param skip_optimization Whether to skip optimization or not
#'   (default).
#' @param Q The sparse precision matrix. It is calculated internally if not specified (default).
#' @param Qinv The dense matrix (M). It is calculated internally if not specified (default).
#' @param globals A named list of objects to pass to new R
#'   sessions when running in parallel and using RTMB. Typically
#'   this is the `data` object for now.
#' @param model_name An optional character giving the model name.
#'   If NULL it will use the DLL name which for RTMB models is
#'   just 'RTMB'. The name is used only for printing.
#' @param refresh How often to print updates to console
#'   (integer). 0 will turn off printing. The default is 100.
#' @param ... Additional arguments to pass to
#'   \code{\link{StanEstimators::stan_sample}}.
#' @return A fitted MCMC object of class 'adfit'
#' @details The TMB metric is used to decorrelate/descale the
#'   posterior before sampling using the Stan algorithms via the
#'   StanEstimator interface. The default is 'auto' which uses an
#'   algorithm to determine the optimal metric to use for a
#'   particular model. The algorithm depends on whether Q and/or
#'   M=Qinv are available, the extent of parameter correlations,
#'   and the speed of gradient calculations. The chosen metric
#'   and reasoning are printed to the console before NUTS
#'   sampling. A sparse and dense matrix will decorrelate and
#'   descale the posterior in the same way (up to numerical
#'   precision), but the sparse one will be more efficient with
#'   high levels of sparsity and larger dimensions. The 'diag'
#'   option is to take the marginal SDs from M and thus only
#'   descales, while the 'unit' option is the default Stan
#'   algorithm and should be used with mass matrix adaptation.
#'   Note that the \code{metric} is the TMB metric and distinct
#'   from the Stan metric which is controlled via the
#'   \code{control} list.
#' @export
sample_sparse_tmb <-
  function(obj, iter=2000, warmup=floor(iter/2),
           chains=4, cores=chains,
           control=NULL, seed=NULL, laplace=FALSE,
           init=c('last.par.best', 'random'),
           metric=c('auto', 'sparse','dense','diag', 'unit'),
           skip_optimization=FALSE, Q=NULL, Qinv=NULL,
           globals=NULL, model_name=NULL, refresh=NULL,
           ...){

  iter <- iter-warmup
  metric <- match.arg(metric)
  init <- match.arg(init)
  obj$env$beSilent()
  if(!is.null(model_name)){
    stopifnot(is.character(model_name))
  } else {
    model_name <- obj$env$DLL
  }
  inputs <- .get_inputs(obj=obj, skip_optimization=skip_optimization,
                        laplace=laplace, metric=metric, Q=Q, Qinv=Qinv)
  mle <- inputs$mle
  Q <- inputs$Q
  Qinv <- inputs$Qinv

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
  rotation <- .rotate_posterior(metric=metric, fn=obj2$fn, gr=obj2$gr, Q=Q, Qinv=Qinv, y.cur=mle$est)
  metric <- rotation$metric # update if using auto for printing later
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
  if(cores>1) message("Preparing parallel workspace...")
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
                     num_chains = chains, seed = seed,
                     refresh=refresh, ...)

  fit2 <- as.tmbfit(fit, mle=mle, invf=finv, metric=metric, model=model_name)
  fit2$time.Q <- inputs$time.Q; fit2$time.Qinv <- inputs$time.Qinv; fit2$time.opt <- inputs$time.opt
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
#' @param metric The metric used
#' @param model A character giving the model name
#' @export
as.tmbfit <- function(x, mle, invf, metric, model='anonymous'){
  parnames <- mle$parnames
  ## move lp__ to end to match order of draws
  mon <- StanEstimators::summary(x)
  mon$variable <- c('lp__', parnames)
  mon <- rbind(mon[-1,], mon[1,])
  mon$n_eff <- mon$ess_bulk
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
            monitor=mon, model=model,
            metric=metric,
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

#' Prepare inputs for sparse sampling
#'  @param obj Object
#'  @param skip_optimization Whether to skip or not
#'  @param laplace Whether to due the LA or not
#'  @param metric Which metric
#'  @param Q Sparse precision
#'  @param Qinv Inverse of Q
#'  @return A list containing Q, Qinv, the mle list, and timings
.get_inputs <- function(obj, skip_optimization, laplace, metric, Q, Qinv) {

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
    if(metric %in% c('unit', 'auto')){
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
 out <- list(Q=Q, Qinv=Qinv, mle=mle, time.opt=time.opt, time.Qinv=time.Qinv, time.Q=time.Q)
 return(out)
}

#' Update algorithm for mass matrix.
#'
#' @param metric The metric to use
#' @param fn The current fn function.
#' @param gr The current gr function
#' @param y.cur The current parameter vector in unrotated (Y) space.
#' @param Q The sparse precision matrix
#' @param Qinv The inverse of Q
.rotate_posterior <- function(metric, fn, gr, Q,  Qinv, y.cur){
  ## Rotation done using choleski decomposition
  ## First case is a dense mass matrix
  M <- as.matrix(Qinv)
  if(metric=='dense'){
      if(!matrixcalc::is.symmetric.matrix(M) ||
         !matrixcalc::is.positive.definite(M))
        warning("Estimated dense matrix was not positive definite so may be unreliable. Try different metric or turn on the laplace if there are random effects if it fails.")
    J <- NULL
    chd <- t(chol(M))               # lower triangular Cholesky decomp.
    chd.inv <- solve(chd)               # inverse
    ## Define rotated fn and gr functions
    fn2 <- function(x) fn(chd %*% x)
    gr2 <- function(x) {as.vector( gr(chd %*% x) %*% chd )}
    ## Now rotate back to "x" space using the new mass matrix M
    x.cur <- as.numeric(chd.inv %*% y.cur)
    finv <- function(x){
      t(chd %*% x)
    }
  } else if(metric=='diag'){
    M <- diag(M)
    ## diagonal but not unit
      J <- NULL
      chd <- sqrt(M)
      fn2 <- function(x) fn(chd * x)
      gr2 <- function(x) as.vector(gr(chd * x) ) * chd
      ## Now rotate back to "x" space using the new mass matrix M. M is a
      ## vector here. Note the big difference in efficiency without the
      ## matrix operations.
      x.cur <- (1/chd) * y.cur
      finv <- function(x) chd*x
  } else if(metric=='unit') {
      ## unit metric, change nothing
      fn2 <- function(x) fn(x)
      gr2 <- function(x) gr(x)
      x.cur <- y.cur
      finv <- function(x) x
      chd <- J <- NULL
  } else if(metric=='sparse'){
    ##  warning( "Use of Q is highly experimental still" )
    stopifnot(require(Matrix))
    if(!is(Q,"Matrix")) stop("Q is not a Matrix object, something went wrong")
        # M is actually Q, i.e., the inverse-mass
    # Antidiagonal matrix JJ = I
    J = Matrix::sparseMatrix( i=1:nrow(Q), j=nrow(Q):1 )
    #chd <- Cholesky(M, super=FALSE, perm=FALSE)
    #chd <- Matrix::Cholesky(M, super=TRUE, perm=FALSE)
    chd <- Matrix::Cholesky(J%*%Q%*%J, super=TRUE, perm=FALSE) # perm
    Linv_times_x = function(chd,x){
      as.numeric(J%*% Matrix::solve(chd, Matrix::solve(chd, J%*%x, system="Lt"), system="Pt"))
    }
    x_times_Linv = function(chd,x){
      #x %*% chol()
      as.numeric(J%*%Matrix::solve(chd, Matrix::solve(chd, Matrix::t(x%*%J), system="L"), system="Pt"))
    }
    fn2 <- function(x){
      Linv_x = Linv_times_x(chd, x)
      fn(Linv_x)
    }
    gr2 <- function(x){
      Linv_x = Linv_times_x(chd, x)
      grad = gr( Linv_x )
      as.vector(  x_times_Linv(chd, grad) )
    }
    ## Now rotate back to "x" space using the new mass matrix M
    #  solve(t(chol(solve(M)))) ~~ IS EQUAL TO ~~ J%*%chol(M)%*%J
    # J%*%chol(J%*%prec%*%J) %*% J%*%x
    x.cur <- as.numeric(J%*%chol(J%*%Q%*%J) %*% J%*%y.cur)
    finv <- function(x){
      t(as.numeric(J%*%Matrix::solve(chd, Matrix::solve(chd, J%*%x, system="Lt"), system="Pt")))
    }
  } else if(metric=='auto'){
    ## use recursion then pick the right one depending on several criteria
    if(!is.null(Q)) rsparse <- .rotate_posterior(metric='sparse', fn=fn, gr=gr, Q=Q, Qinv=Qinv, y.cur=y.cur)
    if(!is.null(Qinv))
      rdiag <- .rotate_posterior(metric='diag', fn=fn, gr=gr, Q=Q, Qinv=Qinv, y.cur=y.cur)
    if(!is.null(Qinv)){
        rdense <- tryCatch(.rotate_posterior(metric='dense', fn=fn, gr=gr, Q=Q, Qinv=Qinv, y.cur=y.cur),
                           error=function(e) "Failed")
    }
    runit <- .rotate_posterior(metric='unit', fn=fn, gr=gr, Q=Q, Qinv=Qinv, y.cur=y.cur)

    if(is.character(rdense)){
      message("unit metric selected b/c Qinv was not positive definite")
      return(runit)
    }

    if(is.null(Q)){
      if(is.null(Qinv)){
        # must be unit since no other option
        message("unit metric selected b/c no Q or Qinv info available")
        return(runit)
      } else {
        # no Q but does have Qinv, e.g., a model w/o RE or using the LA
        ## check for high correlations
        cors <- cov2cor(Qinv)[lower.tri(Qinv, diag=FALSE)]
        if(max(abs(cors))<.3){
          message("diag metric selected b/c no Q available and low correlations")
          return(rdiag)
        } else {
          message("dense metric selected b/c no Q availabile and high correlations")
          return(rdense)
        }
      }
    } else {
      # has a Q
      cors <- cov2cor(Qinv)[lower.tri(Qinv, diag=FALSE)]
      if(max(abs(cors))<.3){
        message("diag metric selected b/c of and low correlations")
        return(rdiag)
      } else {
        if(!require(microbenchmark)){
          message("sparse metric selected b/c no timing available -- please install microbenchmark")
          ## check for speed differences
          return(rsparse)
        } else {
          bench <- microbenchmark::microbenchmark(rdense$gr2(rdense$x.cur),
                                                  rsparse$gr2(rsparse$x.cur),
                                                  times = 500)
          tdense <- summary(bench)$median[1]
          tsparse <- summary(bench)$median[2]
          if(tdense < tsparse){
            message("dense metric selected b/c faster than sparse and high correlations")
            return(rdense)
          } else {
            message("sparse metric selected b/c faster than dense and high correlations")
            return(rsparse)
          }
        }
      }
    }
  }  else {
    stop("Invalid metric")
  }
  ## Redefine these functions
  ## Need to adjust the current parameters so the chain is
  ## continuous. First rotate to be in Y space.
  return(list(gr2=gr2, fn2=fn2, finv=finv, x.cur=x.cur, chd=chd, J=J, metric=metric))
}


