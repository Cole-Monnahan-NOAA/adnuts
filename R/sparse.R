
#' Deprecated version of sample_snuts
#' @export
sample_sparse_tmb <- function(...){
  .Deprecated('sample_snuts', package='adnuts')
  sample_snuts(...)
}


#' Extract posterior samples from a tmbfit object
#' @param x A fitted tmbfit object
#' @param invf The inverse function to decorrelate the parameters
#' @param parnames A vector of parameter names, excluding lp__
#' @param array Whether to return a data.frame (default) or array
#'   which is used in constructing other objects downstream
#' @export
get_post <- function(x, invf, parnames, array=FALSE) {
  p <- x@draws |> as.data.frame()
  q <- subset(p, select=-c(lp__, .iteration, .draw, .chain))
  names(q) <- parnames
  if(ncol(q)==1){
    q <- as.data.frame(apply(q, 1, invf)) |> cbind(p$lp__)
  } else {
    q <- as.data.frame(t(apply(q, 1, invf))) |> cbind(p$lp__)
  }
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
#' @param parnames A character vector of unique par names
#' @param mle A list of MLE parameters
#' @param invf The inverse function for the parameters
#' @param metric The metric used
#' @param model A character giving the model name
#' @export
as.tmbfit <- function(x, parnames, mle, invf, metric, model='anonymous'){
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
  thin <- as.numeric(x@metadata$thin)
  warmup <- ceiling(as.numeric(x@metadata$num_warmup)/thin)
  iter <- ceiling(as.numeric(x@metadata$num_samples)/thin)
  if(dim(post)[1] != warmup + iter){
    stop("Error in output dimensions: iter=",iter, "; warmup=", warmup,
         "; nrow(post)=", nrow(post))
  }
  x <- list(samples=post, sampler_params=spl, mle=mle,
            monitor=mon, model=model,
            metric=metric,
            par_names=parnames,
            max_treedepth=x@metadata$max_depth,
            warmup=warmup, iter=iter, thin=thin,
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
  if(NROW(x)==1) return(NULL) # not a matrix!
  nm <- deparse(substitute(x))
  e <- eigen(x,TRUE)
  mine <- min(e$value); maxe <- max(e$value); ratio <- maxe/mine
  pct.sparsity <- round(100*mean(x[lower.tri(x)] == 0),2)
  message(nm, " is ", pct.sparsity,
          "% zeroes, with condition factor=",round(ratio,0),
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
  if(metric=='stan'){
    parnames <- .make_unique_names(names(obj$env$last.par.best))
    # mle <- list(nopar=length(obj$env$last.par.best),
    #             parnames=parnames)
    out <- list(time.opt=time.opt,
                time.Qinv=time.Qinv, time.Q=time.Q,
                parnames=parnames,
                laplace=laplace, metric=metric)
    return(out)
  }
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
    mle <- obj$env$last.par.best
    ## Make parameter names unique if vectors exist
    parnames <- .make_unique_names(names(mle))
    if(is.null(Q) & hasRE){
      message("Getting Q...")
      time.Q <- as.numeric(system.time(
        sdr <- sdreport(obj, getJointPrecision=TRUE))[3])
      Q <- sdr$jointPrecision
    }
    if(!is.null(Q))  dimnames(Q) <- list(parnames, parnames)
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
    dimnames(Qinv) <- list(parnames,parnames)
    .print.mat.stats(Q)
    #.print.mat.stats(Qinv)
    stopifnot(all.equal(length(mle), nrow(Qinv)))
  } else { #laplace is turned on
    message("Getting M for fixed effects...")
    time.Qinv <- as.numeric(system.time(sdr <- sdreport(obj))[3])
    Qinv <- sdr$cov.fixed
    .print.mat.stats(Qinv)
    if(is.null(opt)){
      mle <- opt$par
    } else {
      mle <- obj$par
    }
    ## Make parameter names unique if vectors exist
    parnames <- .make_unique_names(names(mle))
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
  names(mle) <- parnames
  mle <- list(nopar=length(mle), est=mle, se=ses,
              cor=mycor, Q=Q)#,              Qinv=Qinv)
 out <- list(Q=Q, Qinv=Qinv, mle=mle, time.opt=time.opt,
             time.Qinv=time.Qinv, time.Q=time.Q, parnames=parnames,
             laplace=laplace, metric=metric)
 return(out)
}



#' Get a single initial value vector in untransformed model space
#' @param init The initial value strategy
#' @param obj2 The joint TMB model
#' @param seed RNG seed
#' @param inputs A list as returned by \code{.get_inputs}.
.get_inits <- function(init, metric, obj2, seed, inputs) {
  # only certain combinations of metrics and inputs can work
  metric <- inputs$metric
  if(metric=='stan' & init %in% c('random', 'random-t'))
    stop("'stan' metric not allowed with 'random' or 'random-t' init b/c no Qinv")
  if((is.null(inputs$Qinv | is.null(inputs$mle$est))) &
     init %in% c('random', 'random-t'))
    stop("Cannot use random inits b/c mode or Qinv does not exist, use a different init")
  if(!is.null(seed)) set.seed(seed)
  # this will be the mode if skip_optimization=FALSE in .get_inputs
  lpb <- obj2$env$last.par.best
  if(inputs$laplace){
    # joint vector needs to collapse down to just the fixed
    # effects when doing ELA
    lpb <- lpb[!obj2$env$lrandom()]
  }
  npars <- length(lpb)
  for(ii in 1:10){
    inits <-
      switch(init,
      'last.par.best' = lpb,
      'random-t'      = inputs$mle$est + mvtnorm::rmvt(n=1, sigma=inputs$Qinv, df=2),
      'random'        = inputs$mle$est + mvtnorm::rmvnorm(n=1, sigma=inputs$Qinv),
      'unif'          = runif(n=npars, min=-2, max=2))
    inits <- as.numeric(inits)
    if(length(inits)!=npars)
      stop("Wrong vector length for inits:", length(inits), " when should be", npars)
    # if(init=='last.par.best')
    #   inits <- obj2$env$last.par.best
    # if(init=='random-t')
    #   inits <- as.numeric(inputs$mle$est + mvtnorm::rmvt(n=1, sigma=inputs$Qinv, df=2))
    # if(init=='random')
    #   inits <- as.numeric(inputs$mle$est + mvtnorm::rmvnorm(n=1, sigma=inputs$Qinv))
    # if(init=='unif')
    #   inits <- as.numeric(runif(n=length(rotation$x.cur), min=-2, max=2))
    success <- is.finite(obj2$fn(inits)) & is.finite(sum(obj2$gr(inits)))
    if(success) break
  }
  if(!success)
    stop(init, " inits resulted in NaN log-posterior after 10 tries, try another method or investigate model")
  return(inits)
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
  if(metric=='dense'){
    ## First case is a dense mass matrix
    M <- as.matrix(Qinv)
    # took this out b/c it was warning too often, better way to test?
    #if(!matrixcalc::is.symmetric.matrix(M) ||
    #  !matrixcalc::is.positive.definite(M))
    # warning("Estimated dense matrix was not positive definite so may be unreliable. Try different metric or turn on the laplace if there are random effects if it fails.")
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
    M <- diag(as.matrix(Qinv))
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
  } else if(metric=='unit' | metric=='stan') {
    ## unit metric, change nothing
    fn2 <- function(x) fn(x)
    gr2 <- function(x) gr(x)
    x.cur <- y.cur
    finv <- function(x) x
    chd <- J <- NULL
  } else if(metric=='sparse-naive'){
    # This metric is carefully constructured to match the dense
    # metric up to numerical precision. But as it is slower it is
    # not typically used.
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
  } else if(metric=='sparse'){
    # Do Cholesky on Q permuted directly
    J <- NULL
    chd <- Matrix::Cholesky(Q, super=TRUE, perm=TRUE)
    L <- as(chd, "sparseMatrix")
    perm <- chd@perm + 1L
    iperm <- Matrix::invPerm(perm)
    # Drop all numerical zeros and convert to triangular storage
    L <- Matrix::tril(Matrix::drop0(L)) ## class(L) == "dtCMatrix"
    Lt <- Matrix::t(L) ## class(Lt) == "dtCMatrix"
    x.cur <- as.vector(Matrix::t(L) %*% y.cur[perm])
    fn2 <- function(y)  fn(Matrix::solve(Lt, y)[iperm])
    gr2 <- function(y){
      x <- Matrix::solve(Lt, y)[iperm]
      Matrix::solve(L, as.numeric(gr(x))[perm])
    }
    finv <- function(x)   as.numeric(Matrix::solve(Lt, x)[iperm])
  } else if(metric=='auto'){
    ## use recursion then pick the right one depending on several criteria
    if(!is.null(Q) && NROW(Q)>1) rsparse <- .rotate_posterior(metric='sparse', fn=fn, gr=gr, Q=Q, Qinv=Qinv, y.cur=y.cur)
    if(!is.null(Qinv))
      rdiag <- .rotate_posterior(metric='diag', fn=fn, gr=gr, Q=Q, Qinv=Qinv, y.cur=y.cur)
    if(!is.null(Qinv) && NROW(Qinv)>1){
      rdense <- tryCatch(.rotate_posterior(metric='dense', fn=fn, gr=gr, Q=Q, Qinv=Qinv, y.cur=y.cur),
                         error=function(e) "Failed")
    }
    runit <- .rotate_posterior(metric='unit', fn=fn, gr=gr, Q=Q, Qinv=Qinv, y.cur=y.cur)

    if(NROW(Qinv)==1){
      message("diag metric selected b/c only 1 parameter")
      return(rdiag)
    }

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
        if(NROW(cors)==1) {
          message("diag metric selected b/c only a single parameter")
          return(rdiag)
        } else if(max(abs(cors))<.8){
          message("diag metric selected b/c low correlations")
          return(rdiag)
        } else {
          message("dense metric selected b/c high correlation (",
                  round(max(abs(cors)),4), ")")
          return(rdense)
        }
      }
    } else {
      # has a Q
      cors <- cov2cor(Qinv)[lower.tri(Qinv, diag=FALSE)]
      if(max(abs(cors))<.8){
        message("diag metric selected b/c of low correlations")
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
            message("dense metric selected b/c faster than sparse and high correlation (",
                    round(max(abs(cors)),4), ")")
            return(rdense)
          } else {
            message("sparse metric selected b/c faster than sparse and high correlation (",
                    round(max(abs(cors)),4), ")")
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






#' Make an image plot showing the correlation (lower triangle)
#' and sparsity (upper triangle).
#'
#' @param fit A fitted object
#' @param Q A sparse matrix. If NULL it will be extracted from
#'   \code{fit}.
#'
#' @details This function is used to visualize the sparsity and
#'   correlation patterns of the joint model. The upper triangle
#'   shows whether an element is 0 (white) or not (gray), while
#'   the lower triangle shows the correlation calculated from
#'   \code{cov2cor(solve(Q))}.
#' @return A plot created by \code{Matrix::image}.
#' @export
plot_Q <- function(fit, Q=NULL){
  if(is.null(Q)){
    if(!is.adfit(fit)) stop("fit is not a valid fitted object")
    if(is.null(fit$mle$Q)) return(NULL)
    nn <- length(fit$par_names)
    if(is.null(fit$mle$Qinv)){
      corr <- cov2cor(as.matrix(Matrix::solve(fit$mle$Q)))
    } else {
      corr <- cov2cor(fit$mle$Qinv)
    }
    Q <- fit$mle$Q
  } else {
    corr <- cov2cor(solve(Q))
  }
  Q[Q!=0] <- 1e-10
  Q[lower.tri(Q,TRUE)] <- corr[lower.tri(Q,TRUE)]
  Matrix::image(Q, useRaster=TRUE, at=seq(-1,1, len=50))
}

