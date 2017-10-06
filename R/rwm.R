
#' [Deprecated] Draw MCMC samples from a model posterior using a
#' Random Walk Metropolis (RWM) sampler.
#'
#' @details This algorithm does not yet contain adaptation of \code{alpha}
#' so some trial and error may be required for efficient sampling.
#' @param alpha The amount to scale the proposal, i.e,
#' Xnew=Xcur+alpha*Xproposed where Xproposed is generated from a mean-zero
#' multivariate normal. Varying \code{alpha} varies the acceptance rate.
#' @return A list containing samples and other metadata.
#' @inheritParams sample_tmb_nuts
#' @seealso \code{\link{sample_tmb}}
sample_tmb_rwm <- function(iter, fn, init, alpha=1, chain=1,
                         warmup=floor(iter/2), thin=1,
                         seed=NULL, control=NULL){
  if(!is.null(seed)) set.seed(seed)
  control <- .update_control(control)
  lp <- accepted <- rep(0, length=iter)
  init <- as.vector(unlist(init))
  n.params <- length(init)
  theta.out <- matrix(NA, nrow=iter, ncol=n.params)
  ## If using covariance matrix and Cholesky decomposition, redefine
  ## these functions to include this transformation. The algorithm will
  ## work in the transformed space.
  metric <- control$metric
  if(!is.null(metric)){
    fn2 <- function(theta) fn(chd %*% theta)
    chd <- t(chol(metric))               # lower triangular Cholesky decomp.
    chd.inv <- solve(chd)               # inverse
    theta.cur <- chd.inv %*% init
  } else {
    fn2 <- fn
    theta.cur <- init
  }
  fn.cur <- fn2(theta.cur)
  time.start <- Sys.time()
  message('')
  message(paste('Starting RWM at', time.start))
  for(m in 1:iter){
    ## generate proposal
    theta.new <- theta.cur + alpha*rnorm(n=n.params, mean=0, sd=1)
    fn.new <- fn2(theta.new)
    if(log(runif(1))< fn.new-fn.cur){
      ## accept
      accepted[m] <- 1
      theta.cur <- theta.out[m,] <- theta.new
      fn.cur <- fn.new
    } else {
      ## do not accept
      theta.out[m,] <- theta.cur
    }
    lp[m] <- fn.cur
    if(m==warmup) time.warmup <- difftime(Sys.time(), time.start, units='secs')
    .print.mcmc.progress(m, iter, warmup, chain)
  } ## end of MCMC loop

  ## Back transform parameters if metric is used
  if(!is.null(metric)) {
    theta.out <- t(apply(theta.out, 1, function(x) chd %*% x))
  }
  theta.out <- cbind(theta.out, lp)
  theta.out <- theta.out[seq(1, nrow(theta.out), by=thin),]
  message(paste0("Final acceptance ratio=", sprintf("%.2f", mean(accepted[-(1:warmup)]))))
  time.total <- difftime(Sys.time(), time.start, units='secs')
  .print.mcmc.timing(time.warmup=time.warmup, time.total=time.total)
  return(list(par=theta.out, sampler_params=NULL,
              time.total=time.total, time.warmup=time.warmup, warmup=warmup/thin))
}




