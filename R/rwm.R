
#' [Deprecated] Draw MCMC samples from a model posterior using a
#' Random Walk Metropolis (RWM) sampler.
#'
#' @param iter The number of samples to return.
#' @param fn A function that returns the log of the posterior density.
#' @param init A vector of initial parameter values.
#' @param diagnostic Whether to return a list of diagnostic metrics about
#' the chain. Useful for assessing efficiency and tuning chain.
#' @details This algorithm does not yet contain adaptation of \code{alpha}
#' so some trial and error may be required for efficient sampling.
#' @param covar An optional covariance matrix which can be used to improve
#' the efficiency of sampling. The lower Cholesky decomposition of this
#' matrix is used to transform the parameter space. If the posterior is
#' approximately multivariate normal and \code{covar} approximates the
#' covariance, then the transformed parameter space will be close to
#' multivariate standard normal. In this case the algorithm will be more
#' efficient, but there will be overhead in the matrix calculations which
#' need to be done at each step. The default of NULL specifies to not do
#' this transformation.
#' @param alpha The amount to scale the proposal, i.e,
#' Xnew=Xcur+alpha*Xproposed where Xproposed is generated from a mean-zero
#' multivariate normal. Varying \code{alpha} varies the acceptance rate.
#' @return If \code{diagnostic} is FALSE (default), returns a matrix of
#' \code{iter} samples from the posterior. Otherwise returns a list
#' containing samples ('par'), proposed samples ('par.proposed'), vector of
#' which proposals were accepted ('accepted'), and the total function calls
#' ('n.calls'), which for this algorithm is \code{iter}
#' @seealso \code{\link{run_mcmc}}, \code{\link{run_mcmc.nuts}}, \code{\link{run_mcmc.hmc}}
sample_tmb_rwm <- function(iter, fn, init, alpha=1, chain=1,
                         warmup=floor(iter/2), thin=1,
                         seed=NULL, control){
  if(!is.null(seed)) set.seed(seed)
  control <- update_control(control)
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




