#' Draw MCMC samples from a model posterior using a static
#' HMC sampler.
#'
#' @details This function implements algorithm 5 of Hoffman and Gelman
#'   (2014), which includes adaptive step sizes (\code{eps}) via an
#'   algorithm called dual averaging.
#' @inheritParams sample_tmb_nuts
#' @param L The number of leapfrog steps to take. The NUTS algorithm does
#'   not require this as an input. If \code{L=1} this function will perform
#'   Langevin sampling. In some contexts \code{L} can roughly be thought of
#'   as a thinning rate.
#' @param eps The step size. If a numeric value is passed, it will be used
#'   throughout the entire chain. A \code{NULL} value will initiate
#'   sampler_params of \code{eps} using the dual averaging algorithm during
#'   the first \code{warmup} steps.
#' @references \itemize{ \item{Neal, R. M. (2011). MCMC using Hamiltonian
#'   dynamics. Handbook of Markov Chain Monte Carlo.}  \item{Hoffman and
#'   Gelman (2014). The No-U-Turn sampler: Adaptively setting path lengths
#'   in Hamiltonian Monte Carlo. J. Mach. Learn. Res.  15:1593-1623.}  }
#' @seealso \code{\link{sample_tmb}}
#' @return A list containing samples ('par') and algorithm details such as
#'   step size adaptation and acceptance probabilities per iteration
#'   ('sampler_params').
sample_tmb_hmc <- function(iter, fn, gr, init, L, eps,
                          warmup=floor(iter/2), seed=NULL,
                         chain=1,thin=1, control=NULL){
  warning("NUTS should be prefered to sHMC except in rare, specific cases")
  if(!is.null(seed)) set.seed(seed)
  control <- .update_control(control)
  metric <- control$metric
  adapt_delta <- control$adapt_delta
  ## If using covariance matrix and Cholesky decomposition, redefine
  ## these functions to include this transformation. The algorithm will
  ## work in the transformed space
  if(!is.null(metric)){
    fn2 <- function(theta) fn(chd %*% theta)
    gr2 <- function(theta) as.vector( t( gr(chd %*% theta) ) %*% chd )
    chd <- t(chol(metric))               # lower triangular Cholesky decomp.
    chd.inv <- solve(chd)               # inverse
    theta.cur <- chd.inv %*% init
  } else {
    fn2 <- fn; gr2 <- gr
    theta.cur <- init
  }
  accepted <- divergence <- lp <- rep(NA, iter)
  theta.out <- matrix(NA, nrow=iter, ncol=length(init))
  sampler_params <- matrix(numeric(0), nrow=iter, ncol=4, # holds DA info by iteration
                       dimnames=list(NULL, c("accept_stat__",
                                                "stepsize__", "int_time__", "energy__")))
  ## A NULL value for eps signifies to use the dual averaging algorithm
  useDA <- is.null(eps)
  if(useDA){
    ## Initialize the dual-averaging algorithm.
    epsvec <- Hbar <- epsbar <- rep(NA, length=warmup+1)
    eps <- epsvec[1] <- epsbar[1] <-
      .find.epsilon(theta=theta.cur, fn=fn2, gr=gr2, eps=.1, verbose=FALSE)
    mu <- log(10*eps)
    Hbar[1] <- 0; gamma <- 0.05; t0 <- 10; kappa <- 0.75
  } else {
    ## dummy values to return
    epsvec <- epsbar <- Hbar <- NULL
  }
  ## Start of MCMC chain
  time.start <- Sys.time()
  message('')
  message(paste('Starting static HMC at', time.start))
  for(m in 1:iter){
    ## Jitter step size to mitigate potential negative autocorrelations,
    ## only once fixed though
    if(useDA & m > warmup) eps <- eps*runif(1,.9,1.1)
    r.cur <- r.new <- rnorm(length(init),0,1)
    theta.new <- theta.cur
    theta.leapfrog <- matrix(NA, nrow=L, ncol=length(theta.cur))
    r.leapfrog <- matrix(NA, nrow=L, ncol=length(theta.cur))
    ## Make a half step for first iteration
    r.new <- r.new+eps*gr2(theta.new)/2
    for(i in 1:L){
      theta.leapfrog[i,] <- theta.new
      r.leapfrog[i,] <- r.new
      theta.new <- theta.new+eps*r.new
      ## Full step except at end
      if(i!=L) r.new <- r.new+eps*gr2(theta.new)
      ## If divergence, stop trajectory earlier to save computation
      if(any(!is.finite(r.new)) | any(!is.finite(theta.new))) break
    }
    ## half step for momentum at the end
    r.new <- r.new+eps*gr2(theta.new)/2
    logalpha <- -fn2(theta.cur)+fn2(theta.new)+ sum(r.cur^2)/2-sum(r.new^2)/2
    ## Numerical divergence is registered as a NaN above. In this case we
    ## want to reject the proposal, mark the divergence, and adjust the
    ## step size down if still adapting (see below).
    if(!is.finite(logalpha)){
      divergence[m] <- 1
      logalpha <- -Inf
    } else {
      divergence[m] <- 0
    }
    ## Test whether to accept or reject proposed state
    if(log(runif(1)) < logalpha){
      ## accept the proposed state
      theta.cur <- theta.new
      accepted[m] <- TRUE
    } else {
      ## otherwise reject it and stay there
      accepted[m] <- FALSE
    }
    theta.out[m,] <- theta.cur
    lp[m] <- fn(theta.cur)
    if(useDA){
      ## Do the adapting of eps.
      if(m <= warmup){
        Hbar[m+1] <-
          (1-1/(m+t0))*Hbar[m] + (adapt_delta-min(1,exp(logalpha)))/(m+t0)
        logeps <- mu-sqrt(m)*Hbar[m+1]/gamma
        epsvec[m+1] <- exp(logeps)
        logepsbar <- m^(-kappa)*logeps + (1-m^(-kappa))*log(epsbar[m])
        epsbar[m+1] <- exp(logepsbar)
        eps <- epsvec[m+1]
      } else {
        eps <- epsbar[warmup]
      }
    }
    ## Save adaptation info.
    sampler_params[m,] <- c(min(1,exp(logalpha)), eps, eps*L, fn2(theta.cur))
    if(m==warmup) time.warmup <- difftime(Sys.time(), time.start, units='secs')
    .print.mcmc.progress(m, iter, warmup, chain)
  } ## end of MCMC loop
  ## Back transform parameters if metric is used
  if(!is.null(metric)) {
    theta.out <- t(apply(theta.out, 1, function(x) chd %*% x))
  }
  theta.out <- cbind(theta.out, lp)
  theta.out <- theta.out[seq(1, nrow(theta.out), by=thin),]
  sampler_params <- sampler_params[seq(1, nrow(sampler_params), by=thin),]
  if(sum(divergence[-(1:warmup)])>0)
    message(paste0("There were ", sum(divergence[-(1:warmup)]),
                   " divergent transitions after warmup"))
  message(paste0("Final acceptance ratio=", sprintf("%.2f", mean(accepted[-(1:warmup)])),
                 " and target is ", adapt_delta))
  if(useDA) message(paste0("Final step size=", round(epsbar[warmup], 3),
                           "; after ", warmup, " warmup iterations"))
  time.total <- difftime(Sys.time(), time.start, units='secs')
  .print.mcmc.timing(time.warmup=time.warmup, time.total=time.total)
  return(list(par=theta.out, sampler_params=sampler_params,
              time.total=time.total, time.warmup=time.warmup, warmup=warmup/thin))
}
