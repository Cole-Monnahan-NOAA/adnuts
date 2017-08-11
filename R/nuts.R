#' [BETA VERSION] Draw MCMC samples from a model posterior using the
#' No-U-Turn (NUTS) sampler with dual averaging.
#'
#' @details This function implements algorithm 6 of Hoffman and Gelman
#'   (2014), which includes adaptive step sizes (\code{eps}) via an
#'   algorithm called dual averaging. In theory neither the step length nor
#'   step size needs to be input by the user to obtain efficient sampling
#'   from the posterior.
#' @param iter The number of samples to return.
#' @param eps The length of the leapfrog steps. If a numeric value is
#'   passed, it will be used throughout the entire chain. A \code{NULL}
#'   value will initiate adaptation of \code{eps} using the dual averaging
#'   algorithm during the first \code{warmup} steps.
#' @param warmup An optional argument for how many iterations to adapt
#'   \code{eps} in the dual averaging algorithm. A value of \code{NULL}
#'   results in a default of \code{warmup=iter/2}.
#' @param adapt_delta The target acceptance rate for the dual averaging
#'   algorithm. Defaults to 80\%. NUTS does not include an accept/reject
#'   Metropolis step, so this rate can be understood as the
#'   "average acceptance probability that HMC would give to the position-momentum states explored during the final doubling iteration."
#' @param fn A function that returns the log of the posterior density.
#' @param gr A function that returns a vector of gradients of the log of
#'   the posterior density (same as \code{fn}).
#' @param covar An optional covariance matrix which can be used to improve
#'   the efficiency of sampling. The lower Cholesky decomposition of this
#'   matrix is used to transform the parameter space. If the posterior is
#'   approximately multivariate normal and \code{covar} approximates the
#'   covariance, then the transformed parameter space will be close to
#'   multivariate standard normal. In this case the algorithm will be more
#'   efficient, but there will be overhead in the matrix calculations which
#'   need to be done at each step. The default of NULL specifies to not do
#'   this transformation.
#' @param init A vector of initial parameter values.
#' @param max_treedepth Integer representing the maximum times the path
#'   length should double within an MCMC iteration. Default of 4, so 16
#'   steps. If a U-turn has not occured before this many steps the
#'   algorithm will stop and return a sample from the given tree.
#' @references \itemize{ \item{Neal, R. M. (2011). MCMC using Hamiltonian
#'   dynamics. Handbook of Markov Chain Monte Carlo.}  \item{Hoffman and
#'   Gelman (2014). The No-U-Turn sampler: Adaptively setting path lengths
#'   in Hamiltonian Monte Carlo. J. Mach. Learn. Res.  15:1593-1623.}  }
#' @return A list containing samples ('par') and algorithm details such as
#'   step size adaptation and acceptance probabilities per iteration
#'   ('sampler_params').
#' @seealso \code{sample_tmb} and \code{run_mcmc.rwm}
run_mcmc.nuts <- function(iter, fn, gr, init, warmup=floor(iter/2),
                          chain, thin=1, seed=NULL, control=NULL){
  ## Now contains all required NUTS arguments
  if(!is.null(seed)) set.seed(seed)
  control <- update_control(control)
  eps <- control$stepsize
  init <- as.vector(unlist(init))
  metric <- control$metric
  if(is.null(metric)) metric <- diag(length(init))
  if(is.matrix(metric)){
    M <- metric
  } else {
    stop("Only allowed to pass M matrix as metric")
  }
  max_td <- control$max_treedepth
  adapt_delta <- control$adapt_delta
  adapt_mass <- control$adapt_mass
  if(warmup < 100 & adapt_mass){
    warning("Mass matrix adaptation not allowed if warmup < 100")
    adapt_mass <- FALSE
  }
  ## Mass matrix adapatation algorithm arguments
  w1 <- 75; w2 <- 50; w3 <- 25
  aws <- w2 # adapt window size
  anw <- w1+w2 # adapt next window

  ## Using a mass matrix means redefining what fn and gr do and
  ## backtransforming the initial value.
  if(!is.null(M)){
    temp <- rotate_space(fn=fn, gr=gr, M=M, y.cur=init)
    fn2 <- temp$fn2; gr2 <- temp$gr2
    theta.cur <- temp$theta.cur
    chd <- temp$chd
  } else {
    fn2 <- fn; gr2 <- gr
    theta.cur <- init
  }
  npar <- length(theta.cur)
  sampler_params <- matrix(numeric(0), nrow=iter, ncol=6,
                           dimnames=list(NULL, c("accept_stat__", "stepsize__", "treedepth__",
                                                 "n_leapfrog__", "divergent__", "energy__")))
  ## This holds the rotated but untransformed variables ("y" space)
  theta.out <- matrix(NA, nrow=iter, ncol=npar)
  ## how many steps were taken at each iteration, useful for tuning
  j.results <- lp <- rep(NA, len=iter)
  useDA <- is.null(eps)               # whether to use DA algorithm
  if(useDA){
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
  message(paste('Starting NUTS at', time.start))
  for(m in 1:iter){
    ## Initialize this iteration from previous in case divergence at first
    ## treebuilding. If successful trajectory they are overwritten
    theta.minus <- theta.plus <- theta.cur
    theta.out[m,] <- if(is.null(M)) theta.cur else t(chd %*% theta.cur)
    lp[m] <- if(m==1) fn2(theta.cur) else lp[m-1]
    r.cur <- r.plus <- r.minus <-  rnorm(npar,0,1)
    H0 <- .calculate.H(theta=theta.cur, r=r.cur, fn=fn2)

    ## Draw a slice variable u
    u <- .sample.u(theta=theta.cur, r=r.cur, fn=fn2)
    j <- 0; n <- 1; s <- 1; divergent <- 0
    ## Track steps and divergences; updated inside .buildtree
    info <- as.environment(list(n.calls=0, divergent=0))
    while(s==1) {
      v <- sample(x=c(1,-1), size=1)
      if(v==1){
        ## move in right direction
        res <- .buildtree(theta=theta.plus, r=r.plus, u=u, v=v,
                          j=j, eps=eps, H0=H0,
                          fn=fn2, gr=gr2, info=info)
        theta.plus <- res$theta.plus
        r.plus <- res$r.plus
      } else {
        ## move in left direction
        res <- .buildtree(theta=theta.minus, r=r.minus, u=u, v=v,
                          j=j, eps=eps, H0=H0,
                          fn=fn2, gr=gr2, info=info)
        theta.minus <- res$theta.minus
        r.minus <- res$r.minus
      }
      ## test whether to accept this state
      if(!is.finite(res$s)) res$s <- 0
      if(res$s==1) {
        if(runif(n=1, min=0,max=1) <= res$n/n){
          theta.cur <- res$theta.prime
          lp[m] <- fn2(theta.cur)
          ## Rotate parameters
          theta.out[m,] <- if(is.null(M)) theta.cur else t(chd %*% theta.cur)
        }
      }
      n <- n+res$n
      s <- as.vector(res$s*.test.nuts(theta.plus, theta.minus, r.plus, r.minus))
      if(!is.finite(s)) s <- 0
      j <- j+1
      ## Stop doubling if too many or it's diverged enough
      if(j>=max_td) {
        ## warning("j larger than max_treedepth, skipping to next m")
        break
      }
    }
    j.results[m] <- j-1

    alpha2 <- res$alpha/res$nalpha
    if(!is.finite(alpha2)) alpha2 <- 0
    ## ---------------
    ## Step size adapation with the
    ## Do the adapting of eps.
    if(useDA){
      if(m <= warmup){
        ## Adaptation during warmup:
        Hbar[m+1] <- (1-1/(m+t0))*Hbar[m] +
          (adapt_delta-alpha2)/(m+t0)
        ## If logalpha not defined, skip this updating step and use
        ## the last one.
        ## if(is.nan(Hbar[m+1])) Hbar[m+1] <- abs(Hbar[m])
        logeps <- mu-sqrt(m)*Hbar[m+1]/gamma
        epsvec[m+1] <- exp(logeps)
        logepsbar <- m^(-kappa)*logeps + (1-m^(-kappa))*log(epsbar[m])
        epsbar[m+1] <- exp(logepsbar)
        eps <- epsvec[m+1]
      } else {
        ## Fix eps for sampling period
        eps <- epsbar[warmup]
      }
    }
    ## ---------------
    ## Do the adaptation of mass matrix. The algorithm is working in X
    ## space but I need to calculate the mass matrix in Y space. So need to
    ## do this coversion in the calcs below.
    if(adapt_mass & slow_phase(m, warmup, w1, w3)){
      ## If in slow phase, update running estimate of variances
      ## The Welford running variance calculation, see
      ## https://www.johndcook.com/blog/standard_deviation/
      if(m== w1){
        ## Initialize algorithm from end of first fast window
        m1 <- theta.out[m,]; s1 <- rep(0, len=npar); k <- 1
      } else if(m==anw){
        ## If at end of adaptation window, update the mass matrix to the estimated
        ## variances
        vars <- as.numeric(s1/(k-1)) # estimated variance
        M <- diag(x=vars)
        ## Update density and gradient functions for new mass matrix
        temp <- rotate_space(fn=fn, gr=gr, M=M,  y.cur=theta.out[m,])
        fn2 <- temp$fn2; gr2 <- temp$gr2; chd <- temp$chd;
        theta.cur <- temp$theta.cur
        ## Reset the running variance calculation
        k <- 1; m1 <- theta.out[m,]; s1 <- rep(0, len=npar)
        ## Calculate the next end window. If this overlaps into the final fast
        ## period, it will be stretched to that point (warmup-w3)
        anw <- compute_next_window(m, anw, warmup, w1, aws, w3)
      } else {
        k <- k+1; m0 <- m1; s0 <- s1
        ## Update M and S
        m1 <- m0+(theta.out[m,]-m0)/k
        s1 <- s0+(theta.out[m,]-m0)*(theta.out[m,]-m1)
      }
    }
    ## End of mass matrix adaptation
    ##---------------

    ## Save adaptation info.
    sampler_params[m,] <-
      c(alpha2, eps, j, info$n.calls, info$divergent, fn2(theta.cur))
    if(m==warmup) time.warmup <- difftime(Sys.time(), time.start, units='secs')
    .print.mcmc.progress(m, iter, warmup, chain)
  } ## end of MCMC loop

  ## Process the output for returning
  theta.out <- cbind(theta.out, lp)
  theta.out <- theta.out[seq(1, nrow(theta.out), by=thin),]
  sampler_params <- sampler_params[seq(1, nrow(sampler_params), by=thin),]
  ndiv <- sum(sampler_params[-(1:warmup),5])
  if(ndiv>0)
    message(paste0("There were ", ndiv, " divergent transitions after warmup"))
  msg <- paste0("Final acceptance ratio=", sprintf("%.2f", mean(sampler_params[-(1:warmup),1])))
  if(useDA) msg <- paste0(msg,", and target=", adapt_delta)
  message(msg)
  if(useDA) message(paste0("Final step size=", round(eps, 3),
                           "; after ", warmup, " warmup iterations"))
  time.total <- difftime(Sys.time(), time.start, units='secs')
  .print.mcmc.timing(time.warmup=time.warmup, time.total=time.total)
  return(list(par=theta.out, sampler_params=sampler_params,
              time.total=time.total, time.warmup=time.warmup,
              warmup=warmup/thin, max_treedepth=max_td))
}

#' Draw a slice sample for given position and momentum variables
.sample.u <- function(theta, r, fn)
  runif(n=1, min=0, max=exp(.calculate.H(theta=theta,r=r, fn=fn)))
#' Calculate the log joint density (Hamiltonian) value for given position and
#' momentum variables.
#' @details This function currently assumes iid standard normal momentum
#' variables.
.calculate.H <- function(theta, r, fn) fn(theta)-(1/2)*sum(r^2)
#' Test whether a "U-turn" has occured in a branch of the binary tree
#' created by \ref\code{.buildtree} function. Returns TRUE if no U-turn,
#' FALSE if one occurred
.test.nuts <- function(theta.plus, theta.minus, r.plus, r.minus){
  theta.temp <- theta.plus-theta.minus
  res <- (crossprod(theta.temp,r.minus) >= 0) *
    (crossprod(theta.temp, r.plus) >= 0)
  return(res)
}

#' A recursive function that builds a leapfrog trajectory using a balanced
#' binary tree.
#'
#' @references This is from the No-U-Turn sampler with dual averaging
#' (algorithm 6) of Hoffman and Gelman (2014).
#'
#' @details The function repeatedly doubles (in a random direction) until
#' either a U-turn occurs or the trajectory becomes unstable. This is the
#' 'efficient' version that samples uniformly from the path without storing
#' it. Thus the function returns a single proposed value and not the whole
#' trajectory.
#'
.buildtree <- function(theta, r, u, v, j, eps, H0, fn, gr,
                       delta.max=1000, info = environment() ){
  if(j==0){
    ## ## Useful code for debugging. Returns entire path to global env.
    ## if(!exists('theta.trajectory'))
    ##   theta.trajectory <<- data.frame(step=0, t(theta))
    ## base case, take one step in direction v
    r <- r+(v*eps/2)*gr(theta)
    theta <- theta+(v*eps)*r
    r <- r+(v*eps/2)*gr(theta)
    ## verify valid trajectory. Divergences occur if H is NaN, or drifts
    ## too from from true H.
    H <- .calculate.H(theta=theta, r=r, fn=fn)
    n <- log(u) <= H
    s <- log(u) < delta.max + H
    if(!is.finite(H) | s == 0){
     info$divergent <- 1; s <- 0
    }
    ## Acceptance ratio in log space: (Hnew-Hold)
    logalpha <- H-H0
    alpha <- min(exp(logalpha),1)
    info$n.calls <- info$n.calls + 1
    ## theta.trajectory <<-
    ##   rbind(theta.trajectory, data.frame(step=tail(theta.trajectory$step,1),t(theta)))
    return(list(theta.minus=theta, theta.plus=theta, theta.prime=theta, r.minus=r,
                r.plus=r, s=s, n=n, alpha=alpha, nalpha=1))
  } else {
    ## recursion - build left and right subtrees
    xx <- .buildtree(theta=theta, r=r, u=u, v=v, j=j-1, eps=eps,
                       H0=H0, fn=fn, gr=gr, info=info)
    theta.minus <- xx$theta.minus
    theta.plus <- xx$theta.plus
    theta.prime <- xx$theta.prime
    r.minus <- xx$r.minus
    r.plus <- xx$r.plus
    alpha <- xx$alpha
    nalpha <- xx$nalpha
    s <- xx$s
    if(!is.finite(s)) s <- 0
    nprime <- xx$n
    ## If it didn't fail, update the above quantities
    if(s==1){
      if(v== -1){
        yy <- .buildtree(theta=theta.minus, r=r.minus, u=u, v=v,
                         j=j-1, eps=eps, H0=H0,
                         fn=fn, gr=gr, info=info)
        theta.minus <- yy$theta.minus
        r.minus <- yy$r.minus
      } else {
        yy <- .buildtree(theta=theta.plus, r=r.plus, u=u, v=v,
                         j=j-1, eps=eps, H0=H0,
                         fn=fn, gr=gr, info=info)
        theta.plus <- yy$theta.plus
        r.plus <- yy$r.plus
      }
      ### Update elements:
      ## If both slice variables fail you get 0/0.
      nprime <- yy$n+ xx$n
      if(!is.finite(nprime)) {nprime <- 0}
      ## choose whether to keep this theta
      if(nprime>0)
        if(runif(1) <= yy$n/nprime)
          theta.prime <- yy$theta.prime
      alpha <- xx$alpha+yy$alpha
      nalpha <- xx$nalpha+yy$nalpha
      ## check for valid proposal
      b <- .test.nuts(theta.plus=theta.plus,
                      theta.minus=theta.minus, r.plus=r.plus,
                      r.minus=r.minus)
      s <- yy$s*b
    }
    return(list(theta.minus=theta.minus, theta.plus=theta.plus,
                theta.prime=theta.prime,
                r.minus=r.minus, r.plus=r.plus, s=s, n=nprime,
                alpha=alpha, nalpha=nalpha))
  }
}

#' Estimate a reasonable starting value for epsilon (step size) for a given
#' model, for use with Hamiltonian MCMC algorithms.
                                        #
#' This is Algorithm 4 from Hoffman and Gelman (2010) and is used in the
#' dual-averaging algorithms for both HMC and NUTS to find a reasonable
#' starting value.
#' @title Estimate step size for Hamiltonian MCMC algorithms
#' @param theta An initial parameter vector.
#' @param fn A function returning the log-likelihood (not the negative of
#' it) for a given parameter vector.
#' @param gr A function returning the gradient of the log-likelihood of a
#' model.
#' @param eps A value for espilon to initiate the algorithm. Defaults to
#' 1. If this is far too big the algorithm won't work well and an
#' alternative value can be used.
#' @return Returns the "reasonable" espilon invisible, while printing how
#' many steps to reach it.
#' @details The algorithm uses a while loop and will break after 50
#' iterations.
#'
.find.epsilon <- function(theta,  fn, gr, eps=1, verbose=TRUE){
  r <- rnorm(n=length(theta), mean=0, sd=1)
  ## Do one leapfrog step
  r.new <- r+(eps/2)*gr(theta)
  theta.new <- theta+eps*r.new
  r.new <- r.new+(eps/2)*gr(theta.new)
  H1 <- .calculate.H(theta=theta, r=r, fn=fn)
  H2 <- .calculate.H(theta=theta.new, r=r.new, fn=fn)
  a <- 2*(exp(H2)/exp(H1)>.5)-1
  ## If jumped into bad region, a can be NaN so setup algorithm to keep
  ## halving eps instead of throwing error
  if(!is.finite(a)) a <- -1
  k <- 1
  ## Similarly, keep going if there are infinite values
  while (!is.finite(H1) | !is.finite(H2) | a*H2-a*H1 > -a*log(2)) {
    eps <- (2^a)*eps
    ## Do one leapfrog step
    r.new <- r+(eps/2)*gr(theta)
    theta.new <- theta+eps*r.new
    r.new <- r.new+(eps/2)*gr(theta.new)
    H2 <- .calculate.H(theta=theta.new, r=r.new, fn=fn)
    k <- k+1
    if(k>50) {
      stop("More than 50 iterations to find reasonable eps. Model is likely misspecified or some other issue.")
    }
  }
  if(verbose) message(paste("Reasonable epsilon=", eps, "found after", k, "steps"))
  return(invisible(eps))
}
