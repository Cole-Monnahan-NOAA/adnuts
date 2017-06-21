
#' Check whether adaptation is in the slow phase
#'
#' @param i MCMC iteration number
#' @param warmup Number of warmup iterations
#' @param w1 The first adapation window (usually 75)
#' @param w3 The last adaptation window (usually 50)
compute_next_window <- function(i, anw, warmup, w1, aws, w3){
 ##  if(anw == warmup-w3) stop("Something bad")
  aws <- aws*2
  anw <- i+aws
  if(anw== (warmup-w3) ) return(anw)
  ## Check that the next anw is not too long. This will be the anw for the
  ## next time this is computed.
  nwb <- anw+2*aws
  if(nwb >= warmup-w3){
    message(paste("Extending last slow window from", anw, "to", warmup-w3))
    anw <- warmup-w3
  }
  return(anw)
}

#' Check whether adaptation is in the slow phase
#'
#' @param i MCMC iteration number
#' @param warmup Number of warmup iterations
#' @param w1 The first adapation window (usually 75)
#' @param w3 The last adaptation window (usually 50)
#' @return Bool whether in slow phase
#' @details During the slow phase the mass matrix is updated in a series of
#'   expanding windows. See Stan manual on adaptation.
slow_phase <- function(i, warmup, w1, w3){
  ## After w1, before start of w3
  x1 <- i>= w1 # after initial fast window
  x2 <- i<= (warmup-w3) # but before last fast window
  x3 <- i < warmup # definitely not during sampling
  return(x1 & x2 & x3)
}

#' Update algorithm for mass matrix.
#'
rotate_space <- function(fn, gr, M, theta.cur){
  ## Rotation done using choleski decomposition
  chd <- t(chol(M))               # lower triangular Cholesky decomp.
  chd.inv <- solve(chd)               # inverse
  ## Redefine these functions
  fn2 <- function(theta) fn(chd %*% theta)
  gr2 <- function(theta) as.vector( t( gr(chd %*% theta) ) %*% chd )
  ## Need to adjust the current parameters so the chain is continuous
  theta.cur <- chd.inv %*% theta.cur
  x <- list(gr2=gr2, fn2=fn2, theta.cur=theta.cur, chd=chd)
}



## ## Quick sketch of algorithm used to test.
## w1 <- 75
## w2 <- 50
## w3 <- 25
## aws <- w2
## anw <- w1+w2
## warmup <- 400
## ## Initialize algorithm
## M <- c(.5,.5)
## sd <- .1
## m <- m0 <- rnorm(2, sd=sd)
## s <- s0 <- c(0,0)
## i <- 1
## k <- 2
## X <- matrix(NA, warmup, 2)
## for(i in 1:warmup){
##   X[i,] <- theta <- rnorm(2, sd=sd)
##   ## If in slow phase, update running estimate of variances
##   if(slow_phase(i,warmup, w1, w3)){
##     ## The Welford running variance calculation, see
##     ## https://www.johndcook.com/blog/standard_deviation/
##     ## Save last iteration
##     m0 <- m; s0 <- s
##     ## Update M and S
##     m <- m0+(theta-m0)/k
##     s <- s0+(theta-m0)*(theta-m)
##     k <- k+1
##   }
##   ## If at end of adaptation window, update the mass matrix to the estimated
##   ## variances
##   if(i==anw & slow_phase(i,warmup, w1, w3)){
##     ## Update the mass matrix
##     M <- s/(k-1)
##     message(paste("Updating M:", i, M[1]))
##     ## temp <- rotate_space(fn, gr, M, theta.cur)
##     ## Reset the running variance calculation
##     k <- 2; m0 <- m <- theta; s0 <- s <- 0
##     ## Calculate the next end window. If this overlaps into the final fast
##     ## period, it will be stretched to that point (warmup-w3)
##     aws <- 2*aws
##     anw <- compute_next_window(i, anw, warmup, w1, aws, w3)
##   }
## }
## M
## apply(X[125:375,],2, var)
