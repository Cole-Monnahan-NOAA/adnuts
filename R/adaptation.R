
#' Compute the next window size
#'
#' @param i MCMC iteration number
#' @param warmup Number of warmup iterations
#' @param w1 The first adapation window (usually 75)
#' @param w3 The last adaptation window (usually 50)
#' @details This function calculates the size of the next window for
#'   adapation. If the next window size would be too long then this is
#'   extended to the end of that window.
compute_next_window <- function(i, anw, warmup, w1, aws, w3){
  ##  if(anw == warmup-w3) stop("Something bad")
  aws <- aws
  anw <- i+aws
  if(anw== (warmup-w3) ) return(anw)
  ## Check that the next anw is not too long. This will be the anw for the
  ## next time this is computed.
  nwb <- anw+2*aws
  if(nwb >= warmup-w3){
    ## if(i != warmup-w3)
    ##   message(paste("Extending last slow window from", anw, "to", warmup-w3))
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
#' @param fn The current fn function.
#' @param gr The current gr function
#' @param y.cur The current parameter vector in unrotated (Y) space.
#' @param M The new mass matrix
rotate_space <- function(fn, gr, M,  y.cur){
  ## Rotation done using choleski decomposition
  chd <- t(chol(M))               # lower triangular Cholesky decomp.
  chd.inv <- solve(chd)               # inverse
  ## Redefine these functions
  fn2 <- function(theta) fn(chd %*% theta)
  gr2 <- function(theta) as.vector( t( gr(chd %*% theta) ) %*% chd )
  ## Need to adjust the current parameters so the chain is
  ## continuous. First rotate to be in Y space.
  ## Now rotate back to "x" space using the new mass matrix M
  x.cur <- chd.inv %*% y.cur
  return(list(gr2=gr2, fn2=fn2, theta.cur=x.cur, chd=chd))
}

## ## The basic algorithm
## m <- m0 <- c(1,1)
## s <- s0 <- c(0,0)
## k <- 1
## warmup <- 20
## X <- matrix(NA, warmup, 2)
## results <- matrix(NA, warmup, 7)
## X[1,] <- m
## results[1,] <- c(1, m, m, s)
## for(i in 1:warmup){
##   X[i,] <- theta <- rnorm(2)
##   ## The Welford running variance calculation, see
##   ## https://www.johndcook.com/blog/standard_deviation/
##   m0 <- m; s0 <- s
##   ## Update M and S
##   m <- m0+(theta-m0)/k
##   s <- s0+(theta-m0)*(theta-m)
##   results[i,] <- c(k, theta, m, s)
##   if(i!=warmup) k <- k+1
## }
## sqrt(s/(k-1)) - apply(X,2, sd)


## ## Quick sketch of algorithm used to test.
## w1 <- 3
## w2 <- 5
## w3 <- 4
## aws <- w2
## anw <- w1+w2
## warmup <- 56
## sd <- c(1,.01)
## X <- matrix(NA, warmup, 2)
## X[1,] <- c(1.1,1.1)
## s1 <- m1 <- c(NA,NA)
## k <- 0
## phase <- 'no adapt'
## temp2 <- matrix(NA, nrow=warmup, ncol=6)
## temp2[1,] <- c(1, k,X[1,1], 0, 0, 0)
## for(i in 2:warmup){
##   X[i,] <- theta <- rnorm(2, sd=sd)
##   if(slow_phase(i,warmup, w1, w3)){
##     ## If in slow phase, update running estimate of variances
##     ## The Welford running variance calculation, see
##     ## https://www.johndcook.com/blog/standard_deviation/
##     ## Save last iteration
##     if(i== (w1)){
##       phase <- c(phase,"init")
##       ## Initialize algorithm from end of first fast window
##       m1 <- theta
##       s1 <- c(0,0)
##       k <- 1
##     } else if(i==anw){
##       phase <- c(phase,"update + reset")
##       ## If at end of adaptation window, update the mass matrix to the estimated
##       ## variances
##       vars <- s1/(k-1)
##       message(paste("Updating M:", i, vars[1]))
##       ## Reset the running variance calculation
##       k <- 1; m1 <- theta; s1 <- c(0,0)
##       ## Calculate the next end window. If this overlaps into the final fast
##       ## period, it will be stretched to that point (warmup-w3)
##       aws <- 2*aws
##       anw0 <- anw
##       anw <- compute_next_window(i, anw, warmup, w1, aws, w3)
##       message(c("window=",anw0, "-", anw))
##     } else {
##       phase <- c(phase,"adapt")
##       k <- k+1
##       m0 <- m1; s0 <- s1
##       ## Update M and S
##       m1 <- m0+(theta-m0)/k
##       s1 <- s0+(theta-m0)*(theta-m1)
##     }
##   } else {
##     phase <- c(phase, 'no adapt')
##   }
##   temp2[i,] <- c(i, k, theta[1], m1[1], s1[1], (s1/(k-1))[1])
## }
## phase

## temp2 <- as.data.frame(temp2)
## names(temp2) <- c('i', 'k', 'x', 'm', 's', 'var')
## xx <- 18:51
## apply(X[xx,], 2, var)
## vars


