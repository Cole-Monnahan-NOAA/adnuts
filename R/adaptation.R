## These fucntions were modeled after those in Stan at:
## https://github.com/stan-dev/stan/blob/967f5717edc8a101c8712838d47e6d393ac4a3db/src/stan/mcmc/windowed_adaptation.hpp

## Helper functions for mass matrix adapatation.

## Compute the next window size
##
## @param i MCMC iteration number
## @param warmup Number of warmup iterations
## @param w1 The first adapation window (usually 75)
## @param aws The last adapt window size
## @param w3 The last adaptation window (usually 50)
## @return A list containing the end of the next adaptation window (anw)
## and the current window size (aws).
## @details This function calculates the size of the next window for
##   adapation. If the next window size would be too long then this is
##   extended to the end of that window.
.compute_next_window <- function(i, anw, warmup, w1, aws, w3){
  aws <- aws*2
  anw <- i+aws
  if(anw== (warmup-w3) ) return(anw)
  ## Check that the next anw is not too long. This will be the anw for the
  ## next time this is computed. If the next one is too long, extend this
  ## one to the very end.
  nwb <- anw+aws
  if(nwb >= warmup-w3){
    ## if(i != warmup-w3)
    ##   message(paste("Extending last slow window from", anw, "to", warmup-w3))
    anw <- warmup-w3
  }
  return(list(anw=anw, aws=aws))
}

## Check whether adaptation is in the slow phase
##
## @param i MCMC iteration number
## @param warmup Number of warmup iterations
## @param w1 The first adapation window (usually 75)
## @param w3 The last adaptation window (usually 50)
## @return Bool whether in slow phase
## @details During the slow phase the mass matrix is updated in a series of
##   expanding windows. See Stan manual on adaptation.
.slow_phase <- function(i, warmup, w1, w3){
  ## After w1, before start of w3
  x1 <- i>= w1 # after initial fast window
  x2 <- i<= (warmup-w3) # but before last fast window
  x3 <- i < warmup # definitely not during sampling
  return(x1 & x2 & x3)
}
