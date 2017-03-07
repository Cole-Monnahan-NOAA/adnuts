
#' Combine multiple fits as returned from \code{sample_tmb} or
#' \code{sample_admb} run as a single chain.
#'
#' @param fits A list of fits, each having a single chain
#' @return A merged fit across chains.
combine_parallel_fits <- function(fits){
  z <- list()
  test <- lapply(fits, function(x) x$samples)
  samples <- array(NA, dim=c(nrow(test[[1]]), length(test), dim(test[[1]])[3]))
  dimnames(samples) <- dimnames(fits[[1]]$samples)
  for(i in 1:length(test)) samples[,i,] <- test[[i]]
  z$samples <- samples
  sp <- sapply(fits, function(x) x$sampler_params)
  z$sampler_params <- sp
  z$time.warmup <- unlist(lapply(fits, function(x) x$time.warmup))
  z$time.total <- unlist(lapply(fits, function(x) x$time.total))
  z$algorithm <- fits[[1]]$algorithm
  z$warmup <- fits[[1]]$warmup
  z$model <- fits[[1]]$model
  z$max_treedepth <- fits[[1]]$max_treedepth
  return(z)
}

#' A wrapper for running SS models in parallel
#' @export
sample_admb_parallel <- function(parallel_number, dir, ...){
  library(R2admb)
  olddir <- getwd()
  on.exit(setwd(olddir))
  newdir <- paste0(getwd(),'/model',parallel_number)
  if(dir.exists(newdir)) unlink(newdir, TRUE)
  dir.create(newdir)
  trash <- file.copy(from=list.files(dir, full.names=TRUE), to=newdir)
  ## delay in case indexing ties up files briefly
  Sys.sleep(5)
  fit <- sample_admb(chain=1, dir=newdir, ...)
  fit
}
