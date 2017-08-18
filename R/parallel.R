
#' Combine multiple fits as returned from \code{sample_tmb} or
#' \code{sample_admb} run as a single chain.
#'
#' @param fits A list of fits, each having a single chain
#' @return A merged fit across chains.
combine_fits <- function(fits){
  z <- list()
  test <- lapply(fits, function(x) x$samples)
  samples <- array(NA, dim=c(nrow(test[[1]]), length(test), ncol(test[[1]])))
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

#' A wrapper for running ADMB models in parallel
#' @export
sample_admb_parallel <- function(parallel_number, path, algorithm, ...){
  olddir <- getwd()
  on.exit(setwd(olddir))
  newdir <- paste0(file.path(getwd(),path),"_chain_",parallel_number)
  if(dir.exists(newdir)){
    unlink(newdir, TRUE)
    if(dir.exists(newdir)) stop(paste("Could not remove folder:", newdir))
  }
  dir.create(newdir)
  trash <- file.copy(from=list.files(path, full.names=TRUE), to=newdir)
  if(algorithm=="NUTS")
    fit <- adnuts:::sample_admb_nuts(path=newdir, chain=parallel_number, ...)
  if(algorithm=="RWM")
    fit <- adnuts:::sample_admb_rwm(path=newdir, chain=parallel_number, ...)
  unlink(newdir, TRUE)
  return(fit)
}

#' A wrapper for running TMB models in parallel
sample_tmb_parallel <-  function(parallel_number, obj, init, path,
                                 algorithm, lower, upper, seed, ...){
  ## Each node starts in a random work directory. Rebuild TMB model obj so
  ## can link it in each session.
  setwd(path)
  dyn.load(dynlib(obj$env$DLL))
  ## Use 'shape' attribute to obtain full length of 'map'ped parameters.
  map.index <- which(names(obj$env$parameters) %in% names(obj$env$map))
  new.par <- obj$env$parameters
  new.par[map.index] <- lapply(obj$env$parameters[map.index], function(x) attr(x, "shape"))
  obj <- MakeADFun(data=obj$env$data, parameters=new.par, random=obj$env$random,
                   map=obj$env$map, DLL=obj$env$DLL, silent=TRUE)
  obj$env$beSilent()
  ## Parameter constraints, if provided, require the fn and gr functions to
  ## be modified to account for differents in volume. There are four cases:
  ## no constraints, bounded below, bounded above, or both (box
  ## constraint).
  bounded <- !(is.null(lower) & is.null(upper))
  if(bounded){
    if(is.null(lower)) lower <- rep(-Inf, len=length(upper))
    if(is.null(upper)) upper <- rep(Inf, len=length(lower))
    cases <- .transform.cases(lower, upper)
    fn <- function(y){
      x <- .transform(y, lower, upper, cases)
      scales <- .transform.grad(y, lower, upper, cases)
      -obj$fn(x) + sum(log(scales))
    }
    gr <- function(y){
      x <- .transform(y, lower, upper, cases)
      scales <- .transform.grad(y, lower, upper, cases)
      scales2 <- .transform.grad2(y, lower, upper, cases)
      -as.vector(obj$gr(x))*scales + scales2
    }
    ## Don't need to adjust this b/c init is already backtransformed in
    ## sample_tmb.
    ## init <- .transform.inv(x=unlist(init), a=lower, b=upper, cases=cases)
  } else {
    fn <- function(x) -obj$fn(x)
    gr <- function(x) -as.vector(obj$gr(x))
  }
  if(algorithm=="NUTS")
    fit <- run_mcmc.nuts(chain=parallel_number, fn=fn, gr=gr,
                         init=init, seed=seed, ...)
  if(algorithm=="RWM")
    fit <- run_mcmc.rwm(chain=parallel_number, fn=fn, init=init,
                        seed=seed, ...)
  return(fit)
}
