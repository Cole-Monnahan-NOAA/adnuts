#' Update the control list.
#'
#' @param control A list passed from \code{sample_tmb}.
#' @return A list with default control elements updated by those supplied
#'   in \code{control}
#' @export
update_control <- function(control){
  default <- list(adapt_delta=0.8, metric=NULL, stepsize=NULL,
                  algorithm="NUTS", adapt_engaged=TRUE,
                  max_treedepth=10)
  if(!is.null(control))
    for(i in names(control))  default[[i]] <- control[[i]]
  return(default)
}

#' Print MCMC progress to console.
#'
#' @param iteration The iteration of the MCMC chain.
#' @param iter The total iterations.
#' @param warmup The number of warmup iterations.
#' @param chain The chain being run (bookkeeping only).
#' @return Nothing. Prints to message to console.
#'
#' @details This function was modeled after the functionality provided by
# the R package \link{rstan}.
.print.mcmc.progress <- function(iteration, iter, warmup, chain){
  i <- iteration
  refresh <- max(10, floor(iter/10))
  if(i==1 | i==iter | i %% refresh ==0){
    i.width <- formatC(i, width=nchar(iter))
    out <- paste0('Chain ',chain,', Iteration: ', i.width , "/", iter, " [",
                  formatC(floor(100*(i/iter)), width=3), "%]",
                  ifelse(i <= warmup, " (Warmup)", " (Sampling)"))
    message(out)
  }
}

#' Print MCMC timing to console
#' @param time.warmup Time of warmup in seconds.
#' @param time.total Time of total in seconds.
#' @return Nothing. Prints message to console.
#'
#' @details This function was modeled after the functionality provided by
#'   the R package \link{rstan}.
.print.mcmc.timing <- function(time.warmup, time.total){
  x <- ' Elapsed Time: '
  message(paste0(x, sprintf("%.1f", time.warmup), ' seconds (Warmup)'))
  message(paste0(x, sprintf("%.1f", time.total-time.warmup), ' seconds (Sampling)'))
  message(paste0(x, sprintf("%.1f", time.total), ' seconds (Total)'))
}

#' Convert adnuts fit (named list) into a \code{shinystan} object.
#'
#' @details The shinystan packages provides several conversion functions
#'   for objects of different types, such as stanfit classes (Stan ouput)
#'   and simple arrays. For the latter, option NUTS information, such as
#'   \code{sampler_params} can be passed. This function essentially extends
#'   the functionality of \code{as.shinystan} to work specifically with
#'   fits from adnuts (TMB or ADMB). The user can thus explore their model
#'   with \code{launch_shinystan(as.shinystan.tmb(fit))} in the same way
#'   that Stan models are examined.
#' @param fit Output list from \link{\code{sample_tmb}} or
#'   \link{\code{sample_admb}}.
#' @seealso launch_shinytmb, launch_shinyadmb
#' @return An S4 object of class shinystan. Depending on the algorithm
#'   used, this list will have slight differences.
#' @export
as.shinyadnuts <- function(fit){
  if(fit$algorithm=="NUTS"){
    sso <- with(fit, as.shinystan(samples, warmup=warmup, max_treedepth=max_treedepth,
             sampler_params=sampler_params, algorithm='NUTS', model_name=model))
  } else if(fit$algorithm=="HMC"){
    sso <- with(fit, as.shinystan(samples, warmup=warmup,
             sampler_params=sampler_params, algorithm='HMC', model_name=model))
  } else {
    sso <- with(fit, as.shinystan(samples, warmup=warmup,
             algorithm='RWM', model_name=model))
  }
  return(sso)
}

#' A high level wrapper to launch shinystan for a TMB fit.
#'
#' @details This function simply calls
#'   \code{launch_shinystan(as.shinystan.tmb(tmb.fit))}.
#' @export
launch_shinytmb <- function(fit){
  launch_shinystan(as.shinyadnuts(fit))
}

#' Extract posterior samples from a TMB MCMC fit list.
#'
#' @param fit.tmb A list returned by \code{\link{run_mcmc}}.
#' @param inc_warmup Whether to extract the warmup samples or not
#'   (default). Warmup samples should never be used for inference, but may
#'   be useful for diagnostics.
#' @return An invisible data.frame containing samples (rows) of each
#'   parameter (columns). If multiple chains exist they will be rbinded
#'   together.
#' @export
extract_samples <- function(fit.tmb, inc_warmup=FALSE){
  x <- fit.tmb$samples
  if(!is.array(x)) stop("fit.tmb$samples is not an array -- valid TMB output?")
  ind <- if(inc_warmup) 1:dim(x)[1] else -(1:fit.tmb$warmup)
  y <- do.call(rbind, lapply(1:dim(x)[2], function(i) x[ind, i, -dim(x)[3]]))
  return(invisible(as.data.frame(y)))
}

#' A high level wrapper to launch shinystan for a ADMB fit.
#'
#' @details This function simply calls
#'   \code{launch_shinystan(as.shinystan.tmb(tmb.fit))}.
#' @export
launch_shinyadmb <- function(fit){
  launch_shinystan(as.shinyadnuts(fit))
}

#' Write matrix of samples to a binary .psv file.
#'
#' @details Useful to combine multiple MCMC runs together into a single
#' .psv file which can then be executed with '-mceval'.
#' @param fn Model name
#' @param samples A matrix or data.frame of samples, each column is a
#'   parameter, each row a sample.
write_psv <- function(fn, samples, model.path=getwd()){
  samples <- as.matrix(samples)
  psv <- file.path(getwd(), paste0(fn, '.psv'))
  con <- file(psv, 'wb')
  writeBin(object=ncol(samples), con=con)
  writeBin(object=as.vector(t(samples)), con=con)
  close(con)
}

#' Read in the ADMB covariance file.
#'
#' @export
get.admb.cov <- function(model.path=getwd()){
    wd.old <- getwd(); on.exit(setwd(wd.old))
    setwd(model.path)
    filename <- file("admodel.cov", "rb")
    on.exit(close(filename), add=TRUE)
    num.pars <- readBin(filename, "integer", 1)
    cov.vec <- readBin(filename, "numeric", num.pars^2)
    cov.unbounded <- matrix(cov.vec, ncol=num.pars, nrow=num.pars)
    hybrid_bounded_flag <- readBin(filename, "integer", 1)
    scale <- readBin(filename, "numeric", num.pars)
    cov.bounded <- cov.unbounded*(scale %o% scale)
    result <- list(num.pars=num.pars, cov.bounded=cov.bounded,
                   cov.unbounded=cov.unbounded,
                   hybrid_bounded_flag=hybrid_bounded_flag, scale=scale)
    return(result)
}

write.admb.cov <- function(cov.unbounded, model.path=getwd(), hbf=NULL){
  temp <- file.exists(paste0(model.path, "/admodel.cov"))
  if(!temp) stop(paste0("Couldn't find file ",model.path, "/admodel.cov"))
  temp <- file.copy(from=paste0(model.path, "/admodel.cov"),
                    to=paste0(model.path, "/admodel_original.cov"))
  wd.old <- getwd()
  setwd(model.path)
  ## Read in the output files
  results <- get.admb.cov()
  if(is.null(hbf)) hbf=results$hybrid_bounded_flag
  scale <- results$scale
  num.pars <- results$num.pars
  if(NROW(cov.unbounded) != num.pars)
    stop(paste0("Invalid size of covariance matrix, should be: ", num.pars,
                "instead of ",NROW(cov.unbounded), ". Do you need to rerun the model?"))
  ## Write it to file using original scales, although these are ignored.
  setwd(wd.old)
  file.new <- file(paste0(model.path, "/admodel.cov"),"wb")
  on.exit(close(file.new))
  writeBin(as.integer(num.pars), con=file.new)
  writeBin(as.vector(as.numeric(cov.unbounded)), con=file.new)
  writeBin(as.integer(hbf), con=file.new)
  writeBin(as.vector(scale), con=file.new)
}
