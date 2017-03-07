#' Update the control list.
#'
#' @param control A list passed from \code{sample_tmb}.
#' @return A list with default control elements updated by those supplied
#'   in \code{control}
update_control <- function(control){
  default <- list(adapt_delta=0.8, metric='unit', stepsize=NULL,
                  algorithm="NUTS", adapt_engaged=TRUE, thin=1,
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

#' Convert TMB output from \link{\code{run_mcmc}} into a \code{shinystan}
#' object.
#'
#' @details The shinystan packages provides several conversion functions
#'   for objects of different types, such as stanfit classes (Stan ouput)
#'   and simple arrays. For the latter, option NUTS information, such as
#'   \code{sampler_params} can be passed. This function essentially extends
#'   the functionality of \code{as.shinystan} to work specifically with TMB
#'   MCMC lists. The user can thus explore their TMB model with
#'   \code{launch_shinystan(as.shinystan.tmb(tmb.fit))} in the same way
#'   that Stan models are examined.
#' @param tmb.fit Output list from \link{\code{run_mcmc}} for any of the
#' three algorithms.
#' @seealso launch_shinystan_tmb
#' @return An S4 object of class shinystan. Depending on the algorithm
#'   used, this list will have slight differences.
#' @export
as.shinystan.tmb <- function(tmb.fit){
  if(tmb.fit$algorithm=="NUTS"){
    sso <- with(tmb.fit, as.shinystan(samples, warmup=warmup, max_treedepth=max_treedepth,
             sampler_params=sampler_params, algorithm='NUTS', model_name=model))
  } else if(tmb.fit$algorithm=="HMC"){
    sso <- with(tmb.fit, as.shinystan(samples, warmup=warmup,
             sampler_params=sampler_params, algorithm='HMC', model_name=model))
  } else {
    sso <- with(tmb.fit, as.shinystan(samples, warmup=warmup,
             algorithm='RWM', model_name=model))
  }
  return(sso)
}

#' A high level wrapper to launch shinystan for a TMB MCMC list object.
#'
#' @details This function simply calls
#'   \code{launch_shinystan(as.shinystan.tmb(tmb.fit))}.
#' @export
launch_shinystan_tmb <- function(tmb.fit){
  launch_shinystan(as.shinystan.tmb(tmb.fit))
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

#' @export
as.shinystan.admb <- function(admb.fit){
  sso <- with(admb.fit,
    as.shinystan(samples, warmup=warmup, max_treedepth=max_treedepth,
                 sampler_params=sampler_params, algorithm='NUTS',
                 model_name=model))
  ## pars2 <- array(0, dim=c(nrow(pars), 1, ncol(pars)))
  ## pars2[,1,] <- as.matrix(pars)
  ## dimnames(pars2) <-
  ##   list(iter=1:nrow(pars), chains="chain:1",
  ##        parameters=dimnames(pars)[[2]])
  ## ss <- monitor(sims=pars2)
  ## y <- vector("list", length=length(dimnames(pars2)[[3]]))
  ## names(y) <- dimnames(pars2)[[3]]
  ## z <- lapply(y, function(x) x=numeric(0))
  ## sso <-
  ##   shinystan:::shinystan(
  ##   model_name=model.name, param_names=names(pars), param_dims=z,
  ##   posterior_sample=pars2, sampler_params=list(adapt),
  ##   summary=ss, n_chain=1, n_iter=nrow(pars),
  ##   n_warmup=nrow(pars)/2, model_code='NA',
  ##   misc=list(max_td=12, stan_method='sampling',
  ##             stan_algorithm='NUTS',
  ##             sso_version=utils::packageVersion('shinystan')))

    return(sso)
}

#' @export
launch_shinystan_admb <- function(admb.fit){
  launch_shinystan(as.shinystan.admb(admb.fit))
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

write.admb.cov <- function(cov.unbounded, model.path=getwd()){
  temp <- file.exists(paste0(model.path, "/admodel.cov"))
  if(!temp) stop(paste0("Couldn't find file ",model.path, "/admodel.cov"))
  temp <- file.copy(from=paste0(model.path, "/admodel.cov"),
                    to=paste0(model.path, "/admodel_original.cov"))
  wd.old <- getwd()
  setwd(model.path)
  ## Read in the output files
  results <- get.admb.cov()
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
  writeBin(as.integer(results$hybrid_bounded_flag), con=file.new)
  writeBin(as.vector(scale), con=file.new)
}
