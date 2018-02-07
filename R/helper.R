## Read in PSV file
.get_psv <- function(model){
      if(!file.exists(paste0(model, '.psv'))){
      ## Sometimes ADMB will shorten the name of the psv file for some
      ## reason, so need to catch that here.
      ff <- list.files()[grep(x=list.files(), pattern='psv')]
      if(length(ff)==1){
        warning(paste("No .psv file found, using", ff))
        pars <- R2admb::read_psv(sub('.psv', '', x=ff))
      } else {
        stop(paste("No .psv file found -- did something go wrong??"))
      }
    } else {
      ## If model file exists
      pars <- R2admb::read_psv(model)
    }
  return(pars)
}

## Update algorithm for mass matrix.
##
## @param fn The current fn function.
## @param gr The current gr function
## @param y.cur The current parameter vector in unrotated (Y) space.
## @param M The new mass matrix
.rotate_space <- function(fn, gr, M,  y.cur){
  ## Rotation done using choleski decomposition
  ## First case is a dense mass matrix
  if(is.matrix(M)){
    chd <- t(chol(M))               # lower triangular Cholesky decomp.
    chd.inv <- solve(chd)               # inverse
    ## Define rotated fn and gr functions
    fn2 <- function(x) fn(chd %*% x)
    gr2 <- function(x) {as.vector( gr(chd %*% x) %*% chd )}
    ## Now rotate back to "x" space using the new mass matrix M
    x.cur <- as.numeric(chd.inv %*% y.cur)
  } else if(is.vector(M)){
    chd <- sqrt(M)
    fn2 <- function(x) fn(chd * x)
    gr2 <- function(x) as.vector(gr(chd * x) ) * chd
    ## Now rotate back to "x" space using the new mass matrix M. M is a
    ## vector here. Note the big difference in efficiency without the
    ## matrix operations.
    x.cur <- (1/chd) * y.cur
  } else {
    stop("Mass matrix must be vector or matrix")
  }
  ## Redefine these functions
  ## Need to adjust the current parameters so the chain is
  ## continuous. First rotate to be in Y space.
  return(list(gr2=gr2, fn2=fn2, x.cur=x.cur, chd=chd))
}

## Update the control list.
##
## @param control A list passed from \code{sample_tmb}.
## @return A list with default control elements updated by those supplied
##   in \code{control}
.update_control <- function(control){
  default <- list(adapt_delta=0.8, metric=NULL, stepsize=NULL,
                  adapt_mass=TRUE, max_treedepth=12, w1=75, w2=50, w3=25)
  if(is.matrix(control$metric) & !is.null(control$adapt_mass)){
    if(control$adapt_mass==TRUE){
      warning("Mass matrix adaptation disabled if metric is a matrix")
    }
    control$adapt_mass <- FALSE
  }
  new <- default
  if(!is.null(control))
    for(i in names(control))  new[[i]] <- control[[i]]
  if(is.matrix(new$metric)) new$adapt_mass <- FALSE
  return(new)
}

## Print MCMC progress to console.
##
## @param iteration The iteration of the MCMC chain.
## @param iter The total iterations.
## @param warmup The number of warmup iterations.
## @param chain The chain being run (bookkeeping only).
## @return Nothing. Prints to message to console.
##
## @details This function was modeled after the functionality provided by
## the R package rstan.
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

## Print MCMC timing to console
## @param time.warmup Time of warmup in seconds.
## @param time.total Time of total in seconds.
## @return Nothing. Prints message to console.
##
## @details This function was modeled after the functionality provided by
##   the R package \pkg{rstan}.
.print.mcmc.timing <- function(time.warmup, time.total){
  x <- ' Elapsed Time: '
  message(paste0(x, sprintf("%.1f", time.warmup), ' seconds (Warmup)'))
  message(paste0(x, sprintf("%.1f", time.total-time.warmup), ' seconds (Sampling)'))
  message(paste0(x, sprintf("%.1f", time.total), ' seconds (Total)'))
}

## Convert adnuts fit (named list) into a \code{shinystan} object.
##
## @details The shinystan packages provides several conversion functions
##   for objects of different types, such as stanfit classes (Stan ouput)
##   and simple arrays. For the latter, option NUTS information, such as
##   \code{sampler_params} can be passed. This function essentially extends
##   the functionality of \code{as.shinystan} to work specifically with
##   fits from adnuts (TMB or ADMB). The user can thus explore their model
##   with \code{launch_shinystan(as.shinystan.tmb(fit))} in the same way
##   that Stan models are examined.
## @param fit Output list from \code{sample_tmb} or
##   \code{sample_admb}.
## @seealso launch_shinytmb, launch_shinyadmb
## @return An S4 object of class shinystan. Depending on the algorithm
##   used, this list will have slight differences.
.as.shinyadnuts <- function(fit){
  if(fit$algorithm=="NUTS"){
    sso <- with(fit, shinystan::as.shinystan(samples, warmup=warmup, max_treedepth=max_treedepth,
             sampler_params=sampler_params, algorithm='NUTS', model_name=model))
  } else if(fit$algorithm=="HMC"){
    sso <- with(fit, shinystan::as.shinystan(samples, warmup=warmup,
             sampler_params=sampler_params, algorithm='HMC', model_name=model))
  } else {
    sso <- with(fit, shinystan::as.shinystan(samples, warmup=warmup,
             algorithm='RWM', model_name=model))
  }
  return(invisible(sso))
}

#' Launch shinystan for a TMB fit.
#'
#' @param fit A named list returned by \code{sample_tmb}.
#' @seealso \code{launch_shinyadmb}
#' @export
launch_shinytmb <- function(fit){
  shinystan::launch_shinystan(.as.shinyadnuts(fit))
}

#' Launch shinystan for an ADMB fit.
#'
#' @param fit A named list returned by \code{sample_admb}.
#' @seealso \code{launch_shinytmb}
#' @export
launch_shinyadmb <- function(fit){
  shinystan::launch_shinystan(.as.shinyadnuts(fit))
}


#' Extract posterior samples from a model fit.
#'
#' A helper function to extract posterior samples across multiple chains
#' into a single data.frame.
#'
#' @details This function is loosely based on the \pkg{rstan} function
#'   \code{extract}. Merging samples across chains should only be used for
#'   inference after appropriate diagnostic checks. Do not calculate
#'   diagnostics like Rhat or effective sample size after using this
#'   function, instead, use \code{\link[rstan]{monitor}}. Likewise, warmup
#'   samples are not valid and should never be used for inference, but may
#'   be useful in some cases for diagnosing issues.
#'
#' @param fit A list returned by \code{sample_tmb} or \code{sample_admb}.
#' @param inc_warmup Whether to extract the warmup samples or not
#'   (default). Warmup samples should never be used for inference, but may
#'   be useful for diagnostics.
#' @param inc_lp Whether to include a column for the log posterior density
#'   (last column). For diagnostics it can be useful.
#' @param as.list Whether to return the samples as a list (one element per
#'   chain). This could then be converted to a CODA mcmc object.
#' @return If as.list is FALSE, an invisible data.frame containing samples
#'   (rows) of each parameter (columns). If multiple chains exist they will
#'   be rbinded together, maintaining order within each chain. If as.list
#'   is TRUE, samples are returned as a list of matrices.
#' @export
#' @examples
#' ## A previously run fitted TMB model
#' fit <- readRDS(system.file('examples', 'fit_tmb.RDS', package='adnuts'))
#' post <- extract_samples(fit)
#' tail(apply(post, 2, median))
extract_samples <- function(fit, inc_warmup=FALSE, inc_lp=FALSE, as.list=FALSE){
  x <- fit$samples
  if(!is.array(x)) stop("fit$samples is not an array -- valid fit object?")
  ind <- if(inc_warmup) 1:dim(x)[1] else -(1:fit$warmup)
  ## Drop LP
  if(inc_lp){
    y <-  lapply(1:dim(x)[2], function(i) x[ind, i,])
  } else {
    y <-  lapply(1:dim(x)[2], function(i) x[ind, i, -dim(x)[3]])
  }
  if(as.list){
    return(invisible(y))
  } else {
    return(invisible(as.data.frame(do.call(rbind, y))))
  }
}


#' Extract sampler parameters from a fit.
#'
#' Extract information about NUTS trajectories, such as acceptance ratio
#' and treedepth, from a fitted object.
#'
#' @details Each trajectory (iteration) in NUTS has associated information
#'   about the trajectory: stepsize, acceptance ratio, treedepth, and number of
#'   leapfrog steps. This function extracts these into a data.frame, which
#'   may be useful for diagnosing issues in certain cases. In general, the
#'   user should not need to examine them, or preferably should via
#'   \code{\link{launch_shinytmb}} or \code{\link{launch_shinyadmb}}.
#'
#' @param fit A list returned by \code{sample_admb} or \code{sample_tmb}.
#' @param inc_warmup Whether to extract the warmup samples or not
#'   (default). Warmup samples should never be used for inference, but may
#'   be useful for diagnostics.
#' @return An invisible data.frame containing samples (rows) of each
#'   parameter (columns). If multiple chains exist they will be rbinded
#'   together.
#' @seealso \code{\link{launch_shinytmb}} and \code{\link{launch_shinyadmb}}.
#' @export
#' @examples
#' fit <- readRDS(system.file('examples', 'fit_tmb.RDS', package='adnuts'))
#' ## Examine how step size and treedepth changes as the mass matrix updates
#' ## during warmup
#' sp <- extract_sampler_params(fit, inc_warmup=TRUE)
#' plot(0,0, type='n', xlim=c(0,510), ylim=c(0,1), xlab='Iteration',
#'      ylab='Step size (eps)')
#' for(i in 1:3) lines(1:1000, sp[sp$chain==i,4], col=i)
#' legend('topright', cex=.7, legend=paste("chain1", 1:3), lty=1, col=1:3)
#' plot(0,0, type='n', xlim=c(0,1000), ylim=c(0,10), xlab='Iteration',
#'      ylab='Treedepth')
#' for(i in 1:3) lines(1:1000, sp[sp$chain==i,5], col=i)
#' legend('topright', cex=.7, legend=paste("chain1", 1:3), lty=1, col=1:3)
#'
extract_sampler_params <- function(fit, inc_warmup=FALSE){
  x <- fit$sampler_params
  if(!is.list(x)) stop("fit$sampler_parameters is not a list -- valid fit object?")
  if(inc_warmup){
    ind <- 1:dim(x[[1]])[1]
    its <- 1:length(ind)
  } else{
    ind <- -(1:fit$warmup)
    its <- (1:length(ind)) + fit$warmup
  }
  y <- do.call(rbind, lapply(1:length(x), function(i)
    cbind(chain=i, iteration=its, x[[i]][ind,])))
  return(invisible(as.data.frame(y)))
}

## Write matrix of samples to a binary .psv file.
##
## @details Useful to combine multiple MCMC runs together into a single
## .psv file which can then be executed with '-mceval'.
## @param fn Model name
## @param samples A matrix or data.frame of samples, each column is a
##   parameter, each row a sample.
.write_psv <- function(fn, samples, model.path=getwd()){
  samples <- as.matrix(samples)
  psv <- file.path(model.path, paste0(fn, '.psv'))
  con <- file(psv, 'wb')
  writeBin(object=ncol(samples), con=con)
  writeBin(object=as.vector(t(samples)), con=con)
  close(con)
}

## Read in the ADMB covariance file.
##
## @param model.path Path to model (defaults to working directory)
## @export
.get.admb.cov <- function(model.path=getwd()){
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

## Write a covariance matrix to admodel.cov.
##
## @param cov.unbounded The cov matrix in unbounded space.
## @param hbf The hybrid_bounded_flag value. Use hbf=1 for HMC.
## @param model.path Path to model.
.write.admb.cov <- function(cov.unbounded, model.path=getwd(), hbf=NULL){
  temp <- file.exists(paste0(model.path, "/admodel.cov"))
  if(!temp) stop(paste0("Couldn't find file ",model.path, "/admodel.cov"))
  temp <- file.copy(from=paste0(model.path, "/admodel.cov"),
                    to=paste0(model.path, "/admodel_original.cov"))
  wd.old <- getwd()
  setwd(model.path)
  ## Read in the output files
  results <- .get.admb.cov()
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


## Read maximum likelihood fit for ADMB model
##
## @param model Model name
## @return A list containing, MLE estimates, standard errors, covariance
##   and correlation matrices, and other output from ADMB.
## @details This is based loosely off read.admbFit from r4ss.
##
## @export
.read_mle_fit <- function(model, path=getwd()){
  oldwd <- getwd(); on.exit(setwd(oldwd))
  setwd(path)
  ## Sequentially read .par file which contains model size, minimum NLL,
  ## and maxgrad at the top
  f <- paste(model,'.par', sep='')
  if(!file.exists(f)){
    warning(paste("File", f,
                  "not found so could not read in MLE quantities or parameter names"))
    return(NULL)
  }
  par <- as.numeric(scan(f, what='', n=16, quiet=TRUE)[c(6,11,16)])
  nopar <- as.integer(par[1])
  nll <- par[2] #objective function value
  maxgrad <- par[3]

  ## The .cor file contains parameter (and derived quantity) names,
  ## estimates, and se's. This is more convenient to read in than the .par
  ## file.
  f <- paste(model,'.cor', sep='')
  if(!file.exists(f)){
    warning(paste("File", f,
                  "not found so could not read in MLE quantities or parameter names"))
    return(NULL)
  }
  xx <- readLines(f)
  ## Total parameter including sdreport variables
  totPar <- length(xx)-2
  if(totPar < nopar) {
    warning(paste("File", f,
                  "did not match the .cor file.. maybe hessian failed? MLE object not available"))
    return(NULL)
  }
  ## Log of the determinant of the hessian
  logDetHess <- as.numeric(strsplit(xx[1], '=')[[1]][2])
  sublin <- lapply(strsplit(xx[1:totPar+2], ' '),function(x)x[x!=''])
  names.all <- unlist(lapply(sublin,function(x)x[2]))[1:nopar]
  names.all <- as.vector(do.call(c, sapply(unique(names.all), function(n){
    x <- names.all[names.all==n]
    if(length(x)==1) return(list(x))
    list(paste0(x, '[',1:length(x),']'))})))
  est <- as.numeric(unlist(lapply(sublin,function(x)x[3])))
  std <- as.numeric(unlist(lapply(sublin,function(x)x[4])))
  ## The correlation in the bounded space.
  cor <- matrix(NA, totPar, totPar)
  corvec <- unlist(sapply(1:length(sublin), function(i)sublin[[i]][5:(4+i)]))
  cor[upper.tri(cor, diag=TRUE)] <- as.numeric(corvec)
  cor[lower.tri(cor)]  <-  t(cor)[lower.tri(cor)]
  ## Covariance matrix
  ## cov <- cor*(std %o% std)
  result <- list(nopar=nopar, nll=nll, maxgrad=maxgrad,
                 par.names=names.all[1:nopar],
                 names.all=names.all,
                 est=est, se=std, cor=cor[1:nopar,1:nopar])
  return(result)
}
