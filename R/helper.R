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
#' the R package \link{rstan}.
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
#' @param inc_lp Whether to drop the column of log posterior density (last
#'   column). For diagnostics it should be included.
#' @return An invisible data.frame containing samples (rows) of each
#'   parameter (columns). If multiple chains exist they will be rbinded
#'   together.
#' @export
extract_samples <- function(fit.tmb, inc_warmup=FALSE, inc_lp=FALSE){
  x <- fit.tmb$samples
  if(!is.array(x)) stop("fit.tmb$samples is not an array -- valid TMB output?")
  ind <- if(inc_warmup) 1:dim(x)[1] else -(1:fit.tmb$warmup)
  ## Drop LP
  if(inc_lp){
  y <- do.call(rbind, lapply(1:dim(x)[2], function(i) x[ind, i,]))
  } else {
  y <- do.call(rbind, lapply(1:dim(x)[2], function(i) x[ind, i, -dim(x)[3]]))
  }
  return(invisible(as.data.frame(y)))
}

#' Extract sampler parameters from a fit
#'
#' @param fit A list returned by \code{\link{run_mcmc}}.
#' @param inc_warmup Whether to extract the warmup samples or not
#'   (default). Warmup samples should never be used for inference, but may
#'   be useful for diagnostics.
#' @return An invisible data.frame containing samples (rows) of each
#'   parameter (columns). If multiple chains exist they will be rbinded
#'   together.
#' @export
extract_sampler_params <- function(fit, inc_warmup=FALSE){
  x <- fit$sampler_params
  if(!is.list(x)) stop("fit$sampler_parameters is not a list -- valid output?")
  ind <- if(inc_warmup) 1:dim(x)[1] else -(1:fit$warmup)
  y <- do.call(rbind, lapply(1:length(x), function(i) x[[i]][ind,]))
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
  psv <- file.path(model.path, paste0(fn, '.psv'))
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


#' Plot pairs for MCMC output from an ADMB model.
#'
#' This function is useful for checking covergence and posterior
#' properties.
#'
#' @param posterior Dataframe containing the MCMC output, as read in using
#'   function \link{\code{run.mcmc}}
#' @param mle A list as read in by \code{\link{r4ss::read.admbFit}}. It
#'   uses the parameter estimates and covariance and correlation matrices
#'   as estimated asymptotically.
#' @param diag What type of plot to include on the diagonal, options are
#'   'acf' which plots the autocorrelation function \code{acf}, 'hist'
#'   shows marginal posterior histograms, and 'trace' the trace plot.
#' @param pars A vector of parameter names or integers representing which
#'   parameters to subset. Useful if the model has a larger number of
#'   parameters and you just want to show a few key ones.
#' @param acf.ylim If using the acf function on the diagonal, specify the y
#'   limit. The default is c(-1,1).
#' @param ymult A vector of length ncol(posterior) specifying how much room
#'   to give when using the hist option for the diagonal. For use if the
#'   label is blocking part of the plot. The default is 1.3 for all
#'   parameters.
#' @param limits A list containing the ranges for each parameter to use in
#'   plotting.
#' @return Produces a plot, and returns nothing.
#' @author Cole Monnahan
#' @export
pairs_admb <- function(posterior, mle, divergences=NULL, diag=c("acf","hist", "trace"),
                       acf.ylim=c(-1,1), ymult=NULL, axis.col=gray(.5),
                       pars=NULL,label.cex=.5, limits=NULL, ...){
  ## reset to old par when exiting
  old.par <- par(no.readonly=TRUE)
  on.exit(par(old.par))
  diag <- match.arg(diag)
  if(!(NCOL(posterior) %in% c(mle$nopar, mle$nopar+1)))
    stop("Number of parameters in posterior and mle not the same")
  ## pars will either be NULL, so use all parameters. OR a vector of
  ## indices OR a vector of characters. Want to force lp__ to be the very
  ## last one stylistically, and b/c there is no ellipse for it.
  if(is.null(pars)){
    ## Use all or first 50
    if(NCOL(posterior)>50){
      warning("Only showing first 50 parameters, use 'pars' argument to adjust")
      pars <- 1:50
    } else {
      pars <- 1:NCOL(posterior)
    }
  } else if(is.character(pars[1])){
    ## Drop lp__ here and force below to be the very last column
    pars <- pars[pars!='lp__']
    pars <- match(x=pars, table=names(posterior))
    if(any(is.na(pars))){
      warning("Some par names did not match -- dropped")
      pars <- pars[!is.na(pars)]
    }
  }
  ## pars is now vector of column indices of the parameters, so get the
  ## pieces needed.
  lp <- which(names(posterior)=='lp__')
  pars <- pars[pars!=lp]
  par.names <- mle$names[pars]
  mle.par <- mle$est[pars]
  mle.se <- mle$se[pars]
  mle.cor <- mle$cor[pars, pars]
  ## Now force log density to be last one.
  pars <- c(pars, lp)
  par.names <- c(par.names, 'log-density')
  n <- length(pars)
  if(n==1) stop("This function is only meaningful for >1 parameter")
  posterior <- posterior[,pars]
  if(is.null(ymult)) ymult <- rep(1.3, n)
  ## If no limits given, calculate the max range of the posterior samples and
  ## parameter confidence interval
  if(is.null(limits)){
    limits <- list()
    for(i in 1:n){
      limit.temp <- if(i==n) NULL else mle.par[i]+c(-1,1)*1.96*mle.se[i]
      ## multiplier for the ranges, adjusts the whitespace around the
      ## plots
      min.temp <- min(posterior[,i], limit.temp[1])
      max.temp <- max(posterior[,i], limit.temp[2])
      margin <- .15*(max.temp-min.temp)
      limits[[i]] <- c(min.temp-margin, max.temp+margin)
    }
  }
  ## Change posterior point look depending on how many samples. Makes
  ## easier to read.
  N <- NROW(posterior)
  mypch <- 16; mycol <- 1
  if(N>=1000){
    mycol <- rgb(0,0,0,.5)
  } else if(N>=10000){
    mycol <- rgb(0,0,0,.05)
  }
  if(is.null(divergences)) divergences <- rep(0, N)
  cexs <- ifelse(divergences, .25, .1)
  cols <- ifelse(divergences, rgb(1,0,0), mycol)
  par(mfrow=c(n,n), mar=0*c(.1,.1,.1,.1), yaxs="i", xaxs="i", mgp=c(.25, .25,0),
      tck=-.02, cex.axis=.65, col.axis=axis.col, oma=c(2, 2, 2,.5))
  temp.box <- function() box(col=axis.col, lwd=.5)
  ## Row and col here are not the posterior, but the matrix of pairwise
  ## combinations
  for(row in 1:n){
    for(col in 1:n){
      ## Diagonal, so add user choice
      if(row==col){
        if(diag=="hist"){
          h <- hist(posterior[,row], plot=F)
          ## Annoyingling you can't pass NULL to xlim in hist. So
          ## have to split up for two cases depending on limits.
          if(is.null(limits)){
            hist(posterior[,row], axes=F, freq=FALSE, ann=F,
                 ylim=c(0, ymult[row]*max(h$density)),
                 col=gray(.8), border=gray(.5))
          } else {
            ## Else use the user provided limits
            hist(posterior[,row], axes=F, freq=FALSE, ann=F,
                 ylim=c(0, ymult[row]*max(h$density)),
                 col=gray(.8), border=gray(.5), xlim=limits[[row]])
          }
          temp.box()
        } else if(diag=="acf") {
          acf(posterior[,row], axes=F, ann=F, ylim=acf.ylim)
          ## legend("topright", bty='n', legend=NA,
          ##        title=sprintf("EFS=%.3f", 100*admb_mcmc$diag$efsize[row],2))
          temp.box()
        } else if(diag=="trace") {
          plot(x=posterior[,row], lwd=.5, col=gray(.5), type="l", axes=F,
               ann=F, ylim=limits[[row]])
          temp.box()
        }
      }
      ## If lower triangle add scatterplot
      if(row>col){
        par(xaxs="r", yaxs="r")
        plot(x=posterior[,col], y=posterior[,row], axes=FALSE, ann=FALSE,
             pch=mypch, cex=cexs, col=cols, xlim=limits[[col]], ylim=limits[[row]], ...)
        if(row != n){
          ## Add bivariate 95% normal levels for both the MLE
          ## estimated covariance, but also the user supplied cov.user
          points(x=mle.par[col], y=mle.par[row],
                 pch=16, cex=.1, col=2)
          ## Get points of a bivariate normal 95% confidence contour
          ellipse.temp <- ellipse::ellipse(x=mle.cor[col, row],
                                           scale=mle.se[c(col, row)],
                                           centre= mle.par[c(col, row)], npoints=1000,
                                           level=.95)
          lines(ellipse.temp , lwd=1.5, lty=1, col="red")
        }
        par(xaxs="i", yaxs="i")
        temp.box()
      }
      if(row<col){
        ## If upper triangle add text showing the empirical correlation
        plot(0,0,type="n", xlim=c(0,1), ylim=c(0,1), axes=F,ann=F)
        temp.cor <- round(cor(posterior[,c(row,col)])[1,2],2)
        ## Set a minimum limit for this, so they're still
        ## visible, but still a function of correlation. This
        ## might need to be dynamic with n.
        legend("center", legend=NA, title=temp.cor,
               cex=(3*abs(temp.cor)+.25)*.5, bty='n')
        ## text(.5,.5, labels=temp.cor, cex=(3*abs(temp.cor)+.5)*.9,
        ##      col=1)
        temp.box()
      }
      ## Add special cases of axes on the ends
      if(row==n) {
        par( mgp=c(.05, ifelse(col %% 2 ==0, 0, .5),0) )
        axis(1, col=axis.col, lwd=.5)
      }
      if(col==1 & row >1) {
        par( mgp=c(.05, ifelse(row %% 2 ==1, .15, .65),0) )
        axis(2, col=axis.col, lwd=.5)
      }
      if(col==1 & row ==1){
        par( mgp=c(.05, ifelse(row %% 2 ==1, .15, .65),0) )
        axis(2, col=axis.col, lwd=.5)
      }
      if(row==1) mtext(par.names[col], line=.1, cex=label.cex)
    }
  }
}

#' Read maximum likelihood fit for ADMB model
#'
#' @param model Model name
#' @return A list containing, MLE estimates, standard errors, covariance
#'   and correlation matrices, and other output from ADMB.
#' @details This is based loosely off read.admbFit from r4ss.
#'
#' @export
read_mle_fit <- function(model, path=getwd()){
  oldwd <- getwd(); on.exit(setwd(oldwd))
  setwd(path)
  ## Sequentially read .par file which contains model size, minimum NLL,
  ## and maxgrad at the top
  par <- as.numeric(scan(paste(model,'.par', sep=''),
    what='', n=16, quiet=TRUE)[c(6,11,16)])
  nopar <- as.integer(par[1])
  nll <- par[2] #objective function value
  maxgrad <- par[3]

  ## The .cor file contains parameter (and derived quantity) names,
  ## estimates, and se's. This is more convenient to read in than the .par
  ## file.
  file <- paste(model,'.cor', sep='')
  xx <- readLines(file)
  ## Total parameter including sdreport variables
  totPar <- length(xx)-2
  ## Log of the determinant of the hessian
  logDetHess <- as.numeric(strsplit(xx[1], '=')[[1]][2])
  sublin <- lapply(strsplit(xx[1:totPar+2], ' '),function(x)x[x!=''])
  names.all <- unlist(lapply(sublin,function(x)x[2]))
  names.all <- as.vector(do.call(c, sapply(unique(names.all), function(n){
    x <- names.all[names.all==n]
    if(length(x)==1) return(x)
    list(paste0(x, '[',1:length(x),']'))})))

  est <- as.numeric(unlist(lapply(sublin,function(x)x[3])))
  std <- as.numeric(unlist(lapply(sublin,function(x)x[4])))
  ## The correlation in the bounded space.
  cor <- matrix(NA, totPar, totPar)
  corvec <- unlist(sapply(1:length(sublin), function(i)sublin[[i]][5:(4+i)]))
  cor[upper.tri(cor, diag=TRUE)] <- as.numeric(corvec)
  cor[lower.tri(cor)]  <-  t(cor)[lower.tri(cor)]
  ## Covariance matrix
  cov <- cor*(std %o% std)
  result <- list(nopar=nopar, nll=nll, maxgrad=maxgrad,
                 par.names=names.all[1:nopar],
                 names.all=names.all,
                 est=est, se=std, cor=cor, cov=cov)
  return(result)
}
