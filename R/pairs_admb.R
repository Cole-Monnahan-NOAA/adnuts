#' Plot pairwise parameter posteriors and optionally the MLE points and
#' confidence ellipses.
#'
#' @param fit A list as returned by \code{sample_admb}.
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
#' @param label.cex Control size of labels
#' @param axis.col Color of axes
#' @param ... Arguments to be passed to plot call in lower diagonal panels
#' @param limits A list containing the ranges for each parameter to use in
#'   plotting.
#' @return Produces a plot, and returns nothing.
#' @export
#' @author Cole Monnahan
pairs_admb <- function(fit,
                       diag=c("trace","acf","hist"),
                       acf.ylim=c(-1,1), ymult=NULL, axis.col=gray(.5),
                       pars=NULL, label.cex=.5, limits=NULL, ...){
  mle <- fit$mle
  posterior <- extract_samples(fit, inc_lp=TRUE)
  chains <- rep(1:dim(fit$samples)[2], each=dim(fit$samples)[1]-fit$warmup)
  divs <- if(fit$algorithm=="NUTS")
            extract_sampler_params(fit)$divergent__ else NULL
  ptcex <- .1
  divcex <- .5
  chaincols <- 1:length(unique(chains))
  wp <- function(par.name) {
    ## Temporary function that finds (w)hich (p)arameter position par.name
    ## is in for the MLE. These will not always match due to things like
    ## log-posterior and derived quantities, in this case, MLEs are NA's
    x <- which(mle$par.names == par.name)
    if(length(x)==0) return(NA)
    if(length(x)>1) stop("Par matched multiple??")
    return(x)
  }
  ## reset to old par when exiting
  old.par <- par(no.readonly=TRUE)
  on.exit(par(old.par))
  diag <- match.arg(diag)
  par.names <- names(posterior)
  ## if(!(NCOL(posterior) %in% c(mle$nopar, mle$nopar+1)))
  ##   stop("Number of parameters in posterior and mle not the same")
  ## pars will either be NULL, so use all parameters. OR a vector of
  ## indices OR a vector of characters. Want to force lp__ to be the very
  ## last one stylistically, and b/c there is no ellipse for it.
  if(is.null(pars)){
    ## Use all or first 10
    if(NCOL(posterior)>10){
      warning("Only showing first 10 parameters, use 'pars' argument to adjust")
      pars <- par.names[1:10]
    } else {
      pars <- par.names[1:NCOL(posterior)]
    }
  } else if(is.numeric(pars[1])){
    ## Index can be used instead of character names
    pars <- par.names[pars]
  }
  pars.bad <- match(x=pars, table=names(posterior))
  if(any(is.na(pars.bad))){
    warning("Some par names did not match -- dropped")
    pars <- pars[!is.na(pars.bad)]
  }
  n <- length(pars)
  if(n==1) stop("This function is only meaningful for >1 parameter")
  posterior <- posterior[,pars]
  if(is.null(ymult)) ymult <- rep(1.3, n)
  ## If no limits given, calculate the max range of the posterior samples and
  ## parameter confidence interval
  if(is.null(limits)){
    limits <- list()
    for(i in 1:n){
      pp <- wp(pars[i])
      limit.temp <- if(is.na(pp)) NULL else mle$est[pp]+c(-1,1)*1.96*mle$se[pp]
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
  if(is.null(divs)) divs <- rep(0, N)
  par(mfrow=c(n,n), mar=0*c(.1,.1,.1,.1), yaxs="i", xaxs="i", mgp=c(.25, .25,0),
      tck=-.02, cex.axis=.65, col.axis=axis.col, oma=c(2, 2, 2,2))
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
          temp.box()
        } else if(diag=="trace") {
          ## Trace plots for each chain separately
          xlim <- c(1, length(chains[chains==1]))
          plot(x=0, y=0,  type="n", axes=FALSE,
               ann=FALSE, ylim=limits[[row]], xlim=xlim)
          for(ll in unique(chains)){
            lines(posterior[chains==ll,row], col=chaincols[ll], lwd=.1)
          }
          temp.box()
        }
      }
      ## If lower triangle and covariance known, add scatterplot
      if(row>col){
        par(xaxs="r", yaxs="r")
        plot(x=posterior[,col], y=posterior[,row], axes=FALSE, ann=FALSE,
             pch=mypch, cex=ptcex, col=mycol, xlim=limits[[col]],
             ylim=limits[[row]], ...)
        ## replot divegences on top so they are always visible
        points(x=posterior[which(divs==1),col], y=posterior[which(divs==1),row],
               pch=mypch, cex=divcex, col='red')
        p1 <- wp(pars[row]); p2 <- wp(pars[col])
        if(!is.na(p1) & !is.na(p2)){
          ## Add bivariate 95% normal levels from MLE
          points(x=mle$est[p2], y=mle$est[p1],
                 pch=16, cex=.5, col=2)
          ## Get points of a bivariate normal 95% confidence contour
          if(!requireNamespace("ellipse", quietly=TRUE)){
            warning("ellipse package needs to be installed to show ellipses")
          }
          else {
            ellipse.temp <- ellipse(x=mle$cor[p2, p1],
                                    scale=mle$se[c(p2, p1)],
                                    centre= mle$est[c(p2, p1)], npoints=1000,
                                    level=.95)
            lines(ellipse.temp , lwd=.5, lty=1, col="red")
          }
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
      if(row==1) mtext(pars[col], line=ifelse(col %% 2 ==1, .1, 1.1),
                       cex=label.cex)
      if(col==n)
        mtext(pars[row], side=4, line=ifelse(row %% 2 ==1, 0, 1), cex=label.cex)
    }
  }
}
