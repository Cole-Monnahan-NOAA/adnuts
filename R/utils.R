
#' Constructor for the "adfit" (A-D fit) class
#' @param x Fitted object from \code{\link{sample_admb}}
#' @return An object of class "adfit"
#' @export
adfit <- function(x){
  stopifnot(is.list(x))
  if(is.null(x$samples)) stop("Samples missing from fit")
  if(is.null(x$algorithm)) stop("Algorithm missing from fit")
  class(x) <- 'adfit'
  x
}

#' Check object of class adfit
#' @param x Returned list from \code{\link{sample_admb}}
#' @export
is.adfit <- function(x) inherits(x, "adfit")


#' Plot object of class adfit
#' @param x Fitted object from \code{\link{sample_admb}}
#' @param y Ignored
#' @param ... Ignored
#' @return Plot created
#' @method plot adfit
#' @export
plot.adfit <- function(x, y, ...) plot_marginals(x)

#' Print summary of object of class adfit
#' @param object Fitted object from \code{\link{sample_admb}}
#' @param ... Ignored
#' @return Summary printed to screen
#' @method summary adfit
#' @export
summary.adfit <- function(object, ...) print(object)

#' Print summary of adfit object
#' @param x Fitted object from \code{\link{sample_admb}}
#' @param ... Ignored
#' @return Summary printed to console
#' @method print adfit
#' @export
print.adfit <- function(x, ...){
  iter <- dim(x$samples)[1]
  chains <- dim(x$samples)[2]
  pars <- dim(x$samples)[3]
  samples <- (iter-x$warmup)*chains
  cat(paste0("Model '", x$model,"'"), "has", pars,
      "pars, and was fit using", x$algorithm,
      "with", iter, "iter and", chains,
      "chains\n")
  rt <- sum(x$time.total)/chains
  ru <- 'seconds'
  if(rt>60){
    rt <- rt/60; ru <- 'minutes'
  } else if(rt>60*60) {
    rt <- rt/(60*60); ru <- 'hours'
  }
  cat("Average run time per chain was", round(rt,2),  ru, '\n')
  if(!is.null(x$monitor)){
    minESS <- min(x$monitor$n_eff)
    maxRhat <- round(max(x$monitor$Rhat),3)
    cat(paste0("Minimum ESS=",
                   minESS,
                   " (",
                   round(100*minESS/samples,2),
                   "%), and maximum Rhat=", maxRhat, '\n'))
  }
  if(x$algorithm=='NUTS'){
    ndivs <- sum(extract_sampler_params(x)[,'divergent__'])
    cat(paste0("There were ", ndivs, " divergences after warmup\n"))
  }
}



#' Plot marginal distributions for a fitted model
#'
#' @param fit A fitted object returned by
#'   \code{\link{sample_admb}}.
#' @param pars A numeric or character vector of parameters which
#'   to plot, for plotting a subset of the total (defaults to all)
#' @param mfrow A custom grid size (vector of two) to be called
#'   as \code{par(mfrow)}, overriding the defaults.
#' @param add.mle Whether to add marginal normal distributions
#'   determined from the inverse Hessian file
#' @param add.monitor Whether to add ESS and Rhat information
#' @param breaks The number of breaks to use in \code{hist()},
#'   defaulting to 30
#' @export
#'
#' @details This function plots grid cells of all parameters
#'   in a model, comparing the marginal posterior histogram vs
#'   the asympotitic normal (red lines) from the inverse
#'   Hessian. Its intended use is to quickly gauge differences
#'   between frequentist and Bayesian inference on the same
#'   model.
#'
#' If \code{fit$monitor} exists the effective sample size
#' (ESS) and R-hat estimates are printed in the top right
#' corner. See
#' \url{https://mc-stan.org/rstan/reference/Rhat.html} for more
#' information. Generally Rhat>1.05 or ESS<100 (per chain)
#' suggest inference may be unreliable.
#'
#' This function is customized to work with multipage PDFs,
#' specifically:
#' \code{pdf('marginals.pdf', onefile=TRUE, width=7,height=5)}
#' produces a nice readable file.
#' @examples
#' fit <- readRDS(system.file('examples', 'fit.RDS', package='adnuts'))
#' plot_marginals(fit, pars=1:2)
#'
plot_marginals <- function(fit, pars=NULL, mfrow=NULL,
                           add.mle=TRUE, add.monitor=TRUE,
                           breaks=30){
  if(!is.adfit(fit)) stop("fit is not a valid object")
  if(!is.null(mfrow)) stopifnot(is.vector(mfrow) && length(mfrow)==2)
  stopifnot(add.mle %in% c(TRUE,FALSE))
  if(add.mle & is.null(fit$mle)) {
    add.mle <- FALSE
    warning("No MLE information found in fit$mle so cannot add")
  }
  if(!add.mle) fit$mle <- NULL
  if(!add.monitor) fit$monitor <- NULL
  par.old <- par()
  on.exit(par(mfrow=par.old$mfrow, mar=par.old$mar,
              mgp=par.old$mgp, oma=par.old$oma, tck=par.old$tck))
  posterior <- extract_samples(fit, inc_lp=FALSE)
  par.names <- names(posterior)
  if(is.null(pars)) pars <- par.names
  if(is.character(pars[1])){
    pars.ind <- match(x=pars, table=par.names)
    if(any(is.na(pars.ind))){
      warning("Some par names did not match -- dropped")
      print(pars[is.na(pars.ind)])
      pars.ind <- pars.ind[!is.na(pars.ind)]
    }
    pars <- pars.ind
  } else if(any(pars > NCOL(posterior))){
    warning("Some par numbers too big -- dropped")
    print(pars[pars > NCOL(posterior)])
    pars <- pars[ pars <=NCOL(posterior)]
  }
  n <- length(pars)
  stopifnot(is.numeric(pars[1]))
  stopifnot(ncol(posterior)>1)
  par(mar=c(1.5,0,.1,0), mgp=c(2,.4,0),
      oma=c(.25,.25,.25,.25), tck=-.02)
  if(!is.null(mfrow)){
    par(mfrow=mfrow)
  } else if(n>12){
    par(mfrow=c(4,4))
  } else if(n>9){
    par(mfrow=c(4,3))
  } else if(n>6){
    par(mfrow=c(3,3))
  } else if(n>4){
    par(mfrow=c(3,2))
  } else if(n>3){
    par(mfrow=c(2,2))
  } else {
    par(mfrow=c(1,n))
  }
  for(ii in pars){
    par <- par.names[ii]
    if(!is.null(fit$mle)){
      mle <- fit$mle$est[ii]
      se <-  fit$mle$se[ii]
      x1 <- seq(qnorm(.001, mle, se), qnorm(.999, mle, se), len=100)
      y1 <- dnorm(x1, mle, se)
    } else{
      x1 <- y1 <- NULL
    }
    tmp <- hist(posterior[,ii], plot=FALSE, breaks=breaks)
    x2 <- tmp$mids; y2 <- tmp$density
    plot(0,0, type='n', xlim=range(c(x1,x2)), yaxs='i',
         ylim=c(0, max(c(y1,y2))*1.3), axes=FALSE, ann=FALSE)
    hist(posterior[,ii], breaks=breaks, add=TRUE, yaxs='i', freq=FALSE, col=gray(.8))
    axis(1);  box(col=gray(.5));
    if(!is.null(fit$mle)) lines(x1,y1, col='red', lwd=2)
    if(!is.null(fit$monitor)){
      mon <- fit$monitor
      ## add ESS and Rhat to top right
      tmp <- par("usr"); xy <- c(.85,.88)
      text.x <- tmp[1]+xy[1]*diff(tmp[1:2])
      text.y <- tmp[3]+xy[2]*diff(tmp[3:4])
      label <- paste0('ESS=', mon[ii,'n_eff'], "\nRhat=", round(mon[ii,'Rhat'],3))
      text(x=text.x, y=text.y, labels=label, cex=.8)
    }
    mtext(paste("",par), line=-1.6, adj=0, cex=.9)
  }
}



#' Plot adaptation metrics for a fitted model.
#'
#' @param fit A fitted object returned by
#' \code{\link{sample_admb}}.
#' @param plot Whether to plot the results
#' @return Prints and invisibly returns a ggplot object
#'
#' @details This utility function quickly plots the adaptation output of NUTS
#' chains.
#' @importFrom rlang .data
#' @export
#' @examples
#' fit <- readRDS(system.file('examples', 'fit.RDS', package='adnuts'))
#' plot_sampler_params(fit)
plot_sampler_params <- function(fit, plot=TRUE){
  if(!requireNamespace("ggplot2", quietly=TRUE))
    stop("ggplot2 package not found")
  sp <- adnuts::extract_sampler_params(fit, inc_warmup=TRUE)
  sp.long <-
    data.frame(iteration=sp$iteration, chain=factor(sp$chain),
               value=c(sp$accept_stat__, log(sp$stepsize__),
                       sp$n_leapfrog__, sp$divergent__, sp$energy__),
               variable=rep(c('accept_stat', 'log_stepsize',
                              'n_leapfrog', 'divergent',
                              'energy'), each=nrow(sp)))
  g <- ggplot2::ggplot(sp.long, ggplot2::aes(.data$iteration, y=.data$value, color=.data$chain)) +
    ggplot2::geom_point(alpha=.5) +
    ggplot2::facet_wrap('variable', scales='free_y', ncol=1) + ggplot2::theme_bw()
  if(plot) print(g)
  return(invisible(g))
}

#' Check that the file can be found
#'
#' @param model
#' @param path
.check_model_path <- function(model, path){
  stopifnot(is.character(path))
  stopifnot(is.character(model))
  if(!dir.exists(path))
    stop('Folder ', path, ' does not exist. Check argument \'path\'')
  if (.Platform$OS.type=="windows") {
    ff <- file.path(path, paste(model,".exe",sep=""))
  } else {
    ff <- file.path(path, paste("./",model,sep=""))
  }
  if(!file.exists(ff))
    stop('File ', ff, ' not found in specified folder. Check \'model\' argument')
}

#' Check that the  model is compiled with the right version
#' of ADMB which is 12.0 or later
#'
#' @param model Model name without file extension
#' @param path Path to model folder, defaults to working
#'   directory. NULL value specifies working directory (default).
#' @param min.version Minimum valid version (numeric). Defaults
#'   to 12.0.
#' @param warn Boolean whether to throw warnings or not
#' @return Nothing, errors out if either model could not be run
#'   or the version is incompatible. If compatible nothing
#'   happens.
#' @details Some functionality of packages \pkg{adnuts} is
#'   imbedded in the ADMB source code so that when a model is
#'   compiled it is contained in the model executable. If this
#'   code does not exist adnuts will fail. The solution is to
#'   update ADMB and recompile the model.
.check_ADMB_version <- function(model, path=getwd(),
                                min.version=12, warn=TRUE){
  if(!is.null(path)){
    if(dir.exists(path)){
    wd <- getwd()
    on.exit(setwd(wd))
    setwd(path)
    } else {
      stop("Invalid path, folder does not exist")
    }
  }
  ## Run the model to get the version info
  if (!.Platform$OS.type=="windows") {
    model <- paste0("./", model)
  }
  test <- try(system(paste(model, '-version'), intern=TRUE), silent=TRUE)
  if (inherits(test,"try-error"))
    stop(paste0("Could not detect version of ", model, ". Check executable and path"))
  ## v <- as.numeric(gsub('ADMB-', '', strsplit(test[3], ' ')[[1]][1]))
  v <- as.numeric(gsub('ADMB-', '', substr(strsplit(test[3], ' ')[[1]][1], 1,9)))
  if(is.na(v) | !is.numeric(v)){
    warning("Issue verifying ADMB version. Contact package mantainer")
    return(0)
  }
  if(v < min.version)
    stop(paste(model,"compiled with old version of ADMB. Version >12.0 required, found:\n", v,
               "\nadnuts is incompatible with this version. Update ADMB and try again"))
  if(v < 12.2 & warn){
    warning("This version contains bugs in the NUTS code. Consider updating ADMB to version at least 12.2")
  }
  return(invisible(v))
}


#' Function to generate random initial values from a previous fit using
#' adnuts
#'
#' @param fit An outputted list from \code{\link{sample_admb}}
#' @param chains The number of chains for the subsequent run, which
#'   determines the number to return.
#' @return A list of lists which can be passed back into
#'   \code{\link{sample_admb}}.
#' @export
sample_inits <- function(fit, chains){
  post <- extract_samples(fit)
  ind <- sample(1:nrow(post), size=chains)
  lapply(ind, function(i) as.numeric(post[i,]))
}

#' Read in admodel.hes file
#' @param path Path to folder containing the admodel.hes file
#'
#' @return The Hessian matrix
.getADMBHessian <- function(path){
  ## This function reads in all of the information contained in the
  ## admodel.hes file. Some of this is needed for relaxing the
  ## covariance matrix, and others just need to be recorded and
  ## rewritten to file so ADMB "sees" what it's expecting.
  filename <- file.path(path, "admodel.hes")
  if(!file.exists(filename))
    stop(paste0("admodel.hes not found: ", filename))
  f <- file(filename, "rb")
  on.exit(close(f))
  num.pars <- readBin(f, "integer", 1)
  hes.vec <- readBin(f, "numeric", num.pars^2)
  hes <- matrix(hes.vec, ncol=num.pars, nrow=num.pars)
  hybrid_bounded_flag <- readBin(f, "integer", 1)
  scale <- readBin(f, "numeric", num.pars)
  return(hes)
}


#' Check identifiability from model Hessian
#'
#' @param path Path to model folder, defaults to working directory
#' @param model Model name without file extension
#' @details Read in the admodel.hes file and check the eigenvalues to
#'   determine which parameters are not identifiable and thus cause the
#'   Hessian to be non-invertible. Use this to identify which parameters
#'   are problematic. This function was converted from a version in the
#'   \code{FishStatsUtils} package.
#' @return Prints output of bad parameters and invisibly returns it.
#' @export
check_identifiable <- function(model, path=getwd()){
  ## Check eigendecomposition
  fit <- .read_mle_fit(model, path)
  hes <- .getADMBHessian(path)
  ev  <-  eigen(hes)
  WhichBad <-  which( ev$values < sqrt(.Machine$double.eps) )
  if(length(WhichBad)==0){
    message( "All parameters are identifiable" )
  } else {
    ## Check for parameters
    if(length(WhichBad==1)){
      RowMax <- abs(ev$vectors[,WhichBad])
    } else {
      RowMax  <-  apply(ev$vectors[, WhichBad], MARGIN=1, FUN=function(vec){max(abs(vec))} )
    }
    bad <- data.frame(ParNum=1:nrow(hes), Param=fit$par.names,
                      MLE=fit$est[1:nrow(hes)],
                      Param_check=ifelse(RowMax>0.1, "Bad","OK"))
    row.names(bad) <- NULL
    bad <- bad[bad$Param_check=='Bad',]
    print(bad)
    return(invisible(bad))
  }
}


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
## @param control A list passed from a sampling function
## @return A list with default control elements updated by those supplied
##   in \code{control}
.update_control <- function(control){
  default <- list(adapt_delta=0.8, metric='unit', stepsize=NULL,
                  adapt_mass=TRUE, adapt_mass_dense=FALSE,
                  max_treedepth=12)
  ## Special case if user is doing mle they probably don't want
  ## mass adaptation turned on. They have to override it by
  ## setting TRUE for either adaptation option
  if(is.character(control$metric)| is.matrix(control$metric)){
    if(is.null(control$adapt_mass) &
       is.null(control$adapt_mass_dense)){
      default$adapt_mass <- default$adapt_mass_dense <- FALSE
    }
  }
  new <- default
  if(!is.null(control))
    for(i in names(control))  new[[i]] <- control[[i]]
  if(new$adapt_mass_dense & new$adapt_mass)
    new$adapt_mass <- FALSE
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
##   with \code{launch_shinystan(.as.shinyadnuts(fit))} in the same way
##   that Stan models are examined.
## @param fit Output list from  \code{sample_admb}.
## @seealso launch_shinyadmb
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
#' @param fit A list returned by \code{sample_admb}.
#' @param inc_warmup Whether to extract the warmup samples or not
#'   (default). Warmup samples should never be used for inference, but may
#'   be useful for diagnostics.
#' @param inc_lp Whether to include a column for the log posterior density
#'   (last column). For diagnostics it can be useful.
#' @param as.list Whether to return the samples as a list (one element per
#'   chain). This could then be converted to a CODA mcmc object.
#' @param unbounded Boolean flag whether to return samples in
#'   unbounded (untransformed) space. Will only be differences
#'   when init_bounded types are used in the ADMB template. This
#'   can be useful for model debugging.
#' @return If as.list is FALSE, an invisible data.frame containing samples
#'   (rows) of each parameter (columns). If multiple chains exist they will
#'   be rbinded together, maintaining order within each chain. If as.list
#'   is TRUE, samples are returned as a list of matrices.
#' @export
#' @examples
#' ## A previously run fitted ADMB model
#' fit <- readRDS(system.file('examples', 'fit.RDS', package='adnuts'))
#' post <- extract_samples(fit)
#' tail(apply(post, 2, median))
extract_samples <- function(fit, inc_warmup=FALSE, inc_lp=FALSE,
                            as.list=FALSE, unbounded=FALSE){
  if(!is.adfit(fit)) stop("fit is not a valid object")
  if(unbounded){
    x <- fit$samples_unbounded
    if(is.null(x))
      stop("No unbounded parameters in this fit")
  } else {
    x <- fit$samples
    if(is.null(x)) stop("No posterior samples found")
  }
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
#'   \code{\link{plot_sampler_params}} or  \code{\link{launch_shinyadmb}}.
#'
#' @param fit A list returned by \code{sample_admb}.
#' @param inc_warmup Whether to extract the warmup samples or not
#'   (default). Warmup samples should never be used for inference, but may
#'   be useful for diagnostics.
#' @return An invisible data.frame containing samples (rows) of each
#'   parameter (columns). If multiple chains exist they will be rbinded
#'   together.
#' @seealso \code{\link{launch_shinyadmb}}.
#' @export
#' @examples
#' fit <- readRDS(system.file('examples', 'fit.RDS', package='adnuts'))
#' sp <- extract_sampler_params(fit, inc_warmup=TRUE)
#' str(sp)
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
