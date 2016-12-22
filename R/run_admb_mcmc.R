#' Run an MCMC using an ADMB model, return (1) the posterior draws, MLE
#' fits and covariance/correlation matrices, and some MCMC convergence
#' diagnostics using CODA.
#'
#' @param model.path (Character) A path to the folder containing the model. NULL
#' indicates the current folder.
#' @param mode.name (Character) The name of the model executable. A character string,
#' without '.exe'.
#' @param iter (Integer) The number of draws after thinning and burn in.
#' @param mcsave (Integer) Controls thinning of samples. Save every mcsave
#' value, such that 1 corresponds to keeping all draws, and 100 saving
#' every 100th draw.
#' @param burn.in (Integer) How many samples to discard from the beginning
#' of the chain, *after* thining. The burn in period (i.e., the first
#' burn.in*mcsave draws) should be at least large enough to cover dynamic
#' scaling.
#' @param cov.user (Numeric matrix) A manually defined covariance matrix (in bounded space)
#' to use in the Metropolis-Hastings algorithm.
#' @param init.pin (Numeric vector) A vector of initial values, which are written to file
#' and used in the model via the -mcpin option.
#' @param se.scale (Numeric) A value which scales all of the variances from
#' the MLE fit. A value of 1 indicates to use the estimated variances.
#' @param mcscale (Logical) Whether to use the mcscale option, which
#' dynamically scales the covariance matrix for efficient acceptance
#' ratios.
#' @param mcseed (Integer) Which seed (integer value) to pass ADMB. Used
#' for reproducibility.
#' @param mcrb (Integer) Which value to use in the rescale bounded
#' algorithm. Must be an integer from 1-9. The default NULL value disables
#' this feature. See the vignette for more information on this algorithm
#' and how to best use it.
#' @param mcdiag (Logical) Whether to use the \code{mcdiag} feature. This
#' uses an identity matrix for the covariance matrix.  #' @param mcprobe
#' Which value to use in the probing algorithm. The default NULL value
#' disables this feature. See the vignette for more information on this
#' algorithm and how to best use it.
#' @param hyrbid (Logical) Whether to use the Hamiltonial (hybrid)
#' algorithm. Default is FALSE.
#' @param hyeps (Numeric) The size of the leapfrog jump in the hybrid
#' method, with smaller values leading to smaller but more accurate
#' jumps. Must be a positive value.
#' @param hynstep (Integer) The approximate number of steps used in the
#' leapfrog step of the hybrid algorithm. Steps are randomly generated for
#' each MCMC iteration, centered around \code{hynstep}.
#' @param verbose (Logical) Whether to print ADMB warnings and other
#' information. Useful for testing and troubleshooting.
#' @param extra.args (Character) A string which is passed to ADMB at
#' runtime. Useful for passing additional arguments to the model
#' executable.
#' @export
#' @return Returns a list containing (1) the posterior draws, (2) and
#' object of class 'admb', read in using the results read in using
#' \code{read_admb}, and (3) some MCMC convergence diagnostics using CODA.
run_admb_mcmc <- function(model.path, model.name, iter, mcsave, burn.in,
                          cov.user=NULL, init.pin=NULL, se.scale=NULL,
                          mcscale=FALSE,  mcseed=NULL, mcrb=NULL, mcdiag=FALSE,
                          mcprobe=NULL, verbose=TRUE, extra.args=NULL,
                          hybrid=FALSE, hyeps=NULL, hynstep=NULL,
                          mceval=FALSE, estimate=FALSE){
    ## This function runs an ADMB model MCMC, burns and thins, calculates
    ## effective sizes, and returns stuff depending on verbose.
    ## browser()
    iterations <- (iter+burn.in)*mcsave
    if(iterations <1) stop(paste0("Iterations too low: ", iterations))
    if(verbose) print(paste("Run started at", round(Sys.time())))
    wd.old <- getwd(); on.exit(setwd(wd.old))
    setwd(model.path)
    ## Run to get MLE and covariance matrix
    if(estimate)
        system(model.name, ignore.stdout=T)
    ## Grab original admb fit and metrics
    mle <- read_admb(model.name)
    ## If user provided covar matrix, write it to file and save to results
    if(!is.null(cov.user)){
        cor.user <- cov.user/ sqrt(diag(cov.user) %o% diag(cov.user))
        if(!is.positive.definite(x=cor.user))
            stop("Invalid cov.user matrix, not positive definite")
        write.admb.cov(cov.user)
        mle$cov.user <- cov.user
    } else {
    ## otherwise use the estimated one
        mle$cov.user <-  NULL
    }
    ## Write the starting values to file. Always using a init.pin file b/c
    ## need to use -nohess -noest so that the cov.user can be specified and
    ## not overwritten. HOwever, this feature then starts the mcmc chain
    ## from the initial values instead of the MLEs. So let the user specify
    ## the init values, or specify the MLEs manually
    if(is.null(init.pin)) init.pin <- mle$coefficients[1:mle$npar]
    write.table(file="init.pin", x=init.pin, row.names=F, col.names=F)
    ## Separate the options by algorithm, first doing the shared arguments
    cmd <- paste(model.name,"-mcmc",iterations)
    ## If user written one, make sure not to overwrite it
    if(!is.null(cov.user)) cmd <- paste(cmd, "-nohess")
    cmd <- paste(cmd, "-mcpin init.pin")
    if(!is.null(extra.args)) cmd <- paste(cmd, extra.args)
    if(!is.null(mcseed)) cmd <- paste(cmd, "-mcseed", mcseed)
    if(mcdiag==TRUE) cmd <- paste(cmd, "-mcdiag")
    if(!is.null(mcrb)) cmd <- paste(cmd, "-mcrb",mcrb)
    ## Those options for the standard MH algorithm
    if(!hybrid){
        cmd <- paste(cmd, "-mcsave",mcsave)
        if(mcscale==FALSE) cmd <- paste(cmd, "-mcnoscale")
        if(!is.null(mcprobe)) cmd <- paste(cmd, "-mcprobe",mcprobe)
    } else {
        ## The hybrid options
        if(mcsave!=1)
            stop("mcsave option is incompatible with the hybrid algorithm (fixed at 1 internally)")
        if(hyeps<=0) stop("hyeps must be positive number")
        if(hynstep<=0) stop("hynstep must be positive integer")
        cmd <- paste(cmd, "-hybrid -hyeps",hyeps,"-hynstep",hynstep)
    }
    ## The command is constructed.
    if(verbose) print(cmd)
    ## Scale the covariance matrix
    if(!is.null(se.scale)) SetScale(se.scale) # change the covariance matrix
    ## Run it and get results
    system(cmd, ignore.stdout=T)
    if(mceval)
        system(paste(model.name, "-mceval -noest -nohess"), ignore.stdout=T)
    psv <- file(paste0(model.name, ".psv"), "rb")
    nparams <- readBin(psv, "integer", n=1)
    mcmc <- matrix(readBin(psv, "numeric", n=nparams*(iter+burn.in)), ncol=nparams,
                   byrow=TRUE)
    close(psv)
    mcmc <- as.data.frame(mcmc)
    if(!is.null(mle)){
        names(mcmc) <- names(with(mle, coefficients[1:npar]))
    }
    ## If mcrb is used, read  that in for plotting. It's in the corrtest file.
    if(!is.null(mcrb)){
        L <- readLines('corrtest')
        st <- grep("modified S", L)[1]+1
        en <- grep("S* modified S", L)-1
        L <- gsub("^\\s+|\\s+$", "", L[st:en])
        cov.mcrb <- do.call(rbind, lapply(strsplit(L, " "), as.numeric))
        cor.mcrb <- cov.mcrb/(sqrt(diag(cov.mcrb) %o% diag(cov.mcrb)))
        if(!is.positive.definite(cor.mcrb))
            warning("the modified mcrb covariance matrix read in was not positive definite")
        else mle$cov.user <- cov.mcrb
    }
    ## Remove the 'burn in' specified by user
    if(burn.in>0) mcmc <- mcmc[-(1:burn.in),]
    ## Run effective sample size calcs from CODA, metric of convergence
    efsize <- data.frame(t(effectiveSize(mcmc)/NROW(mcmc)))
    names(efsize) <- paste0(names(mcmc), "_efs")
    results <- list(mcmc=mcmc, mle=mle)
    results$diag <- list(efsize=efsize)
    class(results) <- 'admb_mcmc'
    return(results)
}


