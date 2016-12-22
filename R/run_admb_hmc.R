#' Run an MCMC using an ADMB model, return (1) the posterior
#' draws, MLE fits and covariance/correlation matrices, and some
#' MCMC convergence diagnostics using CODA.
#'
#' @param model.path (Character) A path to the folder containing the model. NULL
#' indicates the current folder.
#' @param mode.name (Character) The name of the model executable. A character string,
#' without '.exe'.
#' @param iter (Integer) The number of draws after thinning and burn in.
#' @param thin (Integer) Controls thinning of samples. Save every thin
#' value, such that 1 corresponds to keeping all draws, and 100 saving
#' every 100th draw.
#' @param warmup (Integer) How many samples to discard from the beginning
#' of the chain, *after* thining. The burn in period (i.e., the first
#' warmup*thin draws) should be at least large enough to cover dynamic
#' scaling.
#' @param cov.user (Numeric matrix) A manually defined covariance matrix (in bounded space)
#' to use in the Metropolis-Hastings algorithm.
#' @param init (Numeric vector) A vector of initial values, which are written to file
#' and used in the model via the -mcpin option.
#' @param mcseed (Integer) Which seed (integer value) to pass ADMB. Used
#' for reproducibility.
#' @param mcdiag (Logical) Whether to use the \code{mcdiag} feature. This
#' uses an identity matrix for the covariance matrix.
#' @param eps (Numeric) The size of the leapfrog jump in the hybrid
#' method, with smaller values leading to smaller but more accurate
#' jumps. Must be a positive value.
#' @param verbose (Logical) Whether to print ADMB warnings and other
#' information. Useful for testing and troubleshooting.
#' @param extra.args (Character) A string which is passed to ADMB at
#' runtime. Useful for passing additional arguments to the model
#' executable.
#' @export
#' @return Returns a list containing (1) the posterior draws, (2) and
#' object of class 'admb', read in using the results read in using
#' \code{read_admb}, and (3) some MCMC convergence diagnostics using CODA.
run_admb_nuts <-
  function(model.path, model.name, iter=2000, thin=1, warmup=ceiling(iter/2),
           init=NULL, eps=NULL, cov.user=NULL,  mcseed=NULL,
           mcdiag=FALSE, verbose=TRUE, extra.args=NULL,
           mceval=TRUE){
  wd.old <- getwd(); on.exit(setwd(wd.old))
  setwd(model.path)
  ## Grab original admb fit and metrics
  if(iter <1)
    stop("Iterations must be >1")
  if(!file.exists(paste0(model.name,'.par'))) {
    system(paste("admb", model.name))
    system(paste(model.name))
    }
  mle <- R2admb::read_admb(model.name)
  ## If user provided covar matrix, write it to file and save to
  ## results
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
  ## Write the starting values to file. Always using a
  ## init file b/c need to use -nohess -noest so
  ## that the cov.user can be specified and not
  ## overwritten. HOwever, this feature then starts the
  ## mcmc chain from the initial values instead of the
  ## MLEs. So let the user specify the init values, or
  ## specify the MLEs manually
  if(is.null(init))
    init <- mle$coefficients[1:mle$npar]
  write.table(file="init.pin", x=init, row.names=F, col.names=F)
  ## Separate the options by algorithm, first doing the shared
  ## arguments
  cmd <- paste(model.name,"-noest -nohess -nuts -mcmc ",iter)
  cmd <- paste(cmd, "-mcpin init.pin")
  if(!is.null(extra.args))
    cmd <- paste(cmd, extra.args)
  if(!is.null(mcseed))
    cmd <- paste(cmd, "-mcseed", mcseed)
  if(!is.null(eps)){
    cmd <- paste(cmd, "-hyeps", eps)
    }
  ## Run it and get results
  system(cmd, ignore.stdout=!verbose)
  if(mceval) system(paste(model.name, "-mceval -noest -nohess"),
                    ignore.stdout=!verbose)
  psv <- file(paste0(model.name, ".psv"), "rb")
  adapt <- as.matrix(read.csv("adaptation.csv"))
  pars <- read_psv(model.name)
  pars[,'log-posterior'] <- adapt[,'energy__']
  pars2 <- array(0, dim=c(nrow(pars), 1, ncol(pars)))
  pars2[,1,] <- as.matrix(pars)
  dimnames(pars2) <-
    list(iter=1:nrow(pars), chains="chain:1",
         parameters=dimnames(pars)[[2]])
  ss <- monitor(sims=pars2)
  y <- vector("list", length=length(dimnames(pars2)[[3]]))
  names(y) <- dimnames(pars2)[[3]]
  z <- lapply(y, function(x) x=numeric(0))
  sso <-
    shinystan:::shinystan(
    model_name=model.name, param_names=names(pars), param_dims=z,
    posterior_sample=pars2, sampler_params=list(adapt),
    summary=ss, n_chain=1, n_iter=nrow(pars),
    n_warmup=nrow(pars)/2, model_code='NA',
    misc=list(max_td=12, stan_method='sampling',
              stan_algorithm='NUTS',
              sso_version=utils::packageVersion('shinystan')))
  return(sso)
  }
