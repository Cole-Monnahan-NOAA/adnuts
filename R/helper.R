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
                "instead of ",NROW(cov.user)))
  ## Write it to file using original scales, although these are ignored.
  setwd(wd.old)
  file.new <- file(paste0(model.path, "/admodel.cov"),"wb")
  on.exit(close(file.new))
  writeBin(as.integer(num.pars), con=file.new)
  writeBin(as.vector(as.numeric(cov.unbounded)), con=file.new)
  writeBin(as.integer(results$hybrid_bounded_flag), con=file.new)
  writeBin(as.vector(scale), con=file.new)
}
