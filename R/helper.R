write.admb.cov.unbounded <- function(cov.unbounded, model.path=getwd()){
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
