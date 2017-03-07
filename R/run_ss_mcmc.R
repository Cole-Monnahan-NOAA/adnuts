
#' A wrapper for running SS models in parallel
#' @export
run_ss_mcmc_parallel <- function(parallel_number,  Nout, mcsave, burn.in, cov.user = NULL,
                                 init.pin = NULL, se.scale = 1, mcscale = TRUE, mcseed = NULL,
                                 mcrb = NULL, mcdiag = FALSE, mcprobe = NULL, verbose = TRUE,
                                 extra.args = NULL, hybrid = FALSE, hyeps = NULL, hynstep = NULL,
                                 mceval = FALSE, estimate = FALSE){
    olddir <- getwd()
    on.exit(setwd(olddir))
    newdir <- paste0(getwd(),'/model',parallel_number)
    if(dir.exists(newdir)) unlink(newdir, TRUE)
    dir.create(newdir)
    trash <- file.copy(from=list.files('model', full.names=TRUE), to=newdir)
    ## delay in case indexing ties up files briefly
    Sys.sleep(runif(1, .5, 1))
    time <- system.time(xx <-
        admbtools::run_ss_mcmc(model.path=newdir, model.name='model',
                               Nout=Nout, mcsave=mcsave, burn.in=burn.in,
                               cov.user=cov.user, init.pin=init.pin,
                               se.scale=se.scale, mcscale=mcscale,
                               mcseed=mcseed, mcrb=mcrb, mcdiag=mcdiag,
                               mcprobe=mcprobe, verbose=verbose,
                               extra.args=extra.args, hybrid=hybrid,
                               hyeps=hyeps, hynstep=hynstep, mceval=mceval,
                               estimate=estimate) )
    cbind(parallel_number=parallel_number, iteration=1:nrow(xx), xx)
}
