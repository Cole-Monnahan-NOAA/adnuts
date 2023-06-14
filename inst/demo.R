### A very quick demonstration of no-U-turn sampling in ADMB
library(adnuts)

### ----------- ADMB example
## This is the packaged simple regression model
path.simple <- system.file('examples', 'simple', package='adnuts')
## It is best to have your ADMB files in a separate folder and provide that
## path, so make a copy of the model folder locally.
path <- 'simple'
dir.create(path)
trash <- file.copy(from=list.files(path.simple, full.names=TRUE), to=path)

## Compile and run model
setwd(path)
## If admb is not in the PATH this line will fail and you need to manually compile the model.
system('admb simple.tpl')
system('simple')
setwd('..')

## 3 options to specify inits: list, NULL (uses MLE), or a function that
## returns a list.
npar <- 2 # number of parameters in the model
init <-  lapply(1:3, function(i) rnorm(npar))
init <- NULL # uses MLEs -- not recommended!
init <- function() rnorm(npar)
fit <- sample_nuts(model='simple', init=init, path=path, cores=1)
fit$cmd[1] # this is the command line arguments used
## Merged chains after discarding warmup phase
post <- extract_samples(fit)
str(post)
## A list with MLE fit
str(fit$mle)

## Can also run in parallel which is the default when
## chains>1. Here we also execute the mceval phase
## afterward. Note that the .psv file only has post-warmup and
## thinned (if used) samples in it. Chains are rbind'ed together,
## and thus the -mceval call runs on all chains and can be used
## directly for inference. You can set mceval=TRUE here, or run
## your model later manually as typically done.
fit <- sample_nuts(model='simple', init=init, path=path, mceval=TRUE)


## Can also specify a duration argument for capping the run
## time. This is useful e.g. if you want to run overnight but
## have results by 8am.  Here we do 0.5 minutes just to
## demonstrate. This period needs to be long enough to do the
## warmup or it'll throw an error. Normally we wouldn't use
## thinning for NUTS but only to demonstrate here since the model
## runs so fast. I recommended setting a reasonable warmup, then
## an unreasonably big iter. It will truncate to the shortest
## chains.
fit <- sample_nuts(model='simple', init=init, path=path,
                   warmup=200, iter=20000000, thin=100,
                   duration=.5)

## The default is use an adaptive mass matrix. If we want to use
## the MLE covariance as the mass matrix, set it using
## control. This *should* make for more efficient sampling in
## most cases.  For technical reasons, must optimize the model
## with flag -hbf 1. Here setting skip_optimization=FALSE will
## rerun the model before starting the chains. Or you could
## manually optimize the model with the flag '-hbf'.
fit <- sample_nuts(model='simple', init=init, path=path,
                   skip_optimization=FALSE,
                   control=list(metric='mle'))
## See the vignette for a discussion on the options for the
## metric (mass matrix adaptation and such).

## Can also use a slightly modified version of the original
## Metropolis MCMC algorithm (Random Walk Metropolis or
## RWM). Here we really want to use the MLE covariance. Since we
## ran above with hbf=1 we need to rerun model to recreate the
## .hes and .cov files. Can do this with skip_optimization
## argument again.
fit.rwm <- sample_rwm(model='simple', init=init, path=path,
                  skip_optimization=FALSE,
                  iter=200000, thin=100, mceval=TRUE)

### ------------------------------------------------------------
### Convergence diagnostics
## Key convergence diagnostics (effective sample size and Rhat)
## are available in the summary of the fit. It is your
## responsibility to check for signs of non-convergence before
## using the output for inference!!
summary(fit)
## or from the rstan::monitor output
str(fit$monitor)

### ------------------------------------------------------------
### Plotting options
pairs_admb(fit) # modified pairs just for ADMB fits like this
## Can also use ShinyStan (make sure to exit it)
## launch_shinyadmb(fit)
plot_sampler_params(fit)                # NUTS adaptation
## Compare MLE and posterior of marginals. See help for
## recommendation for creating multipage PDF for high dimensional
## parameters.
plot_marginals(fit)

### ------------------------------------------------------------
### Extracting posterior samples
## Get post-warmup, merged chains into data frame (these two are
## identical)
str(as.data.frame(fit))
str(extract_samples(fit))
## If you want it in list form, e.g., to put into coda package
str(extract_samples(fit, as.list=TRUE))
## If you want to see warmup samples, and the log-posterior (lp__ column)
str(extract_samples(fit, inc_warmup=TRUE, inc_lp=TRUE))

## Remove folder
unlink(path, TRUE)


### End of demo


