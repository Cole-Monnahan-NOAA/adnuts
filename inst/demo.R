### A very quick demonstration of no-U-turn sampling in ADMB and TMB.
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
system('admb simple.tpl')
system('simple')
setwd('..')

## 3 options to specify inits: list, NULL (uses MLE), or a function that
## returns a list.
init <-  lapply(1:3, function(i) rnorm(2))
init <- NULL
init <- function() rnorm(2)
fit <- sample_admb(model='simple', init=init, path=path)
fit$cmd[1] # this is the command line arguments used
## Merged chains after discarding warmup phase
post <- extract_samples(fit)
str(post)
## A list with MLE stuff
str(fit$mle)
pairs_admb(fit) # modified pairs just for ADMB fits like this
## Can also use ShinyStan (make sure to exit it)
## launch_shinyadmb(fit)

## Can also run in parallel, including executing the mceval phase
## afterward. Note that the .psv file only has post-warmup and thinned (if
## used) samples in it. Chains are rbind'ed together.
library(snowfall)
cores <- parallel::detectCores()-1
fit <- sample_admb(model='simple', init=init, path=path,
                   parallel=TRUE, cores=cores, mceval=TRUE)

## If we want to use the MLE covariance as the metric set it using
## control. This *should* make for more efficient sampling in most cases.
## For technical reasons, must optimize the model with flag -hbf 1
setwd(path); system('simple -hbf 1'); setwd('..')
fit <- sample_admb(model='simple', init=init, path=path,
                    parallel=TRUE, cores=cores, control=list(metric='mle'))


## Can also use the RWM. Here we really want to use the MLE covariance.
fit <- sample_admb(model='simple', init=init, path=path, algorithm='RWM',
                    iter=200000, thin=100, mceval=TRUE,
                    parallel=TRUE, cores=cores, control=list(metric='mle'))

## Can also specify a duration argument for capping the run time. Here we
## do 0.5 minutes. This period needs to be long enough to do the warmup or
## it'll throw an error.
fit <- sample_admb(model='simple', init=init, path=path, algorithm='NUTS',
                    warmup=100, iter=20000000, thin=100, duration=.5,
                    parallel=TRUE, cores=cores)
str(extract_samples(fit))

## Remove folder
unlink(path, TRUE)
### End of ADMB example



### --------------------------------------------------
### TMB example: Note, this is large replaced by tmbstan. See that package
### for more info on using Stan with TMB models.

## Note that this model does not use any box constraints. These can be
## passed to sample_tmb just like with nlminb. Also, there are no explicit
## priors in this model, which can be influential particularly for the
## logsd params. The user is solely responsible for properly defining a
## Bayesian model, including priors and bounds.
library(TMB)
TMB::runExample("simple")

## init can be a function, or a list of lists, or NULL to use starting
## values in obj.
init <- function() list(mu=u, beta=beta, logsdu=0, logsd0=0)
seeds <- 1:3
## The default is to run 3 chains, 50% warmup, and use diagonal mass matrix
## adapation. This is the same as rtan.
fit <- sample_tmb(obj=obj, init=init, seeds=seeds)

## Can also run in parallel, but the DLL needs to be available in this
## folder since each node calls MakeADFun again to remake obj.  This
## usually is your working directory.
path <- system.file("examples", package = "TMB")
## This line is only needed b/c of how TMB::runexample works above.
compile(file.path(path,'simple.cpp'))
library(snowfall)
fit <- sample_tmb(obj=obj, seeds=seeds, init=init,
                   parallel=TRUE, cores=cores, path=path)

## Extract samples like this
post <- extract_samples(fit)
## You can use these however you want. For instance, loop through each
## saved row and call object$report(post1[i,]) to do any calculations or
## projections (equivalent to -mceval in ADMB).

## Use rstan functions to calculate effective sample sizes (ESS) and Rhat
## (potential scale reduction.
mon <- rstan::monitor(fit$samples, print=FALSE)
Rhat <- mon[,"Rhat"]
max(Rhat)
ess <- mon[, 'n_eff']
min(ess)

## Or explore with shinystan. Make sure to click "Save and close" to exit
## properly.
launch_shinytmb(fit)
