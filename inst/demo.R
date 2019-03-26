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
## If admb is not in the PATH this line will fail and you need to manually compile the model.
system('admb simple.tpl')
system('simple -nohess')  # do not need a .hes or .cov file to run!
setwd('..')

## 3 options to specify inits: list, NULL (uses MLE), or a function that
## returns a list.
init <-  lapply(1:3, function(i) rnorm(2))
init <- NULL # uses MLEs -- not recommended!
init <- function() rnorm(2)
fit <- sample_admb(model='simple', init=init, path=path)
fit$cmd[1] # this is the command line arguments used
## Merged chains after discarding warmup phase
post <- extract_samples(fit)
## Since no .cor file we are missing some nice information, like parameter
## names, but it still runs fine. Thus we can run for models without
## defined modes or uninvertible Hessian matrices. However, if available it
## is useful
str(post)
## A list with MLE stuff, after rerunning with Hessian
setwd(path); system('simple'); setwd('..')
fit <- sample_admb(model='simple', init=init, path=path)
str(fit$mle)
pairs_admb(fit) # modified pairs just for ADMB fits like this
## Can also use ShinyStan (make sure to exit it)
## launch_shinyadmb(fit)

## Can also run in parallel, including executing the mceval phase
## afterward. Note that the .psv file only has post-warmup and thinned (if
## used) samples in it. Chains are rbind'ed together, and thus the -mceval
## call runs on all chains and can be used directly for inference. You can
## set mceval=TRUE here, or run your model later manually as typically
## done.
library(snowfall)
cores <- parallel::detectCores()-1
fit <- sample_admb(model='simple', init=init, path=path,
                   parallel=TRUE, cores=cores, mceval=TRUE)

## The default is use an adaptive mass matrix. If we want to use the MLE
## covariance as the mass matrix, set it using control. This *should* make
## for more efficient sampling in most cases.  For technical reasons, must
## optimize the model with flag -hbf 1
setwd(path); system('simple -hbf 1'); setwd('..')
fit <- sample_admb(model='simple', init=init, path=path,
                    parallel=TRUE, cores=cores, control=list(metric='mle'))


## Can also use a slightly modified version of the original Metropolis MCMC
## algorithm (Random Walk Metropolis or RWM). Here we really want to use
## the MLE covariance. Since we ran above with hbf=1 we need to rerun model
## to recreate the .hes and .cov files.
setwd(path); system('simple'); setwd('..')
fit <- sample_admb(model='simple', init=init, path=path, algorithm='RWM',
                    iter=200000, thin=100, mceval=TRUE,
                    parallel=T, cores=cores, control=list(metric='mle'))

## Can also specify a duration argument for capping the run time. This is
## useful e.g. if you want to run overnight but have results by 8am.  Here
## we do 0.5 minutes just to demonstrate. This period needs to be long
## enough to do the warmup or it'll throw an error. Normally we wouldn't
## use thinning for NUTS but only to demonstrate here since the model runs
## so fast.
fit <- sample_admb(model='simple', init=init, path=path, algorithm='NUTS',
                    warmup=200, iter=20000000, thin=100, duration=.5,
                    parallel=TRUE, cores=cores)
## Notice that the variable chain lengths are truncated to the minimum, and
## then merged together.
str(extract_samples(fit))
str(extract_samples(fit, as.list=TRUE))

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
