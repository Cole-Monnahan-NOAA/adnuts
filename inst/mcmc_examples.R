### A very quick demonstration of no-U-turn sampling in ADMB and TMB.

## Note that this model does not use any box constraints. These can be
## passed to sample_tmb just like with nlminb. Also, there are no explicit
## priors in this model, which can be influential particularly for the
## logsd params. The user is solely responsible for properly defining a
## Bayesian model, including priors and bounds.
library(TMB)
library(adnuts)
TMB::runExample("simple")

## init can be a function, or a list of lists, or NULL to use starting
## values in obj.
init <- function() list(mu=u, beta=beta, logsdu=0, logsd0=0)
seeds <- 1:3
## The default is to run 3 chains, 50% warmup, and use diagonal mass matrix
## adapation. This is the same as rtan.
fit1 <- sample_tmb(obj=obj, init=init, seeds=seeds)

## Can also run in parallel, but the DLL needs to be available in this
## folder since each node calls MakeADFun again to remake obj.  This
## usually is your working directory.
path <- system.file("examples", package = "TMB")
## This line is only needed b/c of how TMB::runexample works above.
compile(file.path(path,'simple.cpp'))
fit2 <- sample_tmb(obj=obj, seeds=seeds, init=init, iter=100,
                   parallel=TRUE, cores=3, path=path)

## Extract samples like this
post1 <- extract_samples(fit1)
## You can use these however you want. For instance, loop through each
## saved row and call object$report(post1[i,]) to do any calculations or
## projections (equivalent to -mceval in ADMB).

## Run time, one for each chain. In parallel constrained by longest one.
fit1$time.total
sum(fit1$time.total)
max(fit2$time.total)

## Use rstan functions to calculate effective sample sizes (ESS) and Rhat
## (potential scale reduction.
m1 <- rstan::monitor(fit1$samples, print=FALSE)
Rhat <- m1[,"Rhat"]
max(Rhat)
ess <- m1[, 'n_eff']
min(ess)

## Or explore with shinystan. Make sure to click "Save and close" to exit
## properly.
launch_shinytmb(fit1)

### ----------- ADMB example
### Note: ADMB functionality is still in flux so this will change a bit.
path <- 'simple' # a toy multivariate normal model
## Compile and run model with hbf 1 (needed for NUTS). It is best to have
## your ADMB files in a separate folder and provide that path.
setwd(path)
## note: you need to compile this ADMB fork:
## https://github.com/colemonnahan/admb
system('admb simple.tpl')
system('simple -hbf 1')
setwd('..')
init <-  lapply(1:3, function(i) rnorm(2))
fit3 <- sample_admb(model='simple', init=init, path=path)
## Can also run parallel
library(snowfall)
fit4 <- sample_admb(model='simple', init=init, path=path,
                   parallel=TRUE, cores=3)
launch_shinyadmb(fit4)
fit4$cmd[1] # this is the command line arguments used

