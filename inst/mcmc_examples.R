## A very quick demonstration of the samplers. More to come.

## Run the simple example, so that obj and opt are loaded into workspace
library(TMB)
library(adnuts)
TMB::runExample("simple")
## Remake obj without random components
obj <- MakeADFun(data=list(x=x, B=B, A=A),
                 parameters=list(u=u*0, beta=beta*0, logsdu=1, logsd0=1),
                 random=NULL, DLL="simple", silent=TRUE)

## Note that this model does not use any box constraints. These can be
## passed to sample_tmb just like with nlminb. Also, there are no explicit
## priors in this model, which can be influential particularly for the
## logsd params.

## init can be a function, or a list of lists, or NULL to use starting
## values in obj.
init <- function() list(mu=u, beta=beta, logsdu=0, logsd0=0)
iter <- 2000
seeds <- 1:3

## First use a diagonal mass matrix (not recommended)
fit1 <- sample_tmb(obj=obj, iter=iter, chains=3, init=init, seeds=seeds,
                   control=list(adapt_mass=FALSE))
## Extra posterior samples like this
posterior <- extract_samples(fit1)
apply(posterior, 2, mean)
## Check diagnostics with shinystan:
launch_shinytmb(fit1)
## Or look at the values directly
sp <- extract_sampler_params(fit1)
str(sp)
min(fit1$ess)
max(fit1$Rhat)

## Can also use mass matrix adaptation (default -- diagonal only)
fit2 <- sample_tmb(obj=obj, iter=iter, chains=3, init=init, seeds=seeds,
                   control=list(adapt_mass=TRUE))

## Or pass an estimated one from a previous run (or could be MLE if it
## exists). This will help a lot of there are strong correlations in the
## model.
fit3 <- sample_tmb(obj=obj, iter=iter, chains=3, init=init, seeds=seeds,
                   control=list(metric=fit2$covar.est))

## Parallel works too
library(snowfall)
path <- system.file("examples", package = "TMB")
## DLL needs to be available in this folder since each node calls MakeADFun
## again to remake obj. So recompile here.
compile(file.path(path,'simple.cpp'))
fit4 <- sample_tmb(obj=obj, iter=iter, chains=3, seeds=1:3, init=init,
                   parallel=TRUE, cores=3, path=path)

library(vioplot)
## The mass matrix typically doesn't effect ESS since it runs long enough
## to get nearly independent samples. Instead each trajectory is shorter
## and hence the chains run faster
vioplot(fit1$ess, fit2$ess, fit3$ess, fit4$ess)

## Run time:
sum(fit1$time.total)
sum(fit2$time.total)
sum(fit3$time.total)
max(fit4$time.total)

## Efficiency:
min(fit1$ess)/sum(fit1$time.total)
min(fit2$ess)/sum(fit2$time.total)
min(fit3$ess)/sum(fit3$time.total)
min(fit4$ess)/max(fit4$time.total)


### ----------- ADMB example
### Note: ADMB functionality is still in flux so this will change a bit.
path <- 'mvn' # a toy multivariate normal model
## Compile and run model with hbf 1 (needed for NUTS). It is best to have
## your ADMB files in a separate folder and provide that path.
setwd(path)
## note: you need to compile this ADMB fork:
## https://github.com/colemonnahan/admb
system('admb mvn.tpl')
system('mvn -hbf 1')
setwd('..')
init <-  lapply(1:3, function(i) list(mu=rnorm(50)))
fit5 <- sample_admb(model='mvn', iter=2000, init=init, chains=3, path=path)
## Can also run parallel
library(snowfall)
fit6 <- sample_admb(model='mvn', iter=2000, init=init, chains=3, path=path,
                   parallel=TRUE, cores=3)
launch_shinyadmb(fit6)
fit6$cmd[1] # this is the command line arguments used

## Mass matrix adapation is almost done but not quite ready, check back
## soon!
