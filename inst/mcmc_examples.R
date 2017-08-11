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


