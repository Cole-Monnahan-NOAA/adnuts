## Draw samples with the RWM, HMC and NUTS algorithms and compare.

## Run the simple example, so that obj and opt are loaded into workspace
library(TMB)
library(adnuts)
runExample("simple")
## Remake obj without random components
obj <- MakeADFun(data=list(x=x, B=B, A=A),
                 parameters=list(u=u*0, beta=beta*0, logsdu=1, logsd0=1),
                 random=NULL, DLL="simple", silent=TRUE)
str(obj)

## Run RWM and two gradient based algorithms, using adative step size (eps)
## for each. Start from the MLE.
rwm <- run_mcmc(obj=obj, iter=20000, thin=10, algorithm='RWM', alpha=.05, chains=3)
rwm.samples <- extract_samples(rwm)
rwm.diagnostics <- rstan::monitor(sims=rwm$samples)
## Calculate estimated samples per time (efficiency)
min(rwm.diagnostics[,'n_eff'])/rwm$time.total

nuts <- run_mcmc(obj=obj, iter=2000, algorithm='NUTS', chains=3)
nuts.samples <- extract_samples(nuts)
nuts.diagnostics <- rstan::monitor(sims=nuts$samples)
## Calculate estimated samples per time (efficiency)
min(nuts.diagnostics[,'n_eff'])/nuts$time.total

## Can also use shinystan to look deeper at diagnostics
if(require(shinystan)){
  launch_shinystan_tmb(rwm)
  launch_shinystan_tmb(nuts)
}


