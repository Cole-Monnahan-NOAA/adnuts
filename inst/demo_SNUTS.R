#devtools::load_all()
library(adnuts)
library(TMB)

TMB::runExample('simple')
fit <- sample_snuts(obj)
pairs(fit)
pairs(fit, order='slow')
pairs(fit, order='mismatch')
pairs(fit, order='fast')
print(fit)

# quick plots to explore output
plot_sampler_params(fit)
plot_uncertainties(fit)
plot_marginals(fit, pars=1:9)
pairs_admb(fit, pars=1:5, order='slow')
# launch_shiny(fit)

## It also works for RTMB. Example from the repo, simplified down
## a bit
if('TMB' %in% .packages()) detach(package:TMB)
library(RTMB)
data(ChickWeight)
parameters <- list(mua=0, sda=1, mub=0, sdb=1, sdeps=1,
                   a=rep(0, 50), b=rep(0, 50))
f <- function(parms) {
  RTMB::getAll(ChickWeight, parms, warn=FALSE)
  ## Random slopes
  nll <-  - sum(RTMB::dnorm(a, mean=mua, sd=sda, log=TRUE))
  ## Random intercepts
  nll <- nll - sum(RTMB::dnorm(b, mean=mub, sd=sdb, log=TRUE))
  ## Data
  predWeight <- a[Chick] * Time + b[Chick]
  nll <- nll - sum(RTMB::dnorm(weight, predWeight, sd=sdeps, log=TRUE))
  nll
}
obj <- RTMB::MakeADFun(f, parameters, random=c("a", "b"))
fit <- sample_snuts(obj)
