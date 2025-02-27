#devtools::load_all()
library(adnuts)
library(TMB)

TMB::runExample('simple')
fit <- sample_sparse_tmb(obj)
print(fit)

# quick plots to explore output
plot_sampler_params(fit)
plot_uncertainties(fit)
plot_marginals(fit, pars=1:9)
pairs_admb(fit, pars=1:5, order='slow')

# launch_shinyadmb(fit)
