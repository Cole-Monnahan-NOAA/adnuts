## Note these were tested with ADMB 12.2 on 8/3/2020
test_that("class methods working for ADMB", {
  fit <- readRDS('fit_admb.RDS')
  x <- as.data.frame(fit)
  expect_is(x, 'data.frame')
  summary(fit)
  print(fit)
  plot(fit)
  plot_marginals(fit, pars=1:2)
  plot_uncertainties(fit)
  plot_sampler_params(fit)
  pairs(fit, pars=1:2)
  pairs(fit, pars=1:2, order='slow')
  fit$mle <- NULL
  pairs(fit, pars=1:2)
  pairs(fit, pars=1:2, order='slow')
})


test_that("class methods working for TMB", {
  fit <- readRDS('fit_snuts.RDS')
  x <- as.data.frame(fit)
  expect_is(x, 'data.frame')
  summary(fit)
  print(fit)
  plot(fit)
  plot_marginals(fit, pars=1:2)
  plot_uncertainties(fit)
  plot_sampler_params(fit)
  pairs(fit, pars=1:2)
  pairs(fit, pars=1:2, order='slow')
  fit$mle <- NULL
  pairs(fit, pars=1:2)
  pairs(fit, pars=1:2, order='slow')
})
