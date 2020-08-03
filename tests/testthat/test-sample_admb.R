test_that("simple example works", {
  skip_on_cran()
  inits.fn <- function() list(c(0,0))
  oldwd <- getwd()
  fit <- sample_rwm('simple', path='../simple', chains=1,
                     seeds=1, init=inits.fn,
                     skip_optimization=FALSE,
                     control=list(refresh=0), skip_monitor=TRUE)
  expect_known_output(extract_samples(fit)[1000,],
                      file='_expect_simple_rwm')
  fit <- sample_nuts('simple', path='../simple', chains=1,
                     seeds=1, init=inits.fn,
                     skip_optimization=FALSE,
                     control=list(refresh=0), skip_monitor=TRUE)
  expect_known_output(extract_samples(fit)[1000,],
                      file='_expect_nuts')
  fit <- sample_nuts('simple', path='../simple', chains=1,
                     seeds=1, init=inits.fn,
                     control=list(metric='mle', refresh=0),
                     skip_monitor = TRUE)
  expect_known_output(extract_samples(fit)[1000,],
                      file='_expect_nuts_mle')
})

test_that("reproducibility of algorithms", {
  skip_on_cran()
  ## Check reproducibility given same init and seeds
  inits.fn <- function() list(c(0,0))
  oldwd <- getwd()
  chains <- 10
  fit2 <- sample_rwm('simple', path='../simple', chains=chains,
                     iter=400, cores=1,
                     seeds=rep(45,chains), init=inits.fn,
                     skip_optimization=FALSE,
                     admb_args=' -refresh 0')
  expect_identical(unique(fit2$samples[400,,3]), -16.0439)
  fit1 <- sample_nuts('simple', path='../simple', chains=chains, iter=400,
                      seeds=rep(123,chains), init=inits.fn,
                      skip_optimization=FALSE,
                     cores=1, admb_args=' -refresh 0')
  expect_identical(unique(fit1$samples[400,,3]), -13.2107)
  fit3 <- sample_nuts('simple', path='../simple', chains=chains, iter=400,
                      seeds=rep(123,chains), init=inits.fn, cores=1,
                      control=list(metric='mle'), admb_args=' -refresh 0')
  expect_identical(unique(fit3$samples[400,,3]), -12.2534)
  fit4 <- suppressWarnings(sample_nuts('simple', path='../simple', chains=chains, iter=400,
                      seeds=rep(123,chains), init=inits.fn, cores=1,
                      control=list(metric=diag(2)), admb_args=' -refresh 0'))
  expect_identical(unique(fit4$samples[400,,3]), -13.2107)
  ## Should match fits 1 and 4 since same settings
  fit5 <- sample_nuts('simple', path='../simple', chains=chains, iter=400,
                      seeds=rep(123,chains), init=inits.fn, cores=1,
                      control=list(metric='unit'), admb_args=' -refresh 0')
  expect_identical(unique(fit5$samples[400,,3]), -13.2107)
  ## Uses no adaptation so should be consistent between ADMB versions whereas
  ## the other ones would not be.
  fit6 <- sample_nuts('simple', path='../simple', chains=chains, iter=400,
                      seeds=rep(123,chains), init=inits.fn, cores=1,
                      control=list(metric='unit', stepsize=.1, adapt_mass=FALSE),
                      admb_args=' -refresh 0')
  expect_identical(unique(fit6$samples[400,,3]), -11.6495)
})


test_that("diagnostics and plotting", {
  skip_on_cran()
  inits <- function() list(1,1)
  fit <- sample_nuts('simple', path='../simple', chains=4,
                     iter=2000, cores=1, init=inits, seeds=1,
                     control=list(refresh=-1))
  sp <- extract_sampler_params(fit)
  expect_known_output(tail(sp,1), file='_expect_sp')
  expect_known_output(tail(fit$monitor,1), file='_expect_monitor')
  plot_sampler_params(fit, TRUE)
  pairs_admb(fit)
  pairs_admb(fit, pars=1:3, order='slow')
  pairs_admb(fit, pars=1:3, order='fast')
  pairs_admb(fit, pars=c('a', 'lp__', 'b'), add.monitor=FALSE)
  pairs_admb(fit, add.mle=FALSE)
  pairs_admb(fit, add.mle=FALSE, diag='hist')
  pairs_admb(fit, add.mle=FALSE, diag='acf')
  plot_marginals(fit)
  plot_marginals(fit, add.monitor=FALSE)
  plot_marginals(fit, add.mle=FALSE)
})

test_that("warnings and errors in sample_nuts",{
  skip_on_cran()
  inits <- function() list(1,1)
  test <- expect_warning(sample_nuts('simple', path='../simple',
                                     iter=1000, init=inits,
                                     parallel=TRUE,
                                     chains=1, warmup=500),
                         regexp='parallel is deprecated')
  test <- expect_warning(sample_rwm('simple', path='../simple',
                                     iter=1000, init=inits,
                                     parallel=TRUE,
                                     chains=1, warmup=500),
                         regexp='parallel is deprecated')
  test <- expect_warning(sample_admb('simple', path='../simple',
                                    iter=1000, init=inits,
                                    parallel=TRUE,
                                    chains=1, warmup=500),
                         regexp='parallel is deprecated')
  test <- expect_warning(sample_admb('simple', path='../simple',
                                   iter=1000, init=inits,
                                   chains=1, warmup=500),
                       regexp='sample_admb is deprecated')
  test <- expect_error(sample_nuts('simple', path='../simple',
                                   iter=1000, init=inits,
                                   chains=1, warmup=2000),
                       regexp='warmup <= iter')
  test <- expect_error(sample_nuts('simple', path='../simple',
                                   iter=1000, init=inits,
                                   chains=-1),
                       regexp='chains >= 1')
  test <- expect_error(sample_nuts('simple3', path='../simple',
                                   iter=1000, init=inits,
                                   chains=1),
                       regexp="Check 'path' and 'model'")
  test <- expect_error(sample_nuts('simple', path='../simple',
                                   iter=1000, init=inits,
                                   chains=1,
                                   control=list(adapt_init_buffer=-20)),
                       regexp='NUTS failed to run')
  test <- expect_error(sample_nuts('simple', path='../simple',
                                   iter=1000, init=inits,
                                   chains=3,
                                   control=list(adapt_delta=1.5)),
                       regexp='NUTS failed to run')
  test <- expect_error(sample_nuts('simple', path='../simple',
                                   iter=1000, init=inits,
                                   chains=3, cores=1,
                                   control=list(adapt_delta=1.5)),
                       regexp='NUTS failed to run')
  test <- expect_error(sample_rwm('simple', path='../simple',
                                   iter=1000, init=inits,
                                   chains=1,
                                   admb_args='-refresh -2'),
                       regexp='RWM failed to run')
})
