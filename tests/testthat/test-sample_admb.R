
## Note these were tested with ADMB 12.2 on 8/3/2020
test_that("simple example works", {
  skip_on_cran()
  inits.fn <- function() list(c(0,0))
  fit <- sample_rwm('simple', path='../simple', chains=1,
                     seeds=1, init=inits.fn,
                     skip_optimization=FALSE,
                     control=list(refresh=-1), skip_monitor=TRUE)
  expect_known_output(extract_samples(fit)[1000,],
                      file='_expect_simple_rwm')
  fit <- sample_rwm('simple', path='../simple', chains=1,
                    seeds=1, init=inits.fn,
                    skip_optimization=FALSE,
                    control=list(refresh=-1), skip_monitor=TRUE)
  expect_known_output(extract_samples(fit)[1000,],
                      file='_expect_simple_rwm')
  fit <- sample_nuts('simple', path='../simple', chains=1,
                     seeds=1, init=inits.fn,
                     skip_optimization=FALSE,
                     control=list(refresh=-1), skip_monitor=TRUE)
  expect_known_output(extract_samples(fit)[1000,],
                      file='_expect_nuts')
  fit <- sample_nuts('simple', path='../simple', chains=1,
                     seeds=1, init=inits.fn,
                     control=list(metric='mle', refresh=-1),
                     skip_monitor = TRUE)
  expect_known_output(extract_samples(fit)[1000,],
                      file='_expect_nuts_mle')
  fit <- sample_nuts('simple', path='../simple', chains=1,
                     seeds=1, init=inits.fn,
                     control=list(metric='mle', refresh=-1),
                     skip_monitor = TRUE)
  expect_known_output(extract_samples(fit)[1000,],
                      file='_expect_nuts_mle')
})

test_that("mceval works",{
  skip_on_cran()
  inits.fn <- function() list(c(0,0))
  fit <- sample_nuts('simple', path='../simple', chains=1,
                     seeds=1, init=inits.fn,
                     control=list(metric='mle', refresh=-1),
                     skip_monitor = TRUE,
                     mceval=TRUE)
})

test_that("parallel works",{
  skip_on_cran()
  message("Starting parallel tests")
  inits.fn <- function() list(c(0,0))
  fit <- sample_nuts('simple', path='../simple', chains=3,
                     seeds=1:3, init=inits.fn, iter=1000,
                     control=list(refresh=-1),
                     skip_monitor = TRUE)
  ## expect_equal(extract_samples(fit)[1500,2], 3.483071)
  fit <- sample_rwm('simple', path='../simple', chains=3,
                     seeds=1:3, init=inits.fn, iter=1000,
                     control=list(refresh=-1),
                     skip_monitor = TRUE)
})


test_that("warnings and errors in sample_nuts and sample_rwm",{
  skip_on_cran()
  inits <- function() list(1,1)
  test <- expect_warning(sample_nuts('simple', path='../simple',
                                     iter=1000, init=inits,
                                     extra.args='-test',
                                     control=list(refresh=-1),
                                     chains=1, warmup=500),
                         regexp='extra.args is deprecated')
  test <- expect_warning(sample_rwm('simple', path='../simple',
                                    iter=1000, init=inits,
                                    extra.args='-test',
                                    control=list(refresh=-1),
                                    chains=1, warmup=500),
                         regexp='extra.args is deprecated')
  test <- expect_warning(sample_nuts('simple', path='../simple',
                                     iter=1000, init=inits,
                                     parallel=TRUE,
                                     control=list(refresh=-1),
                                     chains=1, warmup=500),
                         regexp='parallel is deprecated')
  test <- expect_warning(sample_rwm('simple', path='../simple',
                                    iter=1000, init=inits,
                                    parallel=TRUE,
                                    control=list(refresh=-1),
                                    chains=1, warmup=500),
                         regexp='parallel is deprecated')
  test <- expect_warning(sample_rwm('simple', path='../simple',
                                    iter=1000, init=inits,
                                    control=list(refresh=-1, metric='mle'),
                                    chains=1, warmup=500),
                         regexp='Only refresh control argument is used')
  test <- expect_warning(sample_nuts('simple', path='../simple',
                                     iter=1000, init=inits,
                                     chains=1, warmup=500,
                                     control=list(refresh=-1, metric=diag(2))),
                         regexp='admodel.cov overwritten')
  test <- expect_warning(sample_rwm('simple', path='../simple',
                                    iter=1000, init=inits,
                                    parallel=TRUE,
                                    control=list(refresh=-1),
                                    chains=1, warmup=500),
                         regexp='parallel is deprecated')
  test <- expect_warning(sample_admb('simple', path='../simple',
                                     iter=1000, init=inits,
                                     control=list(refresh=-1),
                                     chains=1, warmup=500),
                       regexp='sample_admb is deprecated')
  test <- expect_error(sample_nuts('simple', path='../simple',
                                   iter=1000, init=inits,
                                   control=list(refresh=-1),
                                   chains=1, warmup=2000),
                       regexp='warmup <= iter')
  test <- expect_error(sample_nuts('simple', path='../simple',
                                   iter=1000, init=inits, algorithm='NUTS',
                                   chains=1, warmup=2000),
                       regexp='unused argument \\(algorithm')
  test <- expect_error(sample_rwm('simple', path='../simple',
                                   iter=1000, init=inits, algorithm='RWM',
                                   chains=1, warmup=2000),
                       regexp='unused argument \\(algorithm')
  test <- expect_error(sample_nuts('simple', path='../simple',
                                   iter=1000, init=inits,
                                   chains=-1),
                       regexp='chains >= 1')
  test <- expect_error(sample_nuts('simple', path='../simple55',
                                   iter=1000, init=inits,
                                   chains=1),
                       regexp="does not exist. Check argument \'path\'")
  test <- expect_error(sample_nuts('simple3', path='../simple',
                                   iter=1000, init=inits,
                                   chains=1),
                       regexp="not found in specified folder")
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
  test <- expect_error(sample_rwm('simple', path='../simple',
                                  iter=1000, init=inits,
                                  chains=1,
                                  control=list(refresh='a')),
                                  regexp='Invalid refresh value')
  test <- expect_error(sample_rwm('simple', path='../simple',
                                   iter=1000, init=inits,
                                   chains=3, seeds=1),
                       regexp='Length of seeds must match chains')
})

test_that("long file names work ok on Windows",{
  skip_on_cran()
  inits.fn <- function() list(1,1)
  p <- '../simple_long_filename'
  m <- 'simple_long_filename'
  if(.Platform$OS.type=='windows'){
    ## Should give warning
    test <- expect_warning(sample_nuts(m, path=p, chains=3, cores=1,
                                       seeds=1:3, init=inits.fn, iter=1000,
                                       control=list(refresh=-1),
                                       skip_monitor = TRUE),
                           regexp='It appears a shortened')
    test <- expect_warning(sample_nuts(m, path=p, chains=3, cores=3,
                                       seeds=1:3, init=inits.fn, iter=1000,
                                       control=list(refresh=-1),
                                       skip_monitor = TRUE),
                           regexp='It appears a shortened')
  } else {
    test <- sample_nuts(m, path=p, chains=3, cores=1,
                        seeds=1:3, init=inits.fn, iter=1000,
                        control=list(refresh=-1),
                        skip_monitor = TRUE)
  }
})
