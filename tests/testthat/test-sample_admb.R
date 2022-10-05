
## Note these were tested with ADMB 13.0 on 2022-10-05
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
  ff <- '../simple/mceval.dat'
  if(file.exists(ff)) file.remove(ff)
  fit <- sample_nuts('simple', path='../simple', chains=1,
                     seeds=1, init=inits.fn,
                     control=list(metric='mle', refresh=-1),
                     skip_monitor = TRUE,
                     mceval=TRUE)
  expect_true(file.exists(ff))
})

test_that("parallel works",{
  skip_on_cran()
  message("Starting parallel tests")
  inits.fn <- function() list(c(0,0))
  fit1 <- sample_nuts('simple', path='../simple', chains=3,
                     seeds=1:3, init=inits.fn, iter=1000,
                     control=list(refresh=-1),
                     skip_monitor = TRUE)
  ## expect_equal(extract_samples(fit)[1500,2], 3.483071)
  fit2 <- sample_rwm('simple', path='../simple', chains=3,
                     seeds=1:3, init=inits.fn, iter=1000,
                     control=list(refresh=-1),
                     skip_monitor = TRUE)
  expect_true(exists('fit1'))
  expect_true(exists('fit2'))
})

test_that("duration works",{
  skip_on_cran()
  message("Starting duration tests")
  inits.fn <- function() list(c(0,0))
  duration <- 5/60 # 5 seconds
  expect_warning(fit11 <- sample_nuts('simple', path='../simple',
                                      chains=2, cores=1, seeds=1:2,
                                      init=inits.fn, iter=1e9, warmup=50,
                                      control=list(refresh=-1), duration=duration,
                                      skip_monitor = TRUE, skip_unbounded=FALSE),
                 regexp='Incomplete chain lengths')
  expect_true(exists('fit11'))
  expect_true(dim(fit11$samples)[1] == nrow(fit11$sampler_params[[1]]))
  expect_true(dim(fit11$samples)[1] == dim(fit11$samples_unbounded)[1])

  expect_warning(fit12 <- sample_rwm('simple', path='../simple',
                                     chains=2, cores=1, seeds=1:2,
                                     init=inits.fn, iter=1e9, warmup=500,
                                     duration=duration,
                                     skip_monitor = TRUE, skip_unbounded=FALSE),
                 regexp='Incomplete chain lengths')
  expect_true(exists('fit12'))
  expect_true(dim(fit11$samples)[1] == dim(fit11$samples_unbounded)[1])

 # try again with a weird thin rate
  thin <- 13
  expect_warning(fit13 <- sample_nuts('simple', path='../simple',
                                      chains=2, cores=1, seeds=1:2, thin=thin,
                                      init=inits.fn, iter=1e9, warmup=50,
                                      control=list(refresh=-1), duration=duration,
                                      skip_monitor = TRUE, skip_unbounded=FALSE),
                 regexp='Incomplete chain lengths')
  expect_true(exists('fit13'))
  expect_true(dim(fit13$samples)[1] == nrow(fit13$sampler_params[[1]]))
  expect_true(dim(fit13$samples)[1] == dim(fit13$samples_unbounded)[1])

  expect_warning(fit14 <- sample_rwm('simple', path='../simple',
                                     chains=2, cores=1, seeds=1:2,
                                     init=inits.fn, iter=1e9, warmup=500,
                                     duration=duration,
                                     skip_monitor = TRUE, skip_unbounded=FALSE),
                 regexp='Incomplete chain lengths')
  expect_true(exists('fit14'))
  expect_true(dim(fit14$samples)[1] == dim(fit14$samples_unbounded)[1])
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

test_that("verbose option works", {
  skip_on_cran()
  inits.fn <- function() list(c(0,0))
  message("Should be no console output between here....")
  message("Starting verbose NUTS in parallel..")
  fit <- sample_nuts('simple', path='../simple', chains=3,
                     seeds=1:3, init=inits.fn, iter=800,
                     skip_monitor = TRUE, verbose=FALSE)
  message("Starting verbose NUTS in serial..")
  fit <- sample_nuts('simple', path='../simple', chains=1,
                     seeds=1, init=inits.fn, iter=800,
                     skip_monitor = TRUE, verbose=FALSE)
  message("Starting verbose RWM in parallel..")
  fit <- sample_rwm('simple', path='../simple', chains=3,
                     seeds=1:3, init=inits.fn, iter=800,
                     skip_monitor = TRUE, verbose=FALSE)
  message("Starting verbose RWM in serial..")
  fit <- sample_rwm('simple', path='../simple', chains=1,
                     seeds=1, init=inits.fn, iter=800,
                    skip_monitor = TRUE, verbose=FALSE)
  message("... and here")
  expect_true(exists('fit'))
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
