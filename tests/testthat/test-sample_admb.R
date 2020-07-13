test_that("simple example works", {
  skip_on_cran()
  inits.fn <- function() list(c(0,0))
  oldwd <- getwd()
  fit <- sample_admb('simple', path='../simple', chains=1,
                     seeds=1, init=inits.fn, algorithm='RWM',
                     skip_optimization=FALSE,
                     control=list(refresh=0), skip_monitor=TRUE)
  expect_known_output(extract_samples(fit)[1000,],
                      file='_expect_simple_rwm')
  fit <- sample_admb('simple', path='../simple', chains=1,
                     seeds=1, init=inits.fn,
                     skip_optimization=FALSE,
                     control=list(refresh=0), skip_monitor=TRUE)
  expect_known_output(extract_samples(fit)[1000,],
                      file='_expect_nuts')
  fit <- sample_admb('simple', path='../simple', chains=1,
                     seeds=1, init=inits.fn,
                     control=list(metric='mle', refresh=0),
                     skip_monitor = TRUE)
  expect_known_output(extract_samples(fit)[1000,],
                      file='_expect_nuts_mle')
})

test_that("reproducibility among chains", {
  skip_on_cran()
  ## Check reproducibility given same init and seeds
  inits.fn <- function() list(c(0,0))
  oldwd <- getwd()
  chains <- 10
  fit2 <- sample_admb('simple', path='../simple', chains=chains,
                     iter=400, algorithm='RWM', cores=1,
                     seeds=rep(45,chains), init=inits.fn,
                     skip_optimization=FALSE,
                     admb_args=' -refresh 0')
  expect_identical(unique(fit2$samples[400,,3]), -16.0439)
  fit1 <- sample_admb('simple', path='../simple', chains=chains, iter=400,
                      seeds=rep(123,chains), init=inits.fn,
                      skip_optimization=FALSE,
                     cores=1, admb_args=' -refresh 0')
  expect_identical(unique(fit1$samples[400,,3]), -13.2107)
  fit3 <- sample_admb('simple', path='../simple', chains=chains, iter=400,
                      seeds=rep(123,chains), init=inits.fn, cores=1,
                      control=list(metric='mle'), admb_args=' -refresh 0')
  expect_identical(unique(fit3$samples[400,,3]), -12.2534)
  fit4 <- suppressWarnings(sample_admb('simple', path='../simple', chains=chains, iter=400,
                      seeds=rep(123,chains), init=inits.fn, cores=1,
                      control=list(metric=diag(2)), admb_args=' -refresh 0'))
  expect_identical(unique(fit4$samples[400,,3]), -13.2107)
  ## Should match fit 4 since same settings
  fit5 <- sample_admb('simple', path='../simple', chains=chains, iter=400,
                      seeds=rep(123,chains), init=inits.fn, cores=1,
                      control=list(metric='unit'), admb_args=' -refresh 0')
  expect_identical(unique(fit5$samples[400,,3]), -13.2107)
})



