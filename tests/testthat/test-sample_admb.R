test_that("simple example works", {
  skip_on_cran()
  inits.fn <- function() list(c(0,0))
  oldwd <- getwd()
  print(oldwd)
  setwd('../simple'); system('simple'); setwd(oldwd)
  fit <- sample_admb('simple', path='../simple', chains=1,
                     seeds=1, init=inits.fn, algorithm='RWM')
  expect_known_output(extract_samples(fit)[1000,],
                      file='_expect_simple_rwm')
  setwd('../simple'); system('simple -hbf'); setwd(oldwd)
  fit <- sample_admb('simple', path='../simple', chains=1,
                     seeds=1, init=inits.fn, control=list(metric=NULL))
  expect_known_output(extract_samples(fit)[1000,],
                      file='_expect_nuts')
  fit <- sample_admb('simple', path='../simple', chains=1,
                     seeds=1, init=inits.fn, control=list(metric='mle'))
  expect_known_output(extract_samples(fit)[1000,],
                      file='_expect_nuts_mle')
})

test_that("reproducibility among chains", {
  skip_on_cran()
  ## Check reproducibility given same init and seeds
  inits.fn <- function() list(c(0,0))
  oldwd <- getwd()
  chains <- 50
  setwd('../simple'); system('simple'); setwd(oldwd)
  fit2 <- sample_admb('simple', path='../simple', chains=chains,
                     iter=400, algorithm='RWM',
                     seeds=rep(45,chains), init=inits.fn)
  expect_identical(unique(fit2$samples[400,,3]), -17.0572)
  setwd('../simple'); system('simple -hbf'); setwd(oldwd)
  fit1 <- sample_admb('simple', path='../simple', chains=chains, iter=400,
                     seeds=rep(123,chains), init=inits.fn)
  expect_identical(unique(fit1$samples[400,,3]), -11.5914)
  fit3 <- sample_admb('simple', path='../simple', chains=chains, iter=400,
                      seeds=rep(123,chains), init=inits.fn,
                      control=list(metric='mle'))
  expect_identical(unique(fit3$samples[400,,3]), -12.0708)
  fit4 <- suppressWarnings(sample_admb('simple', path='../simple', chains=chains, iter=400,
                      seeds=rep(123,chains), init=inits.fn,
                      control=list(metric=diag(2))))
  expect_identical(unique(fit4$samples[400,,3]), -12.7999)
  ## Should match fit 4 since same settings
  fit5 <- sample_admb('simple', path='../simple', chains=chains, iter=400,
                      seeds=rep(123,chains), init=inits.fn,
                      control=list(metric='unit'))
  expect_identical(unique(fit5$samples[400,,3]), -12.7999)
})



