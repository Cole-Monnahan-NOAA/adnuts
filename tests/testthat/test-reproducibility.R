test_that("reproducibility of algorithms", {
  skip_on_cran()
  ## Check reproducibility given same init and seeds
  inits.fn <- function() list(c(0,0))
  chains <- 1
  fit <- sample_rwm('simple', path='../simple', chains=chains,
                    iter=400, cores=1,
                    seeds=rep(45,chains), init=inits.fn,
                    skip_optimization=FALSE,
                    control=list(refresh=-1))
  expect_identical(unique(fit$samples[400,,3]), -16.0439)
  ## These correspond to the 6 options in the metric table in the
  ## vignette.
  seeds <- rep(123,chains)
  ## Initialize with diagonal for first three
  ignore <- file.remove('../simple/admodel.cov') # dont need this
  fit1 <- sample_nuts('simple', path='../simple', chains=chains, iter=400,
                      seeds=seeds, init=inits.fn,
                      control=list(refresh=-1, adapt_mass=FALSE),
                      cores=1)
  expect_identical(unique(fit1$samples[400,,3]), -12.9319)
  fit2 <- sample_nuts('simple', path='../simple', chains=chains, iter=400,
                      seeds=seeds, init=inits.fn,
                      control=list(refresh=-1),
                      cores=1)
  expect_identical(unique(fit2$samples[400,,3]), -13.2107)
  fit3 <- sample_nuts('simple', path='../simple', chains=chains, iter=400,
                      seeds=seeds, init=inits.fn,
                      control=list(refresh=-1, adapt_mass_dense=TRUE),
                      cores=1)
 # expect_identical(unique(fit3$samples[400,,3]), -14.2902)
  ## Next three initialize from MLE, need to rerun model to get these
  fit4 <- sample_nuts('simple', path='../simple', chains=chains, iter=400,
                      seeds=seeds, init=inits.fn,
                      skip_optimization=FALSE,
                      control=list(refresh=-1, metric='mle'),
                      cores=1)
  expect_identical(unique(fit4$samples[400,,3]), -12.1684)
  fit5 <- sample_nuts('simple', path='../simple', chains=chains, iter=400,
                      seeds=seeds, init=inits.fn,
                      control=list(refresh=-1, metric='mle', adapt_mass=TRUE),
                      cores=1)
  expect_identical(unique(fit5$samples[400,,3]), -12.2534)
  fit6 <- sample_nuts('simple', path='../simple', chains=chains, iter=400,
                      seeds=seeds, init=inits.fn,
                      control=list(refresh=-1, metric='mle', adapt_mass_dense=TRUE),
                      cores=1)
  expect_identical(unique(fit6$samples[400,,3]), -12.4441)
  ## In addition test passing a user matrix, here unit diag
  fit7 <- suppressWarnings(sample_nuts('simple', path='../simple', chains=chains, iter=400,
                      seeds=seeds, init=inits.fn,
                      control=list(refresh=-1, metric=diag(2)),
                      cores=1))
  expect_identical(unique(fit7$samples[400,,3]), -12.9319)
  fit8 <- suppressWarnings(sample_nuts('simple', path='../simple', chains=chains, iter=400,
                      seeds=seeds, init=inits.fn,
                      control=list(refresh=-1, metric=diag(2), adapt_mass=TRUE),
                      cores=1))
  expect_identical(unique(fit8$samples[400,,3]), -13.2107)
  fit9 <- suppressWarnings(sample_nuts('simple', path='../simple', chains=chains, iter=400,
                      seeds=seeds, init=inits.fn,
                      control=list(refresh=-1, metric=diag(2), adapt_mass_dense=TRUE),
                      cores=1))
  expect_identical(unique(fit9$samples[400,,3]), -14.2902)
  ## All of these test might fail if changes to the adaptation
  ## schemes (stepsize or mass matrix) are done in the ADMB
  ## source. So one last tests which uses no adaptation so should
  ## be consistent between ADMB versions. Also need to reoptimize
  ## since I overwrite the admodel.cov file above
  fit10 <- sample_nuts('simple', path='../simple', chains=chains, iter=400,
                      seeds=seeds, init=inits.fn,
                      skip_optimization = FALSE,
                      control=list(refresh=-1, metric='mle', stepsize=.1),
                      cores=1)
  expect_identical(unique(fit10$samples[400,,3]), -13.6047)

})
