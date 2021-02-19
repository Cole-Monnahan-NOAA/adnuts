test_that("consistency of algorithms within platform", {
  skip_on_cran()

  ## These will work across platforms b/c they compare within
  ## each. Running the same chains repeatedly to make sure the
  ## same answer each time.

  ## Check consistency given same init and seeds
  inits.fn <- function() list(c(0,0))
  chains <- 10
  cores <- NULL # use parallel which probably catches errors better
  iter <- 2000  # slightly longer ones to detect subtle divergences
  myequal <- function(fit) length(unique(fit$samples[iter,,3]))==1

  test <- sample_rwm('simple', path='../simple', iter=10, init=inits.fn,
                     skip_optimization=FALSE, chains=1)
  fit <- sample_rwm('simple', path='../simple', chains=chains,
                    iter=iter*10, cores=cores,
                    seeds=rep(45,chains), init=inits.fn,
                    skip_monitor=TRUE,
                    control=list(refresh=-1))
  expect_identical(myequal(fit), TRUE)
  ## These correspond to the 6 options in the metric table in the
  ## vignette.
  seeds <- rep(123,chains)
  ## Initialize with diagonal for first three
  ignore <- file.remove('../simple/admodel.cov') # dont need this
  fit1 <- sample_nuts('simple', path='../simple', chains=chains, iter=iter,
                      seeds=seeds, init=inits.fn,
                      skip_monitor=TRUE,
                      control=list(refresh=-1, adapt_mass=FALSE),
                      cores=cores)
  expect_identical(myequal(fit1), TRUE)
  fit2 <- sample_nuts('simple', path='../simple', chains=chains, iter=iter,
                      seeds=seeds, init=inits.fn,
                      control=list(refresh=-1),
                      skip_monitor=TRUE,
                      cores=cores)
  expect_identical(myequal(fit2), TRUE)
  fit3 <- sample_nuts('simple', path='../simple', chains=chains, iter=iter,
                      seeds=seeds, init=inits.fn,
                      skip_monitor=TRUE,
                      control=list(refresh=-1, adapt_mass_dense=TRUE),
                      cores=cores)
  expect_identical(myequal(fit3), TRUE)
  ## Next three initialize from MLE, need to rerun model to get these
  test <- sample_nuts('simple', path='../simple',  iter=100,
                      init=inits.fn, chains=1, skip_optimization=FALSE)
  fit4 <- sample_nuts('simple', path='../simple', chains=chains, iter=iter,
                      seeds=seeds, init=inits.fn,
                      skip_monitor=TRUE,
                      control=list(refresh=-1, metric='mle'),
                      cores=cores)
  expect_identical(myequal(fit4), TRUE)
  fit5 <- sample_nuts('simple', path='../simple', chains=chains, iter=iter,
                      seeds=seeds, init=inits.fn,
                      skip_monitor=TRUE,
                      control=list(refresh=-1, metric='mle', adapt_mass=TRUE),
                      cores=cores)
  expect_identical(myequal(fit5), TRUE)
  fit6 <- sample_nuts('simple', path='../simple', chains=chains, iter=iter,
                      seeds=seeds, init=inits.fn,
                      skip_monitor=TRUE,
                      control=list(refresh=-1, metric='mle', adapt_mass_dense=TRUE),
                      cores=cores)
  expect_identical(myequal(fit6), TRUE)
  ## In addition test passing a user matrix, here unit diag
  fit7 <- suppressWarnings(sample_nuts('simple', path='../simple', chains=chains, iter=iter,
                      seeds=seeds, init=inits.fn,
                      skip_monitor=TRUE,
                      control=list(refresh=-1, metric=diag(2)),
                      cores=cores))
  expect_identical(myequal(fit7), TRUE)
  fit8 <- suppressWarnings(sample_nuts('simple', path='../simple', chains=chains, iter=iter,
                      seeds=seeds, init=inits.fn,
                      skip_monitor=TRUE,
                      control=list(refresh=-1, metric=diag(2), adapt_mass=TRUE),
                      cores=cores))
  expect_identical(myequal(fit8), TRUE)
  fit9 <- suppressWarnings(sample_nuts('simple', path='../simple', chains=chains, iter=iter,
                      seeds=seeds, init=inits.fn,
                      skip_monitor=TRUE,
                      control=list(refresh=-1, metric=diag(2), adapt_mass_dense=TRUE),
                      cores=cores))
  expect_identical(myequal(fit9), TRUE)
  fit10 <- sample_nuts('simple', path='../simple', chains=chains, iter=iter,
                      seeds=seeds, init=inits.fn,
                      skip_optimization = FALSE,
                      skip_monitor=TRUE,
                      control=list(refresh=-1, metric='mle', stepsize=.1),
                      cores=cores)
  expect_identical(myequal(fit10), TRUE)
})
