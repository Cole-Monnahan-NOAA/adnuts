test_that("all three metrics work", {
 skip_if(skip_TMB)
  TMB::runExample('simple')
  Q <- sdreport(obj, getJointPrecision = TRUE)$jointPrecision
  M <- as.matrix(solve(Q))
  fits <- list()
  library(StanEstimators)
  for(m in c('dense', 'sparse', 'diag')){
    fits[[m]] <- sample_sparse_tmb(obj, iter=1000,
                                   warmup=200, cores=1, chains=1, seed=1,
                                   metric=m)
  }
 expect_equal(length(fits),3)
})

test_that("Laplace works", {
  skip_if(skip_TMB)
  TMB::runExample('simple')
  fit1 <- sample_sparse_tmb(obj, iter=1000,
                            warmup=200, cores=1, chains=1, seed=1,
                            metric='dense')
   # has RE and uses Laplace
  fit2 <- sample_sparse_tmb(obj, iter=1000, laplace=TRUE,
                            warmup=200, cores=1, chains=1, seed=1,
                            metric='dense')
  expect_equal(ncol(as.data.frame(fit1)),118)
  expect_equal(ncol(as.data.frame(fit2)),4)
})
