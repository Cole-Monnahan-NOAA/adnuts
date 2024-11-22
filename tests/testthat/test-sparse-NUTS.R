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

test_that("metrics are robust to model type",{
  skip_if(skip_TMB)
  TMB::runExample('simple')
  ## rebuild without RE so it fails
  obj <- TMB::MakeADFun(data=obj$env$data, parameters=obj$env$parList(),
                        map=obj$env$map,
                        random=NULL, silent=TRUE,
                        DLL=obj$env$DLL)
  expect_error(sample_sparse_tmb(obj, iter=1000, laplace=TRUE,
                            warmup=200, cores=1, chains=1, seed=1,
                            metric='dense'),
               regexp = 'No random effects found')


  expect_error(sample_sparse_tmb(obj, iter=1000,
                            warmup=200, cores=1, chains=1, seed=1,
                            metric='sparse'),
               regexp='sparse metric only allowed with random effects')
  ## should fail since M not available
  expect_error(sample_sparse_tmb(obj, iter=1000,
                            warmup=200, cores=1, chains=1, seed=1,
                            metric='dense'),
               regexp = 'Some standard errors estimated to be NaN'
              )
  ## should work
  fit4 <- sample_sparse_tmb(obj, iter=1000,
                            warmup=200, cores=1, chains=1, seed=1,
                            metric='unit')
  expect_equal(ncol(as.data.frame(fit4)),118)
  ## should fail since M not available
  expect_error(sample_sparse_tmb(obj, iter=1000,
                            warmup=200, cores=1, chains=1, seed=1,
                            metric='diag'),
               "Some standard errors estimated to be NaN")
  ## should work if variance term is turned off (penalized ML)
  obj2 <- TMB::MakeADFun(data=obj$env$data, parameters=obj$env$parList(),
                         map=list(logsdu=factor(NA)),
                         random=NULL, silent=TRUE,
                         DLL=obj$env$DLL)
  fit6 <- sample_sparse_tmb(obj2, iter=1000,
                            warmup=200, cores=1, chains=1, seed=1,
                            metric='dense')
 expect_equal(ncol(as.data.frame(fit6)), 117)
})


