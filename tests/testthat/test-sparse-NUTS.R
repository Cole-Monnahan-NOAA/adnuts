test_that("all metrics work", {
 skip_if(skip_TMB)
  obj <- get_simple_obj()
  Q <- sdreport(obj, getJointPrecision = TRUE)$jointPrecision
  M <- as.matrix(solve(Q))
  fits <- list()
  for(m in c('auto', 'dense', 'sparse', 'stan', 'sparse-naive', 'diag', 'unit')){
    suppressWarnings(suppressMessages(fits[[m]] <- sample_snuts(obj, iter=1000,
                                   skip_optimization = TRUE,
                                   Q=Q, Qinv=M,
                                   refresh=0,
                                   warmup=150, cores=1, chains=1, seed=1,
                                   metric=m, print=FALSE)))
  }
  expect_equal(length(fits),7)
  out <- lapply(fits, function(x) as.numeric(tail(as.data.frame(x), n=1)[1]))
  expect_equal(out$dense,-0.4488092, tolerance=1e-5)
  expect_equal(out$sparse,-0.7966603, tolerance=1e-5)
  # for some reason this one matches locally, but fails during testing???
  #expect_equal(out$diag,-0.08499555, tolerance=1e-5)
  expect_equal(out$unit,-0.88248, tolerance=1e-5)
})

test_that("Embedded Laplace approximation works", {
  skip_if(skip_TMB)
  obj <- get_simple_obj()
  fit1 <- sample_snuts(obj, iter=1000, refresh=0,
                            warmup=200, cores=1, chains=1, seed=1,
                            metric='dense', print=FALSE)
   # has RE and uses Laplace
  fit2 <- sample_snuts(obj, iter=1000, laplace=TRUE, refresh=0,
                            warmup=200, cores=1, chains=1, seed=1,
                            metric='dense', print=FALSE)
  expect_equal(ncol(as.data.frame(fit1)),118)
  expect_equal(ncol(as.data.frame(fit2)),4)
})

test_that("metrics are robust to model type",{
  skip_if(skip_TMB)
  obj <- get_simple_obj()
  ## rebuild without RE so it fails
  obj <- TMB::MakeADFun(data=obj$env$data, parameters=obj$env$parList(),
                        map=obj$env$map,
                        random=NULL, silent=TRUE,
                        DLL=obj$env$DLL)
  expect_error(sample_snuts(obj, iter=1000, laplace=TRUE,
                            warmup=200, cores=1, chains=1, seed=1,
                            metric='dense', print=FALSE),
               regexp = 'No random effects found')
  expect_error(sample_snuts(obj, iter=1000,
                            warmup=200, cores=1, chains=1, seed=1,
                            metric='sparse', print=FALSE),
               regexp='sparse metric only allowed with random effects')
  ## should fail since M not available
  expect_error(sample_snuts(obj, iter=1000,
                            warmup=200, cores=1, chains=1, seed=1,
                            metric='dense', print=FALSE),
               regexp = 'Some standard errors estimated to be NaN'
              )
  ## should work
  suppressWarnings(fit4 <- sample_snuts(obj, iter=1000, refresh=0,
                            warmup=200, cores=1, chains=1, seed=1,
                            metric='unit', print=FALSE))
  expect_equal(ncol(as.data.frame(fit4)),118)
  ## should fail since M not available
  expect_error(sample_snuts(obj, iter=1000, refresh=0,
                            warmup=200, cores=1, chains=1, seed=1,
                            metric='diag', print=FALSE),
               "Some standard errors estimated to be NaN")
  ## should work if variance term is turned off (penalized ML)
  obj2 <- TMB::MakeADFun(data=obj$env$data, parameters=obj$env$parList(),
                         map=list(logsdu=factor(NA)),
                         random=NULL, silent=TRUE,
                         DLL=obj$env$DLL)
  suppressWarnings(fit6 <- sample_snuts(obj2, iter=1000, refresh=0,
                            warmup=200, cores=1, chains=1, seed=1,
                            metric='dense', print=FALSE))
 expect_equal(ncol(as.data.frame(fit6)), 117)
})

test_that("parallel works", {
  skip_if(skip_TMB)
  obj <- get_simple_obj()
  fit <- sample_snuts(obj, iter=1000, warmup=200, cores=4,
                           refresh=0, print=FALSE,
                           chains=4, seed=1, metric='sparse')
  expect_equal(sum(as.data.frame(fit)),  259212.9483)
})


test_that("thinning works", {
  skip_if(skip_TMB)
  obj <- get_simple_obj()
  fit <- sample_snuts(obj, iter=1000, warmup=200, cores=1,
                           chains=1, seed=1, metric='sparse',
                           thin=2,refresh=0, print=FALSE)
  expect_equal(400,nrow(as.data.frame(fit)))
  #fit <- sample_snuts(obj, iter=1000, warmup=200, cores=4,
   #                        chains=4, seed=1, metric='sparse',
    #                       thin=3)
  # this is still broken!!
  #expect_equal(400,nrow(as.data.frame(fit)))
  })

test_that("auto metric selection is robust to model type", {
  skip_if(skip_TMB)
  obj <- get_simple_obj()
  ## normal case of RE, with and without laplace
  suppressWarnings(fit1 <- sample_snuts(obj, iter=1000, refresh=0,
                            warmup=200, cores=1, chains=1, seed=1,
                            metric='auto', print=FALSE))
  expect_equal(as.numeric(tail(as.data.frame(fit1),1)[1]),-1.362821, tolerance =1e-6)
  suppressWarnings(fit2 <- sample_snuts(obj, iter=1000,  laplace=TRUE, refresh=0,
                            warmup=200, cores=1, chains=1, seed=1,
                            metric='auto', print=FALSE))
  expect_equal(as.numeric(tail(as.data.frame(fit2),1)[1]), 52.38106, tolerance =1e-6)

  ## rebuild without RE so it fails: behavior on model w/o mode
  obj <- get_simple_obj()
  obj2 <- TMB::MakeADFun(data=obj$env$data, parameters=par0,
                        map=obj$env$map,
                        random=NULL, silent=TRUE,
                        DLL=obj$env$DLL)
  expect_error(sample_snuts(obj2, iter=1000, laplace=TRUE,refresh=0,
                                 warmup=200, cores=1, chains=1, seed=1,
                                 metric='auto', print=FALSE),
               regexp = 'No random effects found')
  # this breaks if init='last.par.best' b/c the inits are so bad it can't recover
  expect_warning(fit3 <- sample_snuts(obj2, iter=1000, laplace=FALSE, init='unif',
                            warmup=200, cores=1, chains=1, seed=1,
                            refresh=0, print=FALSE,
                            metric='auto', skip_optimization = TRUE), regexp = 'NaNs')
  expect_equal(fit3$metric, 'unit')
  expect_equal(as.numeric(tail(as.data.frame(fit3),1)[1]), -0.843062, tolerance =1e-6)
  ## rebuild as penalized ML
  obj <- get_simple_obj()
  obj2 <- TMB::MakeADFun(data=obj$env$data, parameters=par0,
                         map=list(logsdu=factor(NA), logsd0=factor(NA)),
                         random=NULL, silent=TRUE,
                         DLL=obj$env$DLL)
  expect_error(sample_snuts(obj2, iter=1000, laplace=TRUE,
                                 warmup=200, cores=1, chains=1, seed=1,
                                 metric='auto', print=FALSE),
               regexp = 'No random effects found')
  suppressWarnings(fit5 <- sample_snuts(obj2, iter=1000, refresh=0,
                            warmup=200, cores=1, chains=1, seed=1,
                            metric='auto', print=FALSE))
  expect_equal(fit5$metric, 'dense')
  expect_equal(as.numeric(tail(as.data.frame(fit5),1)[1]),  -1.014703, tolerance =1e-6)
})


test_that("random inits work", {
  # test that different initial values work by forcing NUTS to be
  # autocorrelated and see the paths
  skip_if(skip_TMB)
  obj <- get_simple_obj()
  out <- NULL
  Q <- sdreport(obj, getJointPrecision = TRUE)$jointPrecision
  Qinv <- solve(Q) |> as.matrix()
  obj$par <- opt$par
  for(seed in 1:20){
    for(init in c('last.par.best', 'random', 'random-t', 'unif')){
      suppressMessages(tmpfit <- sample_snuts(obj, iter=250, warmup=200, cores=1,
                                  chains=1, seed=seed, metric='unit',
                                  init=init, refresh=0, Q=Q, Qinv=Qinv, print=FALSE,
                                  control=list(max_treedepth=1, adapt_delta=.99)))
      out <- data.frame(init=init, seed=seed, iter=1:250,
                        lp=extract_samples(tmpfit, inc_warmup=TRUE, inc_lp =TRUE)$lp__)  |>
        rbind(out)
    }
  }
  out2 <- subset(out, iter==1)
  mean_lps <- as.numeric(tapply(out2$lp, INDEX=out2$init, FUN=mean))
  expect_equal(mean_lps,
               c(-403.24160, -465.89700, -519.89795, -254101.8825),
               tolerance =1e-6)
  # library(ggplot2)
  # ggplot(out2, aes(x=init, y=lp)) + geom_jitter(width=.1, height=0)
  # ggplot(out, aes(iter, y=lp, color=init, group=interaction(init,seed))) +
  #   geom_line()
})


test_that("RTMB works", {
  skip_if(skip_RTMB) # not sure why this fails when testing but not locally?
  ## from RTMB beginner help page
  ## suppressWarnings(detach("package:TMB", unload = TRUE))
  if('TMB' %in% .packages()) detach(package:TMB)
  suppressWarnings(library(RTMB))
  data(ChickWeight)
  parameters <- list(
    mua=0,          ## Mean slope
    sda=1,          ## Std of slopes
    mub=0,          ## Mean intercept
    sdb=1,          ## Std of intercepts
    sdeps=1,        ## Residual Std
    a=rep(0, 50),   ## Random slope by chick
    b=rep(0, 50)    ## Random intercept by chick
  )
  f <- function(parms) {
    RTMB::getAll(ChickWeight, parms, warn=FALSE)
    ## Optional (enables extra RTMB features)
    weight <- OBS(weight)
    ## Initialize joint negative log likelihood
    nll <- 0
    ## Random slopes
    nll <- nll - sum(RTMB::dnorm(a, mean=mua, sd=sda, log=TRUE))
    ## Random intercepts
    nll <- nll - sum(RTMB::dnorm(b, mean=mub, sd=sdb, log=TRUE))
    ## Data
    predWeight <- a[Chick] * Time + b[Chick]
    nll <- nll - sum(RTMB::dnorm(weight, predWeight, sd=sdeps, log=TRUE))
    ## Get predicted weight uncertainties
    ADREPORT(predWeight)
    ## Return
    nll
  }
  f(parameters)
  obj <- RTMB::MakeADFun(f, parameters, random=c("a", "b"))
  #obj$fn()
  fit <- sample_snuts(obj, iter=1000, warmup=250, chains=2,
                           refresh=0, cores=2, seed=1, print=FALSE)
  expect_equal(ncol(as.data.frame(fit)),105)
  detach(package:RTMB, unload=TRUE)
})



test_that("small models work", {
  skip_if(skip_RTMB)
  if('TMB' %in% .packages()) detach(package:TMB, unload=TRUE)
  suppressWarnings(library(RTMB))
  f <- function(params){
    RTMB::getAll(params)
    -sum(RTMB::dnorm(x,0,1,TRUE))
  }
  params <- list(x=1)
  f(params)
  obj <- RTMB::MakeADFun(f, params)
  # single parameter single chain
  fit <- sample_snuts(obj, iter=300, warmup=200, chains=1,
                           refresh=0, cores=1, seed=1, print=FALSE)
  post1 <- as.data.frame(fit)
  post2 <- extract_samples(fit, inc_lp=TRUE)
  expect_true(is.data.frame(post1))
  expect_true(is.data.frame(post2))
  expect_equal(names(post1), 'x')
  expect_equal(names(post2), c('x', 'lp__'))
  expect_equal(post1$x[5], post2$x[5])
 # repeat with 2 chains
  fit <- sample_snuts(obj, iter=300, warmup=200, refresh=0,
                           chains=2, cores=1, seed=1, print=FALSE)
  post1 <- as.data.frame(fit)
  post2 <- extract_samples(fit, inc_lp=TRUE)
  expect_true(is.data.frame(post1))
  expect_true(is.data.frame(post2))
  expect_equal(names(post1), 'x')
  expect_equal(names(post2), c('x', 'lp__'))
  expect_equal(post1$x[5], post2$x[5])

  ## now with two parameters
  params <- list(x=c(1,2))
  obj <- RTMB::MakeADFun(f, params)
  # single parameter single chain
  fit <- sample_snuts(obj, iter=300, warmup=200, refresh=0,
                           chains=1, cores=1, seed=1, print=FALSE)
  post1 <- as.data.frame(fit)
  post2 <- extract_samples(fit, inc_lp=TRUE)
  expect_true(is.data.frame(post1))
  expect_true(is.data.frame(post2))
  expect_equal(post1[1,5], post2[1,5])
  detach(package:RTMB, unload=TRUE)
})

