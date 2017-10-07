#' adnuts: No-U-turn sampling for Template Model Builder and AD Model Builder
#'
#' Draw Bayesian posterior samples from a TMB or ADMB model using the
#' no-U-turn MCMC sampler. Adaptation schemes are used so specifying tuning
#' parameters is not necessary, and parallel execution reduces overall run
#' time.
#'
#' @details
#' The software package Stan pioneered the use of no-U-turn (NUTS) sampling
#' for Bayesian models (Hoffman and Gelman 2014, Carpenter et
#' al. 2017). This algorithm provides fast, efficient sampling across a
#' wide range of models, including hierarchial ones, and thus can be used
#' as a generic modeling tool (Monnahan et al. 2017). The functionality
#' provided by \pkg{adnuts} is based loosely off Stan and \R package
#' \pkg{rstan}
#'
#' \pkg{adnuts} \R package provides NUTS sampling for two existing software
#' platforms: ADMB (Fournier et al. 2011) and TMB (Kristensen et al. 2017,
#' Kristensen 2017). The specific NUTS capabilities include adaptation of
#' step size and metric (mass matrix), parallel execution, and links to
#' diagnostic and inference tools provided by \pkg{rstan} and
#' \pkg{shinystan}.
#'
#' For TMB models, \pkg{adnuts} provides NUTS and other MCMC algorithms
#' written in \R. These can be used with a TMB model by plugging in the
#' \code{obj$fn} and \code{obj$gr} functions from the DLL directly. It is
#' possible to use these functions with models outside TMB, as long as the
#' log density and gradients can be calculated. See
#' \code{\link{sample_tmb}} for more details.
#'
#' The ADMB implementation is different in that the NUTS code is bundled
#' into the ADMB source itself. Thus, when a user builds an ADMB model the
#' NUTS code is incorporated into the model executable. Thus, \pkg{adnuts}
#' simply provides a convenient set of wrappers to more easiy execute,
#' diagnose, and make inference on a model.
#'
#' @references
#' Carpenter, B., Gelman, A., Hoffman, M.D., Lee, D., Goodrich, B.,
#'   Betancourt, M., Riddell, A., Guo, J.Q., Li, P., Riddell, A.,
#'   2017. Stan: A Probabilistic Programming Language.  J Stat
#'   Softw. 76:1-29.
#'
#' Fournier, D.A., Skaug, H.J., Ancheta, J., Ianelli, J., Magnusson, A.,
#'   Maunder, M.N., Nielsen, A., Sibert, J., 2012. AD Model Builder: using
#'   automatic differentiation for statistical inference of highly
#'   parameterized complex nonlinear models.  Optim Method
#'   Softw. 27:233-249.
#'
#' Hoffman, M.D., Gelman, A., 2014. The no-U-turn sampler: adaptively
#'   setting path lengths in Hamiltonian Monte Carlo.  J Mach Learn
#'   Res. 15:1593-1623.
#'
#' Kristensen, K., Nielsen, A., Berg, C.W., Skaug, H., Bell, B.M.,
#'   2016. TMB: Automatic differentiation and Laplace approximation.  J
#'   Stat Softw. 70:21.
#'
#' Kristensen, K., 2017. TMB: General random effect model builder tool
#'   inspired by ADMB. R package version 1.7.11.
#'
#' Monnahan, C.C., Thorson, J.T., Branch, T.A., 2017. Faster estimation of
#'   Bayesian models in ecology using Hamiltonian Monte Carlo.  Methods in
#'   Ecology and Evolution. 8:339-348.
#'
#' Stan Development Team, 2016. Stan modeling language users guide and
#'   reference manual, version 2.11.0.
#'
#' Stan Development Team, 2016. RStan: The R interface to Stan. R package
#' version 2.14.1. http://mc-stan.org.
#'
#' @docType package
#' @name adnuts
#' @importFrom stats rnorm runif cov
#' @importFrom utils read.csv read.table write.table
NULL
