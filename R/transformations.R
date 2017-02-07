#' Determine which transformation is to be used for each parameter.
#'
#' @details The 4 cases are (0) none, (1) lower only [a,Inf], (2) upper
#'   only [-Inf, b], and (3) both lower and upper [a,b]. Each case requires
#'   a different set of transformation functions.
#' @param lower Vector of lower bounds, potentially -infinity for some
#' @param upper Vector of upper bounds, potentially infinity for some
#' @return Vector of cases, in 0:3, to be used in all other transformation
#'   functions. Error checking is only done here, not in other functions.
#' @seealso \code{\link{.transform}}, \code{\link{.transform.inv}},
#'   \code{\link{.transform.grad}}, \code{\link{.transform.grad2}}
#' @export
#'
.transform.cases <- function(lower, upper){
  if(length(lower) != length(upper))
    stop("Lengths of lower and upper do not match")
  if(any(is.na(c(lower, upper))) | any(is.nan(c(lower, upper))))
    stop("Bounds must be finite or -Inf/Inf -- NA and NaN not allowed")
  if(any(lower >= upper))
    stop("Lower bound >= upper bound")
  cases <- rep(NA, length(lower))
  cases[!is.finite(lower) & !is.finite(upper)] <- 0
  cases[is.finite(lower) & !is.finite(upper)] <- 1
  cases[!is.finite(lower) & is.finite(upper)] <- 2
  cases[is.finite(lower) & is.finite(upper)] <- 3
  if(any(is.na(cases)))
    stop("Something unexpected went wrong determining the bounding functions.
 Check lower and upper.")
  return(cases)
}

#' This function returns the transformed variable, x=f(y).
#'
#' @export
.transform <- function(y, a, b, cases){
  x <- y
  ind <- cases==1
  if(length(ind)>0)
    x[ind] <- exp(y[ind])+a[ind]
  ind <- cases==2
  if(length(ind)>0)
    x[ind] <- b[ind]-exp(y[ind])
  ind <- cases==3
  if(length(ind)>0)
   x[ind] <- a[ind]+(b[ind]-a[ind])/(1+exp(-y[ind]))
  return(x)
}

#' The inverse of the transformation, y=f-1(x).
#' @export
.transform.inv <- function(x, a, b, cases){
  if(any(x<a) | any(x>b)) stop("x outside limits provided -- not meaningful")
  y <- sapply(1:length(x), function(i) {
    if(cases[i]==0) return(x[i])
    else if(cases[i]==1) return(log(x[i]-a[i]))
    else if(cases[i]==2) return(log(b[i]-x[i]))
    else if(cases[i]==3) return(-log( (b[i]-x[i])/(x[i]-a[i]) ))
  })
  return(y)
}

#' The absolute value of the derivative of transformation.
#' @export
.transform.grad <- function(y, a, b, cases){
  x <- rep(1, length(y))
  ind <- cases %in% 1:2
  if(length(ind)>0)
    x[ind] <- exp(y[ind])
  ind <- cases==3
  if(length(ind)>0){
    tmp <- exp(-y[ind])
    x[ind] <- (b[ind]-a[ind])*tmp/(1+tmp)^2
  }
  return(x)
}

#' The derivative of the log of the derivative of the transformation. I.e.,
#' d/dy[log(.transform.grad(y,a,b))].
#' @export
.transform.grad2 <- function(y, a, b, cases){
  x <- rep(0, len=length(y))
  ind <- cases %in% 1:2
  if(length(ind)>0)
    x[ind] <- 1
  ind <- cases==3
  if(length(ind)>0)
    x[ind] <- -1+2/(1+exp(y[ind]))
  return(x)
}
