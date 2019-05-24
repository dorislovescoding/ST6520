#' Subsampling regression methods
#'
#' Implements algorithmic leveraging for linear regression using uniform or
#'      leverage score based subsampling of rows, i.e, fit \code{y = X beta + E}.
#'
#' @param X The covariates matrix of size \code{n} by \code{p}.
#' @param y The response vector with length matching to the row number of \code{X}.
#' @param samplesize Number of subsamples to choose for regression.
#' @param seed set the seed
#' @param Type The type of subsampling rows. "Lev" & "Uniform", the default is "Uniform"
#' @return The coefficient vector \code{beta}.
#' @export

algoleverage <- function(X,y,samplesize,Type="Uniform",seed) {
  set.seed(seed)
  n = length(y)

  H=X%*%solve(t(X)%*%X)%*%t(X)
  P_lev=diag(H)
  P_uni=rep(1,n)
  W = diag(rep(1, samplesize))

  if (Type =="lev"){
    subsample=sample(1:n,samplesize,replace=FALSE,prob=P_lev)
    W = diag(1/P_lev[subsample])

  } else{
    subsample=sample(1:n,samplesize,replace=FALSE,prob=P_uni)
  }
  sample_x=X[subsample,]
  sample_y=y[subsample]

  b=solve(t(sample_x)%*%W%*%sample_x)%*%t(sample_x)%*%W%*%sample_y

  return(b)
}
