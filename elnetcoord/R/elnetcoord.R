#' Elastic net regression
#'
#' Solve the elastic net regression problem.
#'
#' @param x Covariate matrix of \code{n} by \code{p}.
#' @param y The depend variable, would be demeaned in the algorithm
#' @param alpha Parameter in \code{[0,1]} specifying the relative weight between L1 penalty
#' and L2 penalty.
#' @param tol Tolerance error between two adjacent coordinate descent iterations.
#' Computed by \code{norm(new.beta - old.beta)}.
#' @param lambdagrid Specify the grid to find the most desirable lambda

#' @export
#
#' @seealso \code{\link[glmnet]{glmnet}}
#' Can also refer to glmnet for a comparison

elnetcoord <- function(x,y,tol,alpha,lambdagrid) {


  y=y-mean(y)
  p=ncol(x)


  sumxy=abs(t(x)%*%y)
  maxl= round(max(sumxy)/alpha,4)+1
  lambda  <- seq(0,min(1000000,maxl),length=lambdagrid)

  beta_cd <- matrix(0,lambdagrid,p)
  iter_cd <- matrix(0,lambdagrid)

  for(l in 1:lambdagrid){
    iter=1
    # Initial values
    b <- rep(0,p)
    r = y - x%*%b
    error=1
    # Coordiante descent
    while (error>tol) {
      b0=b
      for (j in 1:p) {

        # partial residuals
        r = r + x[,j]*b[j]

        # soft-threshold solution
        xr = sum(x[,j]*r)
        xx = sum(x[,j]^2)
        b[j] = (abs(xr)-lambda[l]*alpha)/(xx+lambda[l]*(1-alpha))
        b[j] = sign(xr)*ifelse(b[j]>0,b[j],0)

        # residuals
        r <- r - x[,j]*b[j]
      }
      error=b-b0
      error=norm(error,'2')
      iter=iter+1
    }
    iter_cd[l]=iter
    beta_cd[l,]<-b
  }
  beta_cd
}
