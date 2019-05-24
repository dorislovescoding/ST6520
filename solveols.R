#' Solve linear system
#'
#' Solve linear system \code{Ax = b} using iteration methods Jacobi.
#'
#' @param A Matrix.
#' @param b Vector. Its dimension should match that of matrix \code{A}.
#' @param c Number of cores used in parallel computation.
#' @param type Type of iteration method, 'GS' means Gauss-Seidel,
#' @param iter Number of iterations in the calculation.
#' 'Jacobi' means Jacobi. Default is 'GS'.
#' @param iter Number of iterations in the calculation. Default is 1000.
#'
#' @return The solution vector \code{x} for linear system \code{Ax = b}.
#' @export
#'

solveols <- function(A,b,iter,c=1,type="GS"){

  gs <-function(GS,A,b,iter){
    n_row=nrow(A)
    x_GS=matrix(0,nrow=n_row, ncol=1)
    for (i in 1:iter){
      for (k in 1:n_row){
        x_GS[k]=(b[k]-A[k,]%*%x_GS+A[k,k]*x_GS[k])/A[k,k]
      }
      x_GS[!is.finite(x)] <- -1e+308

    }
    x_GS
  }

  Jacobi <-function(Jacobi,A,b,iter){
  n_row=nrow(A)
  x_J=matrix(0,nrow=n_row, ncol=1)
  for (i in 1:iter){
    for (j in 1:n_row){
      row=matrix(A[j,],nrow=1,ncol=n_row)
      Arow=row[,-j]
      xrow=x_J[-j,]
      x_J[j]=(b[j]-Arow%*%xrow)/A[j,j]
    }
    x_J[!is.finite(x)] <- -1e+308
  }
  x_J
}

#define parallel jacobi


parallel <-function(A,b,c,iter){
  require(doParallel)
  cl <- makeCluster(c)
  registerDoParallel(cl)
  n_row=nrow(A)
  x_JP=matrix(0,nrow=n_row, ncol=1)
  for (i in 1:iter){
    outlist = foreach (j=1:n_row) %dopar%{
      row=matrix(A[j,],nrow=1,ncol=n_row)
      Arow=row[,-j]
      xrow=x_JP[-j]
      x_JP[j]=(b[j]-Arow%*%xrow)/A[j,j]
    }
    x_JP=unlist(outlist)
    x_JP[!is.finite(x)] <- -1e+308

  }
  x_JP
}

if (type=="Jacobi"){
  if (c > 1) {
    x= parallel(A,b,c,iter)
  }
  else{
    x= Jacobi(Jacobi,A,b,iter)
  }
   }else{
   x= gs(GS,A,b,iter)
 }
  x
}
