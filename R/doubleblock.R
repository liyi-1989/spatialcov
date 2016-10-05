#' Generating special class of (covariance) matrices
#'
#' This is the function that used to generate the block bandable matrices \code{S}
#' with \code{p1} by \code{p2} dimensions.
#'
#'  @param p1 The number of x axis (latitude)
#'  @param p2 The number of y axis (longitude)
#'  @param M Multiplicative constant for off-diagonal elements
#'  @param a decay rate for the x axis
#'  @param b decay rate for the y axis
#'  @return The block bandable covariance matrix
#'  @examples
#'  doubleblock(5,5,1,1.5,1.5)
#'  @export



doubleblock=function(p1,p2,M,a,b){
  n=p1*p2
  C=matrix(0,n,n)

  for(i1 in 1:p1){
    for(j1 in 1:p1){
      for(i2 in 1:p2){
        for(j2 in 1:p2){
          i=(i1-1)*p1+j1
          j=(i2-1)*p1+j2
          if((i1==i2)&(j1==j2)){
            C[i,j]=2.5
          }else{
            C[i,j]=M*min(abs(i1-i2)^(-a),abs(j1-j2)^(-b))
            # C[i,j]=exp(-0.4*base::norm(as.matrix(c(i1-i2,j1-j2)),type="2"))
          }
        }
      }
    }
  }

  return(C)
}
