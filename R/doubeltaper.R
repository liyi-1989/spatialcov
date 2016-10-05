#' Estimator for the spatial double block matrices -- double tapering estimator
#'
#' This is the function that used to implement the double tapering estimator for
#' the block bandable matrices \code{S} with \code{p1} by \code{p2} dimensions.
#'
#'
#'  @param S The block bandable matrices
#'  @param p1 The number of x axis (latitude)
#'  @param p2 The number of y axis (longitude)
#'  @param k bandwidth for the x axis
#'  @param l bandwidth for the y axis
#'  @param ck constant in the bandwidth for the x axis
#'  @param cl constant in the bandwidth for the y axis
#'  @param a decay rate for the x axis
#'  @param b decay rate for the y axis
#'  @param n sample size
#'  @param method tapering method, "tapering" for linear tapering, "banding" for banding estimator
#'  @param Axis tapering direction, 0 for tapering for both direction, 1 for x direction only, 2 for y direction only
#'  @return The regularized block bandable covariance matrix with double tapering estimator
#'  @examples
#'  library(MASS)
#'  p1=25; p2=25; p=p1*p2; M=1; ck=cl=10; n=500;
#'  nalpha=1:5
#'  err_hat=err_band=err_taper=rep(0,length(nalpha))
#'  for(i in nalpha){
#'    a=1+i/10
#'    cat("Decay Rate = ",a,"...\n")
#'    S=doubleblock(p1,p2,M,a,a)
#'    data = mvrnorm(n, mu = rep(0,p), Sigma = S)
#'    S_hat=cov(data)
#'    S_band=doubeltaper(S_hat,p1,p2,ck,cl,a,a,n,method="banding",Axis=0)
#'    S_taper=doubeltaper(S_hat,p1,p2,ck,cl,a,a,n,method="tapering",Axis=0)
#'
#'    err_hat[i]=sqrt(sum((S_hat-S)^2))
#'    err_band[i]=sqrt(sum((S_band-S)^2))
#'    err_taper[i]=sqrt(sum((S_taper-S)^2))
#'  }
#'  plot(1+nalpha/10,err_hat/err_hat,col="black",type="o",ylim=c(0,1.5),
#'       xlab="decay rate",ylab="Relative Error",main="Relative Error for Different Estimators")
#'  lines(1+nalpha/10,err_band/err_hat,col="blue",type="o")
#'  lines(1+nalpha/10,err_taper/err_hat,col="red",type="o")
#'  legend("topright",c("Sample Covariance","Banding","Tapering"),col=c("black","blue","red"),
#'         lty=rep(1,3),pch=rep(1,3),cex=0.5)
#'  @export

doubeltaper=function(S,p1,p2,k,l,ck=1,cl=1,a=1.3,b=1.3,n=1000,method="tapering",Axis=0){
  # Input: a covariance matrix S with dimension p=p1p2,
  # based on spatial data X n*p=n*p1*p2
  # Output: S_hat

  # Step 1: check
  d=dim(S)
  if(d[1]!=d[2]){stop("S must be a square matrix!")}
  p=d[1]
  if(d[1]!=p1*p2){stop("Dimension of S must be p1*p2!")}
  # Step 2: parameters
  if(missing(k)){
    pn=(2*b-1)/(4*a*b-1)
    pp1=(2*b)/(4*a*b-1)
    pp2=-1/(4*a*b-1)
    k=n^(pn)*p1^(pp1)*p2^(pp2)
    k=floor(k)
  }
  if(missing(l)){
    pn=(2*a-1)/(4*a*b-1)
    pp1=-1/(4*a*b-1)
    pp2=(2*a)/(4*a*b-1)
    l=n^(pn)*p1^(pp1)*p2^(pp2)
    l=floor(l)
  }
  # if(missing(ck)){
  #   ck=1
  # }
  # if(missing(cl)){
  #   cl=1
  # }
  k=k*ck
  l=l*cl
  if(missing(a)){
    a=1.1
  }
  if(missing(b)){
    b=1.1
  }

  # Step 3: Coefficients

  w=matrix(0,p2,p2)
  v=matrix(0,p1,p1)

  if(method=="banding"){
    # w inside bandable matrix p2*p2
    for(i in 1:p2){
      for(j in 1:p2){
        w[i,j]=ifelse(abs(i-j)<k/2,1,0)
      }
    }
    # v outside bandable matrix p1*p1
    for(i in 1:p1){
      for(j in 1:p1){
        v[i,j]=ifelse(abs(i-j)<l/2,1,0)
      }
    }
  }else if(method=="tapering"){
    # w inside bandable matrix p2*p2
    for(i in 1:p2){
      for(j in 1:p2){
        w[i,j]=ifelse(abs(i-j)<k,ifelse(abs(i-j)<k/2,1,2-2*abs(i-j)/k),0)
      }
    }
    # v outside bandable matrix p1*p1
    for(i in 1:p1){
      for(j in 1:p1){
        v[i,j]=ifelse(abs(i-j)<l,ifelse(abs(i-j)<l/2,1,2-2*abs(i-j)/l),0)
      }
    }
  }else{
    stop("Method must be banding or tapering!")
  }

  # W=matrix(0,p,p)
  # V=matrix(0,p,p)
  # W=matrix(1,p1,p1)%x%w
  # V=v%x%matrix(1,p2,p2)
  # S_hat=V*W*S
  if(Axis==0){
    return((v%x%matrix(1,p2,p2))*(matrix(1,p1,p1)%x%w)*S)
  }else if(Axis==1){
    return((v%x%matrix(1,p2,p2))*S)
  }else if(Axis==2){
    return((matrix(1,p1,p1)%x%w)*S)
  }else{
    stop("Axis must be 0,1,or 2!")
  }

}

