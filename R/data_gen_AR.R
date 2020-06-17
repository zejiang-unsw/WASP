# An autoregressive model of order p, AR(p), can be written as
# yt=c+ϕ1y_t−1+ϕ2y_t−2+⋯+ϕpy_t−p+εt,
# where εt is white noise.
# For an AR(1) model:
# •	when ϕ1=0, yt is equivalent to white noise;
# •	when ϕ1=1 and c=0, yt is equivalent to a random walk;
# •	when ϕ1=1 and c≠0, yt is equivalent to a random walk with drift;
# •	when ϕ1<0, yt tends to oscillate around the mean.
# We normally restrict autoregressive models to stationary data, in which case some constraints on the values of the parameters are required.
# •	For an AR(1) model: −1<ϕ1<1.
# •	For an AR(2) model: −1<ϕ2<1, ϕ1+ϕ2<1, ϕ2−ϕ1<1.
# When p≥3, the restrictions are much more complicated. R takes care of these restrictions when estimating a model.

#' Generate predictor and response data from AR1 model.
#'
#' @param nobs The data length to be generated.
#' @param ndim The number of potential predictors (default is 9).
#'
#' @return A list of 2 elements: a vector of response (x), and a matrix of potential predictors (dp) with each column containing one potential predictor.
#' @export
#'
#' @examples
#' # AR1 model from paper with 9 dummy variables
#' data.ar1<-data.gen.ar1(500)
#' plot.ts(cbind(data.ar1$x,data.ar1$dp))

data.gen.ar1<-function(nobs,ndim=9)
{
  nwarm1=nwarm2=50
  n=nobs+nwarm1+nwarm2
	x<-matrix(0,n,1)
	for (i in 1:nwarm1) {
		x[i]<-rnorm(1,mean=0,sd=1)
	}
	dp<-matrix(0,(nobs),ndim)
	for (i in (nwarm1+1):n){
		eps<-rnorm(1,mean=0,sd=1)
		x[i]<-0.9*x[i-1]+0.866*eps
	}
	for(i in 1:ndim) dp[,i]=x[(n-i-nobs+1):(n-i)]
	x=x[(n-nobs+1):n]
	data_generated<-list(x=x,dp=dp,true.cpy=c(1))
	return(data_generated)
}

#' Generate predictor and response data from AR4 model.
#'
#' @param nobs The data length to be generated.
#' @param ndim The number of potential predictors (default is 9).
#'
#' @return A list of 2 elements: a vector of response (x), and a matrix of potential predictors (dp) with each column containing one potential predictor.
#' @export
#'
#' @examples
#' # AR4 model from paper with total 9 dimensions
#' data.ar4<-data.gen.ar4(500)
#' plot.ts(cbind(data.ar4$x,data.ar4$dp))

data.gen.ar4<-function(nobs,ndim=9)
{
  nwarm1=nwarm2=50
  n=nobs+nwarm1+nwarm2
	x<-matrix(0,n,1)
	for (i in 1:nwarm1) {
		x[i]<-rnorm(1,mean=0,sd=1)
	}
	dp<-matrix(0,(nobs),ndim)
	for (i in (nwarm1+1):n){
		eps<-rnorm(1,mean=0,sd=1)
		x[i]<-0.6*x[i-1]-0.4*x[i-4]+eps
	}
	for(i in 1:ndim) dp[,i]=x[(n-i-nobs+1):(n-i)]
	x=x[(n-nobs+1):n]
	data_generated<-list(x=x,dp=dp,true.cpy=c(1,4))
	return(data_generated)
}

#' Generate predictor and response data from AR9 model.
#'
#' @param nobs The data length to be generated.
#' @param ndim The number of potential predictors (default is 9).
#'
#' @return A list of 2 elements: a vector of response (x), and a matrix of potential predictors (dp) with each column containing one potential predictor.
#' @export
#'
#' @examples
#' # AR9 model from paper with total 9 dimensions
#' data.ar9<-data.gen.ar9(500)
#' plot.ts(cbind(data.ar9$x,data.ar9$dp))

data.gen.ar9<-function(nobs,ndim=9)
{
  nwarm1=nwarm2=50
  n=nobs+nwarm1+nwarm2
  x<-matrix(0,n,1)
  for (i in 1:nwarm1) {
    x[i]<-rnorm(1,mean=0,sd=1)
  }
  dp<-matrix(0,(nobs),ndim)
  for (i in (nwarm1+1):n){
    eps<-rnorm(1,mean=0,sd=1)
    x[i]<-0.3*x[i-1]-0.6*x[i-4]-0.5*x[i-9]+eps
  }
  for(i in 1:ndim) dp[,i]=x[(n-i-nobs+1):(n-i)]
  x=x[(n-nobs+1):n]
  data_generated<-list(x=x,dp=dp,true.cpy=c(4,9,1))
  return(data_generated)
}

#' Generate predictor and response data from TAR1 model.
#'
#' @param nobs  The data length to be generated.
#' @param ndim  The number of potential predictors (default is 9).
#' @param noise The white noise in the data
#'
#' @return A list of 2 elements: a vector of response (x), and a matrix of potential predictors (dp) with each column containing one potential predictor.
#' @export
#'
#' @references Sharma, A. (2000). Seasonal to interannual rainfall probabilistic forecasts for improved water supply management: Part 1—A strategy for system predictor identification. Journal of Hydrology, 239(1-4), 232-239.
#'
#' @examples
#' # TAR1 model from paper with total 9 dimensions
#' data.tar1<-data.gen.tar1(500)
#' plot.ts(cbind(data.tar1$x,data.tar1$dp))

data.gen.tar1<-function(nobs,ndim=9, noise=0.1)
{
  nwarm1=nwarm2=50
  n=nobs+nwarm1+nwarm2
	x<-matrix(0,n,1)
	for (i in 1:nwarm1) {
		x[i]<-rnorm(1,mean=0,sd=1)
	}
	dp<-matrix(0,(nobs),ndim)
	for (i in (nwarm1+1):n){
		eps<-rnorm(1,mean=0,sd=1)
		xi3=x[i-3]
		if(xi3<=0) x[i]<- -0.9*x[i-3]+noise*eps else x[i]=0.4*x[i-3]+noise*eps
	}
	for(i in 1:ndim) dp[,i]=x[(n-i-nobs+1):(n-i)]
	x=x[(n-nobs+1):n]
	data_generated<-list(x=x,dp=dp,true.cpy=c(3))
	return(data_generated)
}

#' Generate predictor and response data from TAR2 model.
#'
#' @param nobs  The data length to be generated.
#' @param ndim  The number of potential predictors (default is 9).
#' @param noise The white noise in the data
#'
#' @return A list of 2 elements: a vector of response (x), and a matrix of potential predictors (dp) with each column containing one potential predictor.
#' @export
#'
#' @references Sharma, A. (2000). Seasonal to interannual rainfall probabilistic forecasts for improved water supply management: Part 1—A strategy for system predictor identification. Journal of Hydrology, 239(1-4), 232-239.
#'
#' @examples
#' # TAR2 model from paper with total 9 dimensions
#' data.tar2<-data.gen.tar2(500)
#' plot.ts(cbind(data.tar2$x,data.tar2$dp))

data.gen.tar2<-function(nobs,ndim=9,noise=0.1)
{
  nwarm1=nwarm2=50
  n=nobs+nwarm1+nwarm2
  x<-matrix(0,n,1)
  for (i in 1:nwarm1){
    x[i]<-rnorm(1,mean=0,sd=1)
  }
  dp<-matrix(0,(nobs),ndim)
  for (i in (nwarm1+1):n){
    eps<-rnorm(1,mean=0,sd=1)
    xi6=x[i-6]
    if(xi6<=0) x[i]<- -0.5*x[i-6]+0.5*x[i-10]+noise*eps else x[i]= 0.8*x[i-10]+noise*eps
  }
  for(i in 1:ndim) dp[,i]=x[(n-i-nobs+1):(n-i)]
  x=x[(n-nobs+1):n]
  data_generated<-list(x=x,dp=dp,true.cpy=c(10,6))
  return(data_generated)
}


