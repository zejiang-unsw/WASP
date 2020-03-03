data.gen.ar1<-function(nobs,ndim=9){
	#nobs<-1000;
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

data.gen.ar4<-function(nobs,ndim=9){
	#nobs<-1000;
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

data.gen.ar9<-function(nobs,ndim=9){
  #nobs<-1000;
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

data.gen.tar1<-function(nobs,ndim=9, noise=0.1){
	#nobs<-1000;
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

# data.gen.tar2<-function(nobs,ndim=9,noise=0.1){
#   #nobs<-1000;
#   nwarm1=nwarm2=500
#   n=nobs+nwarm1+nwarm2
#   x<-matrix(0,n,1)
#   for (i in 1:nwarm1) {
#     x[i]<-rnorm(1,mean=0,sd=1)
#   }
#   dp<-matrix(0,(nobs),ndim)
#   for (i in (nwarm1+1):n){
#     eps<-rnorm(1,mean=0,sd=1)
#     xi2=x[i-2]
#     if(xi2>0) x[i]<- 0.6*x[i-1]-0.1*x[i-2]+noise*eps else x[i]= -1.1*x[i-1]+noise*eps
#   }
#   for(i in 1:ndim) dp[,i]=x[(n-i-nobs+1):(n-i)]
#   x=x[(n-nobs+1):n]
#   data_generated<-list(x=x,dp=dp,true.cpy=c(1,2))
#   return(data_generated)
# }	

data.gen.tar2<-function(nobs,ndim=9,noise=0.1){
  #nobs<-1000;
  nwarm1=nwarm2=50
  n=nobs+nwarm1+nwarm2
  x<-matrix(0,n,1)
  for (i in 1:nwarm1) {
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