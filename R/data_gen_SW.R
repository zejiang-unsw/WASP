#' Generate predictor and response data: Sinewave model
#'
#' @param nobs    The data length to be generated.
#' @param fp      The frequencies in the generated response.
#' @param fd      A vector of frequencies for potential predictors. fd = c(3,5,10,15,25,30,55,70,95) used in the WRR paper.
#' @param sd.x    The noise level in the predictor.
#' @param sd.y    The noise level in the response.
#'
#' @return A list of 3 elements: a vector of response (x), a matrix of potential predictors (dp) with each column containing one potential predictor, and a vector of true predictor numbers.
#' @export
#'
#' @examples
#' ###synthetic example
#' #frequency, sampled from a given range
#' fd <- c(3,5,10,15,25,30,55,70,95)
#'
#' data.SW1 <- data.gen.SW(nobs=512,fp=25,fd=fd)
#' data.SW3 <- data.gen.SW(nobs=512,fp=c(15,25,30),fd=fd)
#'
#' ts.plot(ts(data.SW1$x),ts(data.SW3$x),col=c("black","red"))
#' plot.ts(cbind(data.SW1$x,data.SW1$dp))
#' plot.ts(cbind(data.SW3$x,data.SW3$dp))

data.gen.SW<-function(nobs=512,fp=25,fd,sd.x=0.1,sd.y=0.1)
{

  t <- seq(0,1,length.out = nobs)

  index <- which(fd %in% fp)
  ndim <- length(fd)

  dp<-matrix(0,nobs,ndim)
  for(i in 1:ndim){
    eps<-rnorm(nobs,mean=0,sd=1)
    dp[,i]<- sin(2*pi*fd[i]*t) + sd.x*eps

  }

  x <- rowSums(sapply(1:length(index), function(i) sin(2*pi*fd[index[i]]*t))) + rnorm(nobs,0,sd.y)

  data_generated<-list(x=x,
                       dp=dp,
                       true.cpy=index)

  return(data_generated)
}
