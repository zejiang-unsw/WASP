#--------------------------------------------------------------------------------
#' Variance Transformation Operation - MRA
#' @param data		  A list of response x and dependent variables dp.
#' @param wf		    Name of the wavelet filter to use in the decomposition.
#' @param J      	  Specifies the depth of the decomposition. This must be a number less than or equal to log(length(x),2).
#' @param method	  Either "dwt" or "modwt".
#' @param pad		    The method used for extend data to dyadic size. Use "per", "zero", or "sym".
#' @param boundary  Character string specifying the boundary condition. If boundary=="periodic" the default, then the vector you decompose is assumed to be periodic on its defined interval, if boundary=="reflection", the vector beyond its boundaries is assumed to be a symmetric reflection of itself.
#' @param cov.opt   Options of Covariance matrix sign. Use "pos", "neg", or "auto".
#'
#' @return A list of 8 elements: wf, method, boundary, pad, x (data), dp (data), dp.n (variance trasnformed dp), and S (covariance matrix).
#' @import waveslim
#' @export
#'
#' @references Z Jiang, A Sharma, and F Johnson. WRR
#' @references Percival, D. B. and A. T. Walden (2000) Wavelet Methods for Time Series Analysis, Cambridge University Press.
#'
#' @examples
#' data(rain.mon)
#' data(obs.mon)
#'
#' ## response SPI - calibration
#' SPI.cal <- SPI.calc(window(rain.mon, start=c(1949,1), end=c(1979,12)),sc=12)
#'
#' ## create paired response and predictors dataset for each station
#' data.list <- list()
#' for(id in 1:ncol(SPI.cal)){
#'   x <- window(SPI.cal[,id], start=c(1950,1), end=c(1979,12))
#'   dp <- window(obs.mon, start=c(1950,1), end=c(1979,12))
#'   data.list[[id]] <- list(x=as.numeric(x), dp=matrix(dp, nrow=nrow(dp)))
#' }
#'
#' ## variance transformation
#' dwt.list<- lapply(data.list, function(x) dwt.vt(x, wf="d4", J=7, method="dwt", pad="zero", boundary="periodic"))
#'
#' ## plot original and reconstrcuted predictors for each station
#' for(i in 1:length(dwt.list)){
#'   # extract data
#'   dwt <- dwt.list[[i]]
#'   x <- dwt$x  			      # response
#'   dp <- dwt$dp			      # original predictors
#'   dp.n <- dwt$dp.n       # variance transformed predictors
#'
#'   plot.ts(cbind(x,dp))
#'   plot.ts(cbind(x,dp.n))
#'
#' }

dwt.vt <- function(data, wf, J, method, pad, boundary, cov.opt=c("auto","pos","neg")){

  # initialization
  x= data$x; dp= as.matrix(data$dp)
  mu.dp <- apply(dp,2,mean)

  # variance transfrom
  ndim=ncol(dp);n=nrow(dp)
  S <- matrix(nrow=J+1,ncol=ndim)
  dp.n <- matrix(nrow=n,ncol=ndim)
  idwt.dp <- vector("list", ndim)

  for(i in 1:ndim){
    # center and padding
    dp.c <- scale(dp[,i],scale=F)
    dp.p <- padding(dp.c, pad=pad)

    # Multiresolution Analysis
    idwt.dp[[i]] <- waveslim::mra(dp.p, wf = wf, J = J, method = method, boundary = boundary)
    B <- matrix(unlist(lapply(idwt.dp[[i]], function(z) z[1:n])), ncol=J+1, byrow=FALSE)

    Bn <- scale(B)
    V <- as.numeric(apply(B,2,sd))

    dif <- sum(abs(Bn%*%V-dp.c))
    if(dif>10^-10) warning(paste0("Difference between Reconstructed and original:",dif))

    # variance transformation
    cov <- cov(x, Bn[1:length(x),])
    if(cov.opt=="pos") cov <- cov else if(cov.opt=="neg") cov <- -cov
    S[,i] <- as.vector(cov)

    Vr <- as.numeric(cov/norm(cov,type="2")*sd(dp.c))

    dp.n[,i] <- Bn%*%Vr + mu.dp[i]

    #check the correlation after vt then decide the direction of C
    if(cov.opt=="auto"){
      #if(cor(dp.n[,i],dp[,i])<0) cov <- -cov
      #cat(cor.test(dp.n[,i],dp[,i])$p.value,"&",cor.test(dp.n[,i],dp[,i])$estimate,"\n")
      if(cor.test(dp.n[,i],dp[,i])$estimate<0&cor.test(dp.n[,i],dp[,i])$p.value<0.05) cov <- -cov
      S[,i] <- as.vector(cov)

      Vr <- as.numeric(cov/norm(cov,type="2")*sd(dp.c))
      dp.n[,i] <- Bn%*%Vr + mu.dp[i]

    }

    dif.var <- (var(dp[,i])-var(dp.n[,i]))/var(dp[,i])
    if(dif.var>0.15) warning(paste0("Variance difference between Transformed and original(percentage):",dif.var*100))

  }

  dwt <- list(wavelet = wf,
              method = method,
              boundary = boundary,
              pad = pad,

              x=x,
              dp=dp,
              dp.n=dp.n,

              S=S
              )
  class(dwt) <- "dwt"

  return(dwt)

}

#--------------------------------------------------------------------------------
#' Variance Transformation Operation for Validation
#' @param data		  A list of response x and dependent variables dp.
#' @param J      	  Specifies the depth of the decomposition. This must be a number less than or equal to log(length(x),2).
#' @param dwt       A class of "dwt" data. Output from dwt.vt().
#'
#' @return          A list of 8 elements: wf, method, boundary, pad, x (data), dp (data), dp.n (variance trasnformed dp), and S (covariance matrix).
#' @export
#' @references Jiang, Z., Sharma, A., & Johnson, F. (2020). Refining Predictor Spectral Representation Using Wavelet Theory for Improved Natural System Modeling. Water Resources Research, 56(3), e2019WR026962. doi:10.1029/2019wr026962
#'
#' @examples
#' data(rain.mon)
#' data(obs.mon)
#'
#' ##response SPI - calibration
#' SPI.cal <- SPI.calc(window(rain.mon, start=c(1949,1), end=c(1979,12)),sc=12)
#'
#' ## create paired response and predictors dataset for each station
#' data.list <- list()
#' for(id in 1:ncol(SPI.cal)){
#'   x <- window(SPI.cal[,id], start=c(1950,1), end=c(1979,12))
#'   dp <- window(obs.mon, start=c(1950,1), end=c(1979,12))
#'   data.list[[id]] <- list(x=as.numeric(x), dp=matrix(dp, nrow=nrow(dp)))
#' }
#'
#' ## variance transformation - calibration
#' dwt.list<- lapply(data.list, function(x) dwt.vt(x, wf="d4", J=7, method="dwt", pad="zero", boundary="periodic"))
#'
#' ##response SPI - validation
#' SPI.val <- SPI.calc(window(rain.mon, start=c(1979,1), end=c(2009,12)),sc=12)
#'
#' ## create paired response and predictors dataset for each station
#' data.list <- list()
#' for(id in 1:ncol(SPI.val)){
#'   x <- window(SPI.val[,id], start=c(1980,1), end=c(2009,12))
#'   dp <- window(obs.mon, start=c(1980,1), end=c(2009,12))
#'   data.list[[id]] <- list(x=as.numeric(x), dp=matrix(dp, nrow=nrow(dp)))
#' }
#'
#' #variance transformation - validation
#' dwt.list.val<- lapply(1:length(data.list), function(i) dwt.vt.val(data.list[[i]], J=7, dwt.list[[i]]))
#'
#' ## plot original and reconstrcuted predictors for each station
#' for(i in 1:length(dwt.list.val)){
#'   # extract data
#'   dwt <- dwt.list.val[[i]]
#'   x <- dwt$x  			        # response
#'   dp <- dwt$dp			        # original predictors
#'   dp.n <- dwt$dp.n         # variance transformed predictors
#'
#'   plot.ts(cbind(x,dp))
#'   plot.ts(cbind(x,dp.n))
#'
#' }

dwt.vt.val <- function(data, J, dwt){

  # initialization
  x= data$x; dp= as.matrix(data$dp)
  wf <- dwt$wavelet; method <- dwt$method; boundary <- dwt$boundary; pad=dwt$pad
  mu.dp <- apply(dp,2,mean)

  # variance transfrom
  ndim=ncol(dp);n=nrow(dp)
  dp.n <- matrix(nrow=n,ncol=ndim)
  idwt.dp <- vector("list", ndim)

  for(i in 1:ndim){
    # center and padding
    dp.c <- scale(dp[,i],scale=F)
    dp.p <- padding(dp.c, pad=pad)

    # Multiresolution Analysis
    idwt.dp[[i]] <- waveslim::mra(dp.p, wf = wf, J = J, method = method, boundary = boundary)
    B <- matrix(unlist(lapply(idwt.dp[[i]], function(z) z[1:n])), ncol=J+1, byrow=FALSE)

    Bn <- scale(B)
    V <- as.numeric(apply(B,2,sd))

    dif <- sum(abs(Bn%*%V-dp.c))
    if(dif>10^-10) warning(paste0("Difference between Reconstructed and original:",dif))

    # in case different J
    cov <- rep(0,J+1)
    if(length(dwt$S[,i])>(J+1)) {
      cov <- dwt$S[,i][1:(J+1)]
    } else { cov[1:length(dwt$S[,i])] <- dwt$S[,i] }

    Vr <- as.numeric(cov/norm(cov,type="2")*sd(dp.c))

    dp.n[,i] <- Bn%*%Vr + mu.dp[i]

    dif.var <- (var(dp[,i])-var(dp.n[,i]))/var(dp[,i])
    if(dif.var>0.15) warning(paste0("Variance difference between Transformed and original(percentage):",dif.var*100))

  }

  dwt <- list(wavelet = wf,
              method = method,
              boundary = boundary,
              pad = pad,

              x=x,
              dp=dp,
              dp.n=dp.n,

              S=dwt$S
              )
  class(dwt) <- "dwt"

  return(dwt)

}

#-------------------------------------------------------------------------------
#' Padding data to dyadic sample size
#' @param x     A vector or time series containing the data be to decomposed.
#' @param pad   Method for padding, including periodic, zero and symetric padding.
#'
#' @return      A dyadic length (power of 2) vector or time series.
#' @export
#'
#' @examples
#' x <- rnorm(360)
#' x1 <- padding(x, pad="per")
#' x2 <- padding(x, pad="zero")
#' x3 <- padding(x, pad="sym")
#' ts.plot(cbind(x,x1,x2,x3), col=1:4)

padding <- function(x,pad=c("per","zero","sym")){
  n <- length(x)
  N <- 2^(ceiling(log(n, 2)))
  if(pad=="per")
    xx <- c(x, x)[1:N]
  else if(pad=="zero")
    xx <- c(x, rep(0, N - n))
  else
    xx <- c(x, rev(x))[1:N]
}
