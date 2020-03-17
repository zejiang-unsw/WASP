#--------------------------------------------------------------------------------
#' Variance Transformation Operation - AT(a trous)
#' @param data		  A list of response x and dependent variables dp.
#' @param wf		    Name of the wavelet filter to use in the decomposition.
#' @param J      	  Specifies the depth of the decomposition. This must be a number less than or equal to log(length(x),2).
#' @param pad		    The method used for extend data to dyadic size. Use "per", "zero", or "sym".
#' @param boundary  Character string specifying the boundary condition. If boundary=="periodic" the default, then the vector you decompose is assumed to be periodic on its defined interval, if boundary=="reflection", the vector beyond its boundaries is assumed to be a symmetric reflection of itself.
#' @param cov.opt   Options of Covariance matrix sign. Use "pos", "neg", or "auto".
#'
#' @return A list of 8 elements: wf, J, boundary, pad, x (data), dp (data), dp.n (variance trasnformed dp), and S (covariance matrix).
#' @import wavethresh
#' @export
#'
#' @references Z Jiang, A Sharma, and F Johnson. WRR
#' @references Percival, D. B. and A. T. Walden (2000) Wavelet Methods for Time Series Analysis, Cambridge University Press.
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
#' ## variance transformation
#' dwt.list<- lapply(data.list, function(x) at.vt(x, wf="d4", pad="zero", boundary="periodic"))
#'
#' ## plot original and reconstrcuted predictors for each station
#' for(i in 1:length(dwt.list)){
#'   # extract data
#'   dwt <- dwt.list[[i]]
#'   x <- dwt$x  			      # response
#'   dp <- dwt$dp			      # original predictors
#'   dp.n <- dwt$dp.n        # variance transformed predictors
#'
#'   plot.ts(cbind(x,dp))
#'   plot.ts(cbind(x,dp.n))
#'
#' }

at.vt <- function(data, wf, J, pad, boundary, cov.opt=c("auto","pos","neg")){

  # initialization
  x= data$x; dp= data$dp
  mu.dp <- apply(dp,2,mean)
  if(wf!="haar") v <- as.integer(as.numeric(substr(wf,2,3))/2) else v <- 1

  # variance transfrom
  ndim=ncol(dp); n=nrow(dp); J <- ceiling(log2(n))
  S <- matrix(nrow=J+1, ncol=ndim)
  dp.n <- matrix(nrow=n,ncol=ndim)

  for(i in 1:ndim){

    # Padding
    dp.c <- scale(dp[,i],scale=F)
    pp <- padding(dp.c, pad=pad)

    # Multiresolution Analysis
    coefs <- at.wd(pp, v, nthresh=J, boundary=boundary)
    B <- coefs[1:n,]

    Bn <- scale(B)
    Bn[is.na(Bn)] <- 0 # for haar
    V <- as.numeric(apply(B,2,sd))

    dif <- sum(abs(Bn%*%V+mu.dp[i]-dp[,i]));
    if(dif>10^-10) warning(paste0("Difference between Reconstructed and original:",dif))

    cov <- cov(x, Bn)
    if(cov.opt=="pos") cov <- cov else if(cov.opt=="neg") cov <- -cov
    S[,i] <- as.vector(cov)

    #Vr <- as.numeric(cov/norm(cov,type="2")*sd(dp.c)) # For variance decompoistion
    Vr <- as.numeric(cov/norm(cov,type="2")*norm(V,type="2")) #For a trous!

    dp.n[,i] <- scale(Bn%*%Vr)*sd(dp[,i]) + mu.dp[i]

    #check the correlation after vt then decide the direction of C
    if(cov.opt=="auto"){
      #if(cor(dp.n[,i],dp[,i])<0) cov <- -cov
      #cat(cor.test(dp.n[,i],dp[,i])$p.value,"&",cor.test(dp.n[,i],dp[,i])$estimate,"\n")
      if(cor.test(dp.n[,i],dp[,i])$estimate<0&cor.test(dp.n[,i],dp[,i])$p.value<0.05) cov <- -cov
      S[,i] <- as.vector(cov)

      Vr <- as.numeric(cov/norm(cov,type="2")*norm(V,type="2")) #For a trous!
      dp.n[,i] <- scale(Bn%*%Vr)*sd(dp[,i]) + mu.dp[i]

    }

    dif.var <- (var(dp[,i])-var(dp.n[,i]))/var(dp[,i])
    if(dif.var>0.15) warning(paste0("Variance difference between Reconstructed and original(percentage):",dif.var*100))

  }

  dwt <- list(wavelet = wf,
              J = J,
              boundary = boundary,
              pad = pad,

              x=x,
              dp=dp,
              dp.n=dp.n,

              S=S
              )
  class(dwt) <- "atrous"

  return(dwt)

}

#--------------------------------------------------------------------------------
#' Variance Transformation Operation for Validation
#' @param data		  A list of response x and dependent variables dp.
#' @param J      	  Specifies the depth of the decomposition. This must be a number less than or equal to log(length(x),2).
#' @param dwt       A class of "at" data. Output from at.vt().
#'
#' @return A list of 8 elements: wf, J, boundary, pad, x (data), dp (data), dp.n (variance trasnformed dp), and S (covariance matrix).
#' @export
#' @references Z Jiang, A Sharma, and F Johnson. WRR
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
#' dwt.list<- lapply(data.list, function(x) at.vt(x, wf="d4", pad="zero", boundary="periodic"))
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
#' dwt.list.val<- lapply(1:length(data.list), function(i) at.vt.val(data.list[[i]], dwt.list[[i]]))
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

at.vt.val <- function(data, J, dwt){

  # initialization
  x= data$x; dp= data$dp
  wf <- dwt$wavelet; boundary <- dwt$boundary; pad=dwt$pad
  mu.dp <- apply(dp,2,mean)
  if(wf!="haar") v <- as.integer(as.numeric(substr(wf,2,3))/2) else v <- 1

  # variance transfrom
  ndim=ncol(dp); n=nrow(dp); #J <- ceiling(log2(n))
  dp.n <- matrix(nrow=n,ncol=ndim)

  for(i in 1:ndim){

    # Padding
    dp.c <- scale(dp[,i],scale=F)
    pp <- padding(dp.c, pad=pad)

    # Multiresolution Analysis
    coefs <- at.wd(pp, v, nthresh=J, boundary=boundary)
    B <- coefs[1:n,]

    Bn <- scale(B)
    Bn[is.na(Bn)] <- 0 # for haar
    V <- as.numeric(apply(B,2,sd))

    dif <- sum(abs(Bn%*%V+mu.dp[i]-dp[,i]))
    if(dif>10^-10) warning(paste0("Difference between Reconstructed and original:",dif))

    # in case different J
    cov <- rep(0,J+1)
    if(length(dwt$S[,i])>(J+1)) {
      cov <- dwt$S[,i][1:(J+1)]
    } else { cov[1:length(dwt$S[,i])] <- dwt$S[,i] }

    Vr <- as.numeric(cov/norm(cov,type="2")*sd(dp.c))

    dp.n[,i] <- scale(Bn%*%Vr)*sd(dp[,i]) + mu.dp[i]

    dif.var <- (var(dp[,i])-var(dp.n[,i]))/var(dp[,i])
    if(dif.var>0.15) warning(paste0("Variance difference between Reconstructed and original(percentage):",dif.var*100))

  }

  dwt <- list(wavelet = wf,
              J = J,
              boundary = boundary,
              pad = pad,

              x=x,
              dp=dp,
              dp.n=dp.n,

              S=dwt$S
              )
  class(dwt) <- "atrous"

  return(dwt)

}

#--------------------------------------------------------------------------------
#' a trous (AT) based additive decompostion using Daubechies family wavelet
#'
#' @param xx    The input time series with dyadic number size.
#' @param v     The number of vanishing moments for a given wavelet (e.g., haar(v=1) and d4 (v=2)).
#'
#' @return A matrix of decomposed sub-time series.
#' @export
#' @references Nason, G. P. (1996). Wavelet shrinkage using cross‐validation. Journal of the Royal Statistical Society: Series B (Methodological), 58(2), 463-479.
#'
#' @examples
#' data(SPI.12)
#' data(obs.mon)
#'
#' n <- nrow(SPI.12)
#' SPI <- padding(SPI.12[,1],pad="zero")
#' at.SPI <- at.wd(SPI, v=2)
#'
#' plot.ts(cbind(SPI[1:n],at.SPI[1:n,1:9]))
#' print(sum(abs(SPI[1:n]-rowSums(at.SPI[1:n,]))))
#'
#' names <- colnames(obs.mon)
#' at.atm <- vector("list", ncol(obs.mon))
#' for(i in 1:ncol(obs.mon)){
#'   tmp <- padding(obs.mon[,i], pad="zero")
#'   at.atm <- at.wd(tmp, v=2)
#'
#'   plot.ts(cbind(obs.mon[1:n,i],at.atm[1:n,1:9]), main=names[i])
#'   print(sum(abs(obs.mon[1:n,i]-rowSums(at.atm[1:n,]))))
#'
#' }

at.wd <- function(xx, v, nthresh, boundary,...){

    #DaubExPhase for Daubechies' extremal phase wavelets
    #DaubLeAsymm for Daubechies' “least-asymmetric” wavelets
    at.x <- wavethresh::wd(xx, type="station", filter.number=v, family="DaubExPhase", bc=boundary,...)

    # #nthresh = nlevelsWT(at.x)-1
    # Dj <- matrix(NA, nrow=length(xx), ncol=nthresh)
    # for (j in 1:(nthresh)) {
    #   Dj[,j] <- accessC(at.x,j) - accessC(at.x,j-1)
    # }
    #
    # at.wd <- cbind(accessC(at.x,0),Dj)
    # output <- at.wd[,ncol(at.wd):1]    #reverse the order

    max = nlevelsWT(at.x)
    Dj <- NULL
    for (j in (max-nthresh+1):max) {
      Dj <- cbind(Dj, accessC(at.x,j) - accessC(at.x,j-1))
    }

    at.wd <- cbind(accessC(at.x,max-nthresh),Dj)
    output <- at.wd[,ncol(at.wd):1]    #reverse the order

    return(output)
}


