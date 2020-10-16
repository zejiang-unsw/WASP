#--------------------------------------------------------------------------------
#' Variance Transformation Operation - AT(a trous)
#' @param data		  A list of response x and dependent variables dp.
#' @param wf		    Name of the wavelet filter to use in the decomposition.
#' @param J      	  Specifies the depth of the decomposition. This must be a number less than or equal to log(length(x),2).
#' @param boundary  Character string specifying the boundary condition. If boundary=="periodic" the default, then the vector you decompose is assumed to be periodic on its defined interval, if boundary=="reflection", the vector beyond its boundaries is assumed to be a symmetric reflection of itself.
#' @param cov.opt   Options of Covariance matrix sign. Use "pos", "neg", or "auto".
#' @param flag      Biased or Unbiased variance transformation, c("biased","unbiased").
#' @param detrend   Detrend the input time series or just center, default (F)
#'
#' @return A list of 8 elements: wf, J, boundary, x (data), dp (data), dp.n (variance transformed dp), and S (covariance matrix).
#' @import waveslim
#' @export
#'
#' @references Jiang, Z., Sharma, A., & Johnson, F. (2020). Refining Predictor Spectral Representation Using Wavelet Theory for Improved Natural System Modeling. Water Resources Research, 56(3), e2019WR026962. doi:10.1029/2019wr026962
#'
#' @examples
#' data(rain.mon)
#' data(obs.mon)
#'
#' ##response SPI - calibration
#' # SPI.cal <- SPI.calc(window(rain.mon, start=c(1949,1), end=c(1979,12)),sc=12)
#' SPI.cal <- SPEI::spi(window(rain.mon, start=c(1949,1), end=c(1979,12)),scale=12)$fitted
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
#' dwt.list<- lapply(data.list, function(x) at.vt(x, wf="d4", J=7, boundary="periodic", cov.opt="auto"))
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

at.vt <- function(data, wf, J, boundary, cov.opt=c("auto","pos","neg"), flag="biased", detrend=F)
{
  # initialization
  x= data$x; dp= as.matrix(data$dp)
  mu.dp <- apply(dp,2,mean)

  # variance transform
  ndim=ncol(dp); n=nrow(dp);
  S <- matrix(nrow=J+1, ncol=ndim)
  dp.n <- matrix(nrow=n,ncol=ndim)

  for(i in 1:ndim){
    # center or detrend
    if(!detrend){
      dp.c <- scale(dp[,i],scale=F)
    } else {
      #dp.c <- lm(dp[,i]~c(1:n))$resid
      dp.c <- dp[,i]-smooth.spline(1:n, dp[,i], spar=1)$y
    }

    # AT - additive decomposition
    at.dp <- at.wd(dp.c, wf=wf, J=J, boundary=boundary)
    B <- matrix(unlist(at.dp), ncol=J+1, byrow=FALSE)

    Bn <- scale(B)
    V <- as.numeric(apply(B,2,sd))

    dif <- sum(abs(Bn%*%V-dp.c))
    if(dif>10^-10) warning(paste0("Difference between Reconstructed and original:",dif))

    # variance transformation
    cov <- cov(x, Bn[1:length(x),])
    #cat("Biased: ", round(cov,3),"\n")

    if(flag=="unbiased"){ ###unbiased wavelet variance - only change cov
      at.dp.n <- non.bdy(at.dp, wf=wf, method="modwt")

      B.n <- matrix(unlist(at.dp.n), ncol=J+1, byrow=FALSE)
      cov <- cov(x, scale(B.n)[1:length(x),], use="pairwise.complete.obs")

      #if unbiased cov is not available then use biased
      cov[is.na(cov)] <- cov(x, Bn[1:length(x),])[is.na(cov)]
      #cat("Unbiased: ",round(cov,3),"\n")
    }

    if(cov.opt=="pos") cov <- cov else if(cov.opt=="neg") cov <- -cov
    S[,i] <- as.vector(cov)

    Vr <- as.numeric(cov/norm(cov,type="2")*sd(dp.c))
    if(!detrend){
      dp.n[,i] <- Bn%*%Vr + mu.dp[i]
    } else {
      #dp.n[,i] <- Bn%*%Vr + lm(dp[,i]~c(1:n))$fitted
      dp.n[,i] <- Bn%*%Vr + smooth.spline(1:n, dp[,i], spar=1)$y
    }

    #check the correlation after vt then decide the direction of C
    if(cov.opt=="auto"){

      if(cor.test(dp.n[,i],dp[,i])$estimate<0&cor.test(dp.n[,i],dp[,i])$p.value<0.05) cov <- -cov
      S[,i] <- as.vector(cov)

      Vr <- as.numeric(cov/norm(cov,type="2")*sd(dp.c))
      if(!detrend){
        dp.n[,i] <- Bn%*%Vr + mu.dp[i]
      } else {
        #dp.n[,i] <- Bn%*%Vr + lm(dp[,i]~c(1:n))$fitted
        dp.n[,i] <- Bn%*%Vr + smooth.spline(1:n, dp[,i], spar=1)$y
      }

    }

    #dif.var <- abs(var(dp[,i])-var(dp.n[,i]))/var(dp[,i])
    #if(dif.var>0.15) warning(paste0("Variance difference between Transformed and original(percentage):",dif.var*100))

  }

  dwt <- list(wavelet = wf,
              J = J,
              boundary = boundary,

              x=x,
              dp=dp,
              dp.n=dp.n,

              S=S
              )
  class(dwt) <- "at"

  return(dwt)

}

#--------------------------------------------------------------------------------
#' Variance Transformation Operation for Validation
#' @param data		  A list of response x and dependent variables dp.
#' @param J      	  Specifies the depth of the decomposition. This must be a number less than or equal to log(length(x),2).
#' @param dwt       A class of "at" data. Output from at.vt().
#' @param detrend   Detrend the input time series or just center, default (F)
#'
#' @return A list of 8 elements: wf, J, boundary, x (data), dp (data), dp.n (variance transformed dp), and S (covariance matrix).
#' @export
#' @references Jiang, Z., Sharma, A., & Johnson, F. (2020). Refining Predictor Spectral Representation Using Wavelet Theory for Improved Natural System Modeling. Water Resources Research, 56(3), e2019WR026962. doi:10.1029/2019wr026962
#'
#' @examples
#' data(rain.mon)
#' data(obs.mon)
#'
#' ##response SPI - calibration
#' # SPI.cal <- SPI.calc(window(rain.mon, start=c(1949,1), end=c(1979,12)),sc=12)
#' SPI.cal <- SPEI::spi(window(rain.mon, start=c(1949,1), end=c(1979,12)),scale=12)$fitted
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
#' dwt.list<- lapply(data.list, function(x) at.vt(x, wf="d4", J=7, boundary="periodic", cov.opt="auto"))
#'
#' ##response SPI - validation
#' #SPI.val <- SPI.calc(window(rain.mon, start=c(1979,1), end=c(2009,12)),sc=12)
#' SPI.val <- SPEI::spi(window(rain.mon, start=c(1979,1), end=c(2009,12)),scale=12)$fitted
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
#' dwt.list.val<- lapply(1:length(data.list), function(i) at.vt.val(data.list[[i]], J=7, dwt.list[[i]]))
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

at.vt.val <- function(data, J, dwt, detrend=F){

  # initialization
  x= data$x; dp= as.matrix(data$dp)
  wf <- dwt$wavelet; boundary <- dwt$boundary
  mu.dp <- apply(dp,2,mean)

  # variance transform
  ndim=ncol(dp); n=nrow(dp)
  dp.n <- matrix(nrow=n,ncol=ndim)

  for(i in 1:ndim){
    # center or detrend
    if(!detrend){
      dp.c <- scale(dp[,i],scale=F)
    } else {
      #dp.c <- lm(dp[,i]~c(1:n))$resid
      dp.c <- dp[,i]-smooth.spline(1:n, dp[,i], spar=1)$y
    }

    # AT - additive decomposition
    at.dp <- at.wd(dp.c, wf=wf, J=J, boundary=boundary)
    B <- matrix(unlist(at.dp), ncol=J+1, byrow=FALSE)

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

    if(!detrend){
      dp.n[,i] <- Bn%*%Vr + mu.dp[i]
    } else {
      #dp.n[,i] <- Bn%*%Vr + lm(dp[,i]~c(1:n))$fitted
      dp.n[,i] <- Bn%*%Vr + smooth.spline(1:n, dp[,i], spar=1)$y
    }

    #dif.var <- abs(var(dp[,i])-var(dp.n[,i]))/var(dp[,i])
    #if(dif.var>0.15) warning(paste0("Variance difference between Transformed and original(percentage):",dif.var*100))

  }

  dwt <- list(wavelet = wf,
              J = J,
              boundary = boundary,

              x=x,
              dp=dp,
              dp.n=dp.n,

              S=dwt$S
              )
  class(dwt) <- "at"

  return(dwt)

}

#--------------------------------------------------------------------------------
#' a trous (AT) based additive decompostion using Daubechies family wavelet
#'
#' @param x         The input time series.
#' @param wf		    Name of the wavelet filter to use in the decomposition.
#' @param J      	  Specifies the depth of the decomposition. This must be a number less than or equal to log(length(x),2).
#' @param boundary  Character string specifying the boundary condition. If boundary=="periodic" the default, then the vector you decompose is assumed to be periodic on its defined interval, if boundary=="reflection", the vector beyond its boundaries is assumed to be a symmetric reflection of itself.
#'
#' @return A matrix of decomposed sub-time series.
#' @export
#' @references Nason, G. P. (1996). Wavelet shrinkage using cross‐validation. Journal of the Royal Statistical Society: Series B (Methodological), 58(2), 463-479.
#'
#' @examples
#' data(obs.mon)
#'
#' n <- nrow(obs.mon); v=1
#' J <- floor(log(n/(2*v-1))/log(2)) #(Kaiser, 1994)
#'
#' names <- colnames(obs.mon)
#' at.atm <- vector("list", ncol(obs.mon))
#' for(i in 1:ncol(obs.mon)){
#'   tmp <- as.numeric(scale(obs.mon[,i],scale=FALSE))
#'   at.atm <- do.call(cbind, at.wd(tmp, wf="haar", J = J, boundary = "periodic"))
#'
#'   plot.ts(cbind(obs.mon[1:n,i],at.atm[1:n,1:9]), main=names[i])
#'   print(sum(abs(scale(obs.mon[1:n,i],scale=FALSE)-rowSums(at.atm[1:n,]))))
#'
#' }
at.wd <- function(x, wf, J, boundary="periodic"){
  s <- NULL
  for(i in 1:J){
    s<- cbind(s, waveslim::modwt(x, wf=wf, n.levels = i, boundary = boundary)[[i+1]])
  }
  s <- s[1:length(x),]

  at <- x-s[,1]
  for(i in 1:(J-1)) at <- cbind(at, s[,i]-s[,i+1])
  at <- cbind(at, s[,J])

  at.df <- as.data.frame(at)
  colnames(at.df) <- c(paste0("d",1:J), paste0("s",J))

  out <- as.list(at.df)
  attr(out, "class") <- "at"
  attr(out, "wavelet") <- wf
  attr(out, "boundary") <- boundary
  return(out)

}


