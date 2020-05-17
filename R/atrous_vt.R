#--------------------------------------------------------------------------------
#' Variance Transformation Operation - AT(a trous)
#' @param data		  A list of response x and dependent variables dp.
#' @param wf		    Name of the wavelet filter to use in the decomposition.
#' @param J      	  Specifies the depth of the decomposition. This must be a number less than or equal to log(length(x),2).
#' @param boundary  Character string specifying the boundary condition. If boundary=="periodic" the default, then the vector you decompose is assumed to be periodic on its defined interval, if boundary=="reflection", the vector beyond its boundaries is assumed to be a symmetric reflection of itself.
#' @param cov.opt   Options of Covariance matrix sign. Use "pos", "neg", or "auto".
#'
#' @return A list of 8 elements: wf, J, boundary, x (data), dp (data), dp.n (variance trasnformed dp), and S (covariance matrix).
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
#' dwt.list<- lapply(data.list, function(x) at.vt(x, wf="d4", J=7, boundary="periodic"))
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

at.vt <- function(data, wf, J, boundary, cov.opt=c("auto","pos","neg"), flag=c("biased","unbiased")){

  # initialization
  x= data$x; dp= as.matrix(data$dp)
  mu.dp <- apply(dp,2,mean)

  # variance transfrom
  ndim=ncol(dp); n=nrow(dp);
  S <- matrix(nrow=J+1, ncol=ndim)
  dp.n <- matrix(nrow=n,ncol=ndim)

  for(i in 1:ndim){

    # center
    dp.c <- scale(dp[,i],scale=F)

    # AT - additive decomposition
    at.dp <- at.wd(dp.c, wf=wf, J=J, boundary=boundary)
    B <- matrix(unlist(at.dp), ncol=J+1, byrow=FALSE)

    Bn <- scale(B)
    V <- as.numeric(apply(B,2,sd))

    dif <- sum(abs(Bn%*%V-dp.c))
    if(dif>10^-10) warning(paste0("Difference between Reconstructed and original:",dif))

    # variance transformation
    cov <- cov(x, Bn[1:length(x),])
    cat("Biased: ", round(cov,3),"\n")

    if(flag=="unbiased"){ ###unbiased wavelet variance - only change cov
      at.dp.n <- non.bdy(at.dp, wf=wf, method="modwt")

      B.n <- matrix(unlist(at.dp.n), ncol=J+1, byrow=FALSE)
      cov <- cov(x, scale(B.n)[1:length(x),], use="pairwise.complete.obs")
      cat("Unbiased: ",round(cov,3),"\n")

      #cov[is.na(cov)] <- 0
    }

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
  class(dwt) <- "atrous"

  return(dwt)

}

#--------------------------------------------------------------------------------
#' Variance Transformation Operation for Validation
#' @param data		  A list of response x and dependent variables dp.
#' @param J      	  Specifies the depth of the decomposition. This must be a number less than or equal to log(length(x),2).
#' @param dwt       A class of "at" data. Output from at.vt().
#'
#' @return A list of 8 elements: wf, J, boundary, x (data), dp (data), dp.n (variance trasnformed dp), and S (covariance matrix).
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
#' dwt.list<- lapply(data.list, function(x) at.vt(x, wf="d4", J=7, boundary="periodic"))
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

at.vt.val <- function(data, J, dwt){

  # initialization
  x= data$x; dp= as.matrix(data$dp)
  wf <- dwt$wavelet; boundary <- dwt$boundary
  mu.dp <- apply(dp,2,mean)

  # variance transfrom
  ndim=ncol(dp); n=nrow(dp);
  dp.n <- matrix(nrow=n,ncol=ndim)

  for(i in 1:ndim){
    # center
    dp.c <- scale(dp[,i],scale=F)

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

    dp.n[,i] <- Bn%*%Vr + mu.dp[i]

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
  class(dwt) <- "atrous"

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
#' n <- nrow(obs.mon);v=2
#' J <- ceiling(log(n/(2*v-1))/log(2)) #(Kaiser, 1994)
#'
#' names <- colnames(obs.mon)
#' at.atm <- vector("list", ncol(obs.mon))
#' for(i in 1:ncol(obs.mon)){
#'   tmp <- padding(scale(obs.mon[,i],scale=F), pad="zero")
#'   at.atm <- at.wd(tmp, v=2, nthresh = J, boundary = "periodic")
#'
#'   plot.ts(cbind(obs.mon[1:n,i],at.atm[1:n,1:9]), main=names[i])
#'   print(sum(abs(scale(obs.mon[1:n,i],scale=F)-rowSums(at.atm[1:n,]))))
#'
#' }
at.wd <- function(x, wf, J, boundary="periodic"){
  s <- NULL
  for(i in 1:J){
    s<- cbind(s, waveslim::modwt(x, wf=wf, n.levels = i, boundary = boundary)[[i+1]])
  }

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


