#--------------------------------------------------------------------------------
#' Variance Transformation Operation - MODWT
#' @param data		  A list of response x and dependent variables dp.
#' @param wf		    Name of the wavelet filter to use in the decomposition.
#' @param J      	  Specifies the depth of the decomposition. This must be a number less than or equal to log(length(x),2).
#' @param pad		    The method used for extend data to dyadic size. Use "per", "zero", or "sym".
#' @param boundary  Character string specifying the boundary condition. If boundary=="periodic" the default, then the vector you decompose is assumed to be periodic on its defined interval, if boundary=="reflection", the vector beyond its boundaries is assumed to be a symmetric reflection of itself.
#' @param vt.opt    Options of variance transformation matrix, either "Sxx" or "Cov".
#' @param cov.opt   Options of Covariance matrix sign. Use "pos", "neg", or "auto".
#'
#' @return A list of 8 elements: wf, vt.opt, boundary, pad, x (data), dp (data), dp.n (variance trasnformed dp), and S (covariance matrix).
#' @import waveslim
#' @export
#'
#' @references Z Jiang, A Sharma, and F Johnson. WRR
#' @references Percival, D. B. and A. T. Walden (2000) Wavelet Methods for Time Series Analysis, Cambridge University Press.
#'
#' @examples
#' ###real-world example
#' x <- window(SPI.12,start=c(1950,1),end=c(2009,12))
#' dp <- window(obs.mon,start=c(1950,1),end=c(2009,12))
#' vt.opt = ifelse(0,"Sxx","Cov")
#'
#' for(id in 5){
#'
#'   data <- list(x=x[,id],dp=dp)
#'   dwt <- modwt.vt(data, wf="d4", J=7, pad="zero", boundary="periodic", vt.opt)
#'
#'   par(mfrow=c(ncol(dp),1),pty="m",mar=c(1,4,1,2))
#'   for(i in 1:ncol(dp))
#'     ts.plot(dwt$dp[,i],dwt$dp.n[,i],xlab=NA,col=c("black","red"),lwd=c(2,1))
#' }
#' ###synthetic example
#' #frequency, sampled from a given range
#' fd <- c(3,5,10,15,25,30,55,70,95)
#'
#' data.SW1 <- data.gen.SW(nobs=512,fp=25,fd=fd)
#' dwt.SW1 <- modwt.vt(data.SW1, wf="d4", J=7, pad="zero", boundary="periodic",vt.opt)
#'
#' x.modwt <- waveslim::modwt(dwt.SW1$x, wf = "d4", n.levels = 7,  boundary = "periodic")
#' dp.modwt <- waveslim::modwt(dwt.SW1$dp[,1], wf = "d4", n.levels = 7,  boundary = "periodic")
#' dp.vt.modwt <- waveslim::modwt(dwt.SW1$dp.n[,1], wf = "d4", n.levels = 7, boundary = "periodic")
#'
#' sum(sapply(dp.modwt,var)); var(dwt.SW1$dp[,1])
#' sum(sapply(dp.vt.modwt,var)); var(dwt.SW1$dp.n[,1])
#'
#' data <- rbind(sapply(dp.modwt,var)/sum(sapply(dp.modwt,var)),
#'               sapply(dp.vt.modwt,var)/sum(sapply(dp.vt.modwt,var)))
#'
#' par(mfrow=c(1,1))
#' bar <- barplot(data, beside = T, col=c("red","blue"))
#' lines(x = bar[2,], y = sapply(x.modwt,var)/sum(sapply(x.modwt,var)))
#' points(x = bar[2,], y = sapply(x.modwt,var)/sum(sapply(x.modwt,var)))

modwt.vt <- function(data, wf, J, pad, boundary, vt.opt=c("Sxx","Cov"), cov.opt=c("auto","pos","neg")){

  # initialization
  x= data$x; dp= data$dp
  mu.dp <- apply(dp,2,mean)

  # variance transfrom
  ndim=ncol(dp);n=nrow(dp)
  S <- matrix(nrow=J+1,ncol=ndim)
  dp.n <- matrix(nrow=n,ncol=ndim)
  modwt.dp <- vector("list", ndim)

  modwt.x <- waveslim::modwt(x, wf = wf, n.levels = J, boundary = boundary)
  Sxx <- unlist(lapply(modwt.x, sd))

  for(i in 1:ndim){
    #center
    dp.c <- scale(dp[,i], scale=F)

    # MODWT - variance decomposition
    modwt.dp[[i]] <- waveslim::modwt(dp.c, wf = wf, n.levels = J, boundary = boundary)

    B <- matrix(unlist(modwt.dp[[i]]), ncol=J+1, byrow=FALSE)

    Bn <- scale(B)
    V <- as.numeric(apply(B,2,sd))

    dif <- sum(abs(imodwt(modwt.dp[[i]])-dp.c))
    if(dif>10^-10) warning(paste0("Difference between Reconstructed and original:",dif))

    # variance transformation
    if(vt.opt!="Sxx") cov <- cov(x, Bn) else cov <- Sxx
    if(cov.opt=="pos") cov <- cov else if(cov.opt=="neg") cov <- -cov
    S[,i] <- as.vector(cov)

    Vr <- as.numeric(cov/norm(cov,type="2")*sd(dp.c))
    B.n <- sapply(1:length(Vr), function(i) Bn[,i]*Vr[i])

    B.n.ls <- class.modwt(B.n, J, wf, boundary)
    dp.n[,i] <- imodwt(B.n.ls) + mu.dp[i]

    # modwt.x <- waveslim::modwt(imodwt(B.n.ls), wf = wf, n.levels = J, boundary = boundary)
    # Sxx.dp.n <- unlist(lapply(modwt.x, sd))
    # barplot(cbind(Sxx, Vr, Sxx.dp.n), beside=TRUE, ylim=c(0,3))

    #check the correlation after vt then decide the direction of C
    if(cov.opt=="auto"){
      if(cor(dp.n[,i],dp[,i])<0) cov <- -cov
      #cat(cor.test(dp.n[,i],dp[,i])$p.value,"&",cor.test(dp.n[,i],dp[,i])$estimate,"\n")
      #if(cor.test(dp.n[,i],dp[,i])$estimate<0&cor.test(dp.n[,i],dp[,i])$p.value<0.05) cov <- -cov
      S[,i] <- as.vector(cov)

      Vr <- as.numeric(cov/norm(cov,type="2")*sd(dp.c))
      B.n <- sapply(1:length(Vr), function(i) Bn[,i]*Vr[i])

      B.n.ls <- class.modwt(B.n, J, wf, boundary)
      dp.n[,i] <- imodwt(B.n.ls) + mu.dp[i]

    }

    #cat(sum(unlist(apply(B,2,var))), sum(unlist(apply(B.n,2,var))), var(dp.c))

    #statistics check
    dif.var <- (var(dp.c)-sum(unlist(apply(B,2,var))))/var(dp.c)
    if(dif.var>0.15) warning(paste0("df.var.origin between Reconstructed and original(percentage):",dif.var*100))

    dif.var <- (var(dp.c)-sum(unlist(apply(B.n,2,var))))/var(dp.c)
    if(dif.var>0.15) warning(paste0("df.var.vt between Reconstructed and original(percentage):",dif.var*100))

  }

  dwt <- list(wavelet = wf,
              vt.opt = vt.opt,
              boundary = boundary,
              pad = pad,

              x=x,
              dp=dp,
              dp.n=dp.n,

              S=S
  )
  class(dwt) <- "modwt"

  return(dwt)

}

#--------------------------------------------------------------------------------
#' Variance Transformation Operation for Validation
#' @param data		  A list of response x and dependent variables dp.
#' @param J      	  Specifies the depth of the decomposition. This must be a number less than or equal to log(length(x),2).
#' @param dwt       A class of "dwt" data. Output from dwt.vt().
#'
#' @return          A list of 8 elements: wf, vt.opt, boundary, pad, x (data), dp (data), dp.n (variance trasnformed dp), and S (covariance matrix).
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
#' dwt.list<- lapply(data.list, function(x) modwt.vt(x, wf="d4", J=7, pad="zero", boundary="periodic"))
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
#' dwt.list.val<- lapply(1:length(data.list), function(i) modwt.vt.val(data.list[[i]], J=7, dwt.list[[i]]))
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

modwt.vt.val <- function(data, J, dwt){

  # initialization
  x= data$x; dp= data$dp
  mu.dp <- apply(dp,2,mean)
  wf <- dwt$wavelet; boundary <- dwt$boundary; vt.opt <- dwt$vt.opt; pad <- dwt$pad

  # variance transfrom
  ndim=ncol(dp);n=nrow(dp);
  dp.n <- matrix(nrow=n,ncol=ndim)
  modwt.dp <- vector("list", ndim)
  for(i in 1:ndim){
    #center
    dp.c <- scale(dp[,i], scale=F)

    # MODWT - variance decomposition
    modwt.dp[[i]] <- waveslim::modwt(dp.c, wf = wf, n.levels = J, boundary = boundary)

    B <- matrix(unlist(modwt.dp[[i]]), ncol=J+1, byrow=FALSE)

    Bn <- scale(B)
    V <- as.numeric(apply(B,2,sd))

    dif <- sum(abs(imodwt(modwt.dp[[i]])-dp.c))
    if(dif>10^-10) warning(paste0("Difference between Reconstructed and original:",dif))

    # in case different J
    cov <- rep(0,J+1)
    if(length(dwt$S[,i])>(J+1)) {
      cov <- dwt$S[,i][1:(J+1)]
    } else { cov[1:length(dwt$S[,i])] <- dwt$S[,i] }

    Vr <- as.numeric(cov/norm(cov,type="2")*sd(dp.c))
    B.n <- sapply(1:length(Vr), function(i) Bn[,i]*Vr[i])

    B.n.ls <- class.modwt(B.n, J, wf, boundary)
    dp.n[,i] <- imodwt(B.n.ls) + mu.dp[i]

    #statistic check
    dif.var <- (var(dp.c)-sum(unlist(apply(B,2,var))))/var(dp.c)
    if(dif.var>0.15) warning(paste0("df.var.origin between Reconstructed and original(percentage):",dif.var*100))

    dif.var <- (var(dp.c)-sum(unlist(apply(B.n,2,var))))/var(dp.c)
    if(dif.var>0.15) warning(paste0("df.var.vt between Reconstructed and original(percentage):",dif.var*100))


  }

  dwt <- list(wavelet = wf,
              vt.opt = vt.opt,
              boundary = boundary,
              pad = pad,

              x=x,
              dp=dp,
              dp.n=dp.n,

              S=dwt$S
              )
  class(dwt) <- "modwt"

  return(dwt)

}

#-------------------------------------------------------------------------------
#' Convert matrix to modwt class object
#'
#' @param B.modwt   A matrix of modwt decomposition.
#' @param J         Specifies the depth of the decomposition. This must be a number less than or equal to log(length(x),2).
#' @param wf        Name of the wavelet filter to use in the decomposition.
#' @param boundary  Character string specifying the boundary condition. If boundary=="periodic" the default, then the vector you decompose is assumed to be periodic on its defined interval, if boundary=="reflection", the vector beyond its boundaries is assumed to be a symmetric reflection of itself.
#'
#' @return          A class of "modwt" object
#' @export
#'

class.modwt <- function(B.modwt, J, wf, boundary){

  B.modwt.ls <- as.list(as.data.frame(B.modwt))
  names(B.modwt.ls) <- c(paste0("d", 1:J),paste0("s",J))

  attributes(B.modwt.ls)$class <- "modwt"
  attributes(B.modwt.ls)$wavelet <- wf
  attributes(B.modwt.ls)$boundary <- boundary

  return(B.modwt.ls)
}
