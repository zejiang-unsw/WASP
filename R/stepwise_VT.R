#' Calculate stepwise high order VT in calibration
#'
#' @param data    A list of data, including response and predictors
#' @param alpha   The significance level used to judge whether the sample estimate in Equation \deqn{\hat{PIC} = sqrt(1-exp(-2\hat{PI})} is significant or not. A default alpha value is 0.1.
#' @param mode    A mode of variance transfomration, i.e., MRA, MODWT, or AT
#' @param wf      Wavelet family
#'
#' @return A list of 2 elements: the column numbers of the meaningful predictors (cpy), and partial informational correlation (cpyPIC).
#' @export
#'
#' @references Sharma, A., Mehrotra, R., 2014. An information theoretic alternative to model a natural system using observational information alone. Water Resources Research, 50(1): 650-660.
#' @examples
#' ###Real-world example
#' mode <- switch(1,"MRA", "MODWT","a trous")
#' wf="d4"
#' station.id = 5 # station to investigate
#' SPI.12 <- SPEI::spi(rain.mon,scale=12)$fitted
#' lab.names <- colnames(obs.mon)
#' #plot.ts(SPI.12[,1:10])
#'
#' x <- window(SPI.12[,station.id],start=c(1950,1),end=c(1979,12))
#' dp <- window(obs.mon[,lab.names],start=c(1950,1),end=c(1979,12))
#'
#' data <- list(x=x,dp=matrix(dp, ncol=ncol(dp)))
#'
#' dwt = stepwise.VT(data, mode=mode, wf=wf)
#'
#' ###plot transformed predictor before and after
#' par(mfrow=c(ncol(dp),1), mar=c(0,3,2,1))
#' for(i in 1:ncol(dp))
#' {
#'   ts.plot(cbind(dwt$dp[,i], dwt$dp.n[,i]), xlab="NA", col=1:2)
#' }
stepwise.VT <- function (data, alpha = 0.1, mode=c("MRA","MODWT","AT"), wf)
{
  x = as.matrix(data$x)
  py= as.matrix(data$dp)

  n = nrow(x)
  npy = ncol(py)
  cpy = cpyPIC = NULL
  icpy = 0
  z = NULL;z.vt=NULL
  S=NULL
  isig = T
  icoloutz = 1:npy
  #cat("calc.PIC-----------","\n")
  while (isig) {
    npicmax = npy - icpy
    pictemp = rep(0, npicmax)
    y = py[, icoloutz]

    temp = pic.calc(x,y,z, mode, wf)
    pictemp = temp$pic

    pytmp = temp$py
    Stmp = temp$S

    ctmp = order(-pictemp)[1]
    cpytmp = icoloutz[ctmp]
    picmaxtmp = pictemp[ctmp]
    if (!is.null(z)) {
      z = as.matrix(z)
      df = n - ncol(z)
    } else {
      df = n
    }

    t <-  qt(1-alpha, df=df)
    picthres <- sqrt(t^2/(t^2+df))
    #cat("picthres",picthres,"\n")

    if (picmaxtmp > picthres) {
      cpy = c(cpy, cpytmp)
      cpyPIC = c(cpyPIC, picmaxtmp)
      z = cbind(z, py[, cpytmp])

      S = cbind(S, Stmp[,ctmp])
      z.vt = cbind(z.vt, pytmp[,ctmp])

      icoloutz = icoloutz[-ctmp]
      icpy = icpy + 1
      if ((npy - icpy) == 0) isig = F
    } else {
      isig = F
    }
  }
  #cat("calc.PW------------","\n")
  if (!is.null(z)) {
    out = pw.calc(x, z.vt, cpyPIC)

    outwt = out$pw
    lstwt = abs(lsfit(z.vt[1:length(x),], x)$coef)

    z.n = z; ncpy=length(cpy)
    if(ncpy>1){
      for(i in 2:ncpy){
        tmp=z[,1:(i-1)]
        z.n[,i] <- z[,i]-knnregl1cv(z[,i], tmp)
        #z.n[,i] <- lm.fit(as.matrix(tmp), z[,i])$residuals

      }
    }

    return(list(cpy = cpy, cpyPIC = cpyPIC, wt = outwt,lstwet = lstwt,
                x = data$x, py = py,
                dp=z.n, dp.n=z.vt, S=S,
                wavelet=wf))
  } else {
    message("None of the provided predictors is related to the response variable")
  }
}

#' Calculate stepwise high order VT in validation
#'
#' @param data    A list of data, including response and predictors
#' @param dwt     Output from dwt.vt(), including the transformation covariance
#' @param mode    A mode of variance transfomration, i.e., MRA, MODWT, or AT
#'
#' @return        A list of objects, including transformed predictors
#' @export
#'
#' @examples
#' ###Real-world example
#' mode <- switch(1,"MRA", "MODWT","a trous")
#' wf="d4"
#' station.id = 5 # station to investigate
#' SPI.12 <- SPEI::spi(rain.mon,scale=12)$fitted
#' lab.names <- colnames(obs.mon)
#' #plot.ts(SPI.12[,1:10])
#'
#' #--------------------------------------
#' ###calibration
#' x <- window(SPI.12[,station.id],start=c(1950,1),end=c(1979,12))
#' dp <- window(obs.mon[,lab.names],start=c(1950,1),end=c(1979,12))
#'
#' data <- list(x=x,dp=matrix(dp, ncol=ncol(dp)))
#' dwt = stepwise.VT(data, mode=mode, wf=wf)
#' #--------------------------------------
#' ###validation
#' x <- window(SPI.12[,station.id],start=c(1980,1),end=c(2009,12))
#' dp <- window(obs.mon[,lab.names],start=c(1980,1),end=c(2009,12))
#'
#' data.n <- list(x=x,dp=matrix(dp, ncol=ncol(dp)))
#' dwt.val = stepwise.VT.val(data.n, dwt, mode)
#'
#' ###plot transformed predictor before and after
#' par(mfrow=c(ncol(dp),1), mar=c(0,3,2,1))
#' for(i in 1:ncol(dp))
#' {
#'   ts.plot(cbind(dwt.val$dp[,i], dwt.val$dp.n[,i]), xlab="NA", col=1:2)
#' }
stepwise.VT.val <- function (data, dwt, mode){

  # initialization
  x= data$x; py= as.matrix(data$dp)
  cpy=dwt$cpy; ncpy=length(cpy)
  wf=dwt$wavelet; method <- "dwt"; boundary <- "periodic"; pad="zero"

  if(wf!="haar") v <- as.integer(readr::parse_number(wf)/2) else v <- 1
  #Maximum decomposition level J
  n <- length(x)
  J <- ceiling(log(n/(2*v-1))/log(2)) - 1 #(Kaiser, 1994)

  dwt.n = c(dwt, method=method, boundary=boundary, pad=pad)
  dp.n = py[,cpy]
  if(ncpy>1){
    for(i in 2:ncpy){
      Z=py[,cpy[1:(i-1)]]
      dp.n[,i] <- py[,cpy[i]]-knnregl1cv(py[,cpy[i]], Z)
      #dp.n[,i] <- lm.fit(as.matrix(Z), py[,cpy[i]])$residuals

    }
  }

  data.n = list(x=x, dp=dp.n)

  #variance transform
  if(mode=="MRA"){
    dwt.val<- dwt.vt.val(data.n, J, dwt.n)
  } else if(mode=="MODWT"){
    dwt.val<- modwt.vt.val(data.n, J, dwt.n)
  } else {
    dwt.val<- at.vt.val(data.n, J, dwt.n)
  }

  return(c(dwt.val,py=list(py)))

}

# Calculate the ratio of conditional error standard deviations
#
# @param x     A vector of response.
# @param zin   A matrix containing the meaningful predictors selected from a large set of possible predictors (z).
# @param zout  A matrix containing the remaining possible predictors after taking out the meaningful predictors (zin).
#
# @return The STD ratio.
#
# @references Sharma, A., Mehrotra, R., 2014. An information theoretic alternative to model a natural system using observational information alone. Water Resources Research, 50(1): 650-660.
calc.scaleSTDratio <- function (x, zin, zout)
{
  if(!missing(zout)){
    zout = as.matrix(zout)
    xhat = knnregl1cv(x, zout[1:length(x),])
    stdratxzout = sqrt(var(x - xhat)/var(x))
    zinhat = knnregl1cv(zin, zout)
    stdratzinzout = sqrt(var(zin - zinhat)/var(zin))
    return(0.5 * (stdratxzout + stdratzinzout))
  } else {
    return(1)
  }

}

#-------------------------------------------------------------------------------
kernel.est.uvn <- function(Z) {

  N <- length(Z)
  d <- 1
  # compute sigma & constant
  sigma <- 1.5*bw.nrd0(Z)
  constant <- sqrt(2*pi) * sigma * N

  # Commence main loop
  dens <- vector()
  for(h in 1:N) {
    dis.Z <- (Z - Z[h])^2
    exp.dis <- exp(-dis.Z / (2*sigma^2))
    dens[h] <- sum(exp.dis) / constant
  }
  return(dens)
}
#-------------------------------------------------------------------------------
kernel.est.mvn <- function(Z) {

  # Compute covariance and determinant of cov
  N <- nrow(Z)
  d <- ncol(Z)
  Cov <- cov(Z)
  det.Cov <- det(Cov)

  # Compute sigma & constant
  sigma <- 1.5 * (4/(d + 2))^(1/(d + 4)) * N^(-1/(d + 4))
  constant <- (sqrt(2*pi)*sigma)^d * sqrt(det.Cov) * N

  # Commence main loop
  dens <- vector()
  for(h in 1:N) {
    dist.val <- mahalanobis(Z, center = Z[h,], cov = Cov)
    exp.dis <- exp(-dist.val / (2*sigma^2))
    dens[h] <- sum(exp.dis) / constant
  }

  return(dens)
}
#-------------------------------------------------------------------------------
pmi.calc <- function(X, Y) {

  N <- length(X)

  pdf.X <- kernel.est.uvn(X)
  pdf.Y <- kernel.est.uvn(Y)
  pdf.XY <- kernel.est.mvn(cbind(X, Y))

  calc <- log(pdf.XY / (pdf.X * pdf.Y))
  return(sum(calc) / N)

}
#-------------------------------------------------------------------------------
# Calculate PIC
#
# @param X       A vector of response.
# @param Y       A matrix of new predictors.
# @param Z       A matrix of pre-existing predictors that could be NULL if no prior predictors exist.
# @param mode    A mode of variance transfomration, i.e., MRA, MODWT, or AT
# @param wf      Wavelet family
#
# @return A list of 2 elements: the partial mutual information (pmi), and partial informational correlation (pic).
# @export
#
# @references Sharma, A., Mehrotra, R., 2014. An information theoretic alternative to model a natural system using observational information alone. Water Resources Research, 50(1): 650-660.
# @references Galelli S., Humphrey G.B., Maier H.R., Castelletti A., Dandy G.C. and Gibbs M.S. (2014) An evaluation framework for input variable selection algorithms for environmental data-driven models, Environmental Modelling and Software, 62, 33-51, DOI: 10.1016/j.envsoft.2014.08.015.
pic.calc <- function(X, Y, Z, mode, wf) {

  Y=as.matrix(Y)

  if(wf!="haar") v <- as.integer(readr::parse_number(wf)/2) else v <- 1
  #Maximum decomposition level J
  n <- length(X)
  J <- ceiling(log(n/(2*v-1))/log(2)) - 1 #(Kaiser, 1994)


  if(is.null(Z)){
    x.in <- X
    y.in <- Y
  } else {
    Z=as.matrix(Z)
    x.in <- X-knnregl1cv(X, Z[1:length(X),])
    y.in <- apply(Y, 2, function(i) i-knnregl1cv(i, Z))

    # x.in <- lm.fit(as.matrix(Z[1:length(X),]), X)$residuals
    # y.in <- apply(Y, 2, function(i) lm.fit(Z, i)$residuals)
  }

  data.list <- list(x=x.in, dp=y.in)

  #variance transform
  if(mode=="MRA"){
    dwt.list<- dwt.vt(data.list, wf, J, "dwt", "zero", "periodic", "auto")
  } else if(mode=="MODWT") {
    dwt.list<- modwt.vt(data.list, wf, J, "periodic", "auto")
  } else {
    dwt.list<- at.vt(data.list, wf, J, "periodic", "auto")
  }

  y.in = dwt.list$dp.n

  pmi <- apply(y.in, 2, function(i) pmi.calc(x.in,i))
  pmi[pmi<0] <- 0
  pic <- sqrt(1-exp(-2*pmi))

  return(list(pic=as.numeric(pic),
              py=dwt.list$dp.n, S=dwt.list$S
              ))
}
#-------------------------------------------------------------------------------
# Calculate Partial Weight
#
# @param x       A vector of response.
# @param z       A matrix containing identified predictors of x.
# @param cpyPIC  Partial informational correlation (cpyPIC).
#
# @return A vector of partial weights(pw) of the same length of z.
#
# @references Sharma, A., Mehrotra, R., 2014. An information theoretic alternative to model a natural system using observational information alone. Water Resources Research, 50(1): 650-660.
pw.calc <- function(x, z, cpyPIC){
  wt <- NA
  Z = as.matrix(z)
  if(ncol(Z)==1) {
    wt <- calc.scaleSTDratio(x, Z)*cpyPIC
  } else {
    for(i in 1:ncol(Z)) wt[i] <- calc.scaleSTDratio(x, Z[,i], Z[,-i])*cpyPIC[i]
  }
  return(list(pw=wt))
}


