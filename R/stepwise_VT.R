#' Calculate stepwise high order VT
#'
#' @param x       A vector of response.
#' @param py      A matrix containing possible predictors of x.
#' @param alpha   The significance level used to judge whether the sample estimate in Equation \deqn{\hat{PIC} = sqrt(1-exp(-2\hat{PI})} is significant or not. A default alpha value is 0.1.
#' @param mode    A mode of variance transfomration, i.e., MRA, MODWT, or AT
#' @param wf      Wavelet family
#'
#' @return A list of 2 elements: the column numbers of the meaningful predictors (cpy), and partial informational correlation (cpyPIC).
#' @export
#'
#' @references Sharma, A., Mehrotra, R., 2014. An information theoretic alternative to model a natural system using observational information alone. Water Resources Research, 50(1): 650-660.
#' @examples
#' data(data1) # AR9 model   x(i)=0.3*x(i-1)-0.6*x(i-4)-0.5*x(i-9)+eps
#' x<-data1[,1]   # response
#' py<-data1[,-1]  # possible predictors
#' stepwise.PIC(x,py)
#'
#' data(data2) # AR4 model:  x(i)=0.6*x(i-1)-0.4*x(i-4)+eps
#' x<-data2[,1]   # response
#' py<-data2[,-1]  # possible predictors
#' stepwise.PIC(x,py)
#'
#' data(data3) # AR1 model  x(i)=0.9*x(i-1)+0.866*eps
#' x<-data3[,1]   # response
#' py<-data3[,-1]  # possible predictors
#' stepwise.PIC(x,py)

stepwise.VT <- function (x, py, alpha = 0.1, mode=c("MRA","MODWT","AT"), wf)
{
  x = as.matrix(x)
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
    return(list(cpy = cpy, cpyPIC = cpyPIC, wt = outwt,lstwet = lstwt,
                x = x, py = py,
                dp=z, dp.n=z.vt, S=S,
                wavelet=wf))
  } else {
    message("None of the provided predictors is related to the response variable")
  }
}



#' Calculate stepwise high order VT
#'
#' @param data    A list of response and identifed predictors
#' @param J       The decomposition level
#' @param dwt     Output from dwt.vt(), including the transformation covariance
#' @param mode    A mode of variance transfomration, i.e., MRA, MODWT, or AT
#'
#' @return        A list of objects, including transformed predictors
#' @export
#'
#' @examples
#'
stepwise.VT.val <- function (data, J, dwt, mode){

  # initialization
  x= data$x; dp= as.matrix(data$dp)
  cpy=dwt$cpy; ncpy=length(cpy)
  method <- "dwt"; boundary <- "periodic"; pad="zero"

  dwt.n = c(dwt, method=method, boundary=boundary, pad=pad)
  dp.n = dp[,cpy]
  if(ncpy>1){
    for(i in 2:ncpy){
      Z=dp[,cpy[1:(i-1)]]
      dp.n[,i] <- knnregl1cv(dp[,cpy[i]], Z)-dp[,cpy[i]]
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

  return(c(dwt.val,py=list(dp)))

}
#' Calculate the ratio of conditional error standard deviations
#'
#' @param x     A vector of response.
#' @param zin   A matrix containing the meaningful predictors selected from a large set of possible predictors (z).
#' @param zout  A matrix containing the remaining possible predictors after taking out the meaningful predictors (zin).
#'
#' @return The STD ratio.
#' @export
#'
#' @references Sharma, A., Mehrotra, R., 2014. An information theoretic alternative to model a natural system using observational information alone. Water Resources Research, 50(1): 650-660.
calc.scaleSTDratio <- function (x, zin, zout)
{
  if(!missing(zout)){
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
#' Calculate PIC
#'
#' @param X       A vector of response.
#' @param Y       A matrix of new predictors.
#' @param Z       A matrix of pre-existing predictors that could be NULL if no prior predictors exist.
#' @param mode    A mode of variance transfomration, i.e., MRA, MODWT, or AT
#' @param wf      Wavelet family
#'
#' @return A list of 2 elements: the partial mutual information (pmi), and partial informational correlation (pic).
#' @export
#'
#' @references Sharma, A., Mehrotra, R., 2014. An information theoretic alternative to model a natural system using observational information alone. Water Resources Research, 50(1): 650-660.
#' @references Galelli S., Humphrey G.B., Maier H.R., Castelletti A., Dandy G.C. and Gibbs M.S. (2014) An evaluation framework for input variable selection algorithms for environmental data-driven models, Environmental Modelling and Software, 62, 33-51, DOI: 10.1016/j.envsoft.2014.08.015.
pic.calc <- function(X, Y, Z, mode, wf) {

  Y=as.matrix(Y)

  if(wf!="haar") v <- as.integer(readr::parse_number(wf)/2) else v <- 1
  #Maximum decomposition level J
  n <- length(X)
  J <- ceiling(log(n/(2*v-1))/log(2)) #(Kaiser, 1994)
  if(mode=="MODWT"&&wf=="haar") J=J-1

  if(is.null(Z)){
    x.in <- X
    y.in <- Y
  } else {
    x.in <- knnregl1cv(X, Z[1:length(X),])-X
    y.in <- apply(Y, 2, function(i) knnregl1cv(i, Z)-i)
  }

  data.list <- list(x=x.in, dp=y.in)

  #variance transform
  if(mode=="MRA"){
    dwt.list<- dwt.vt(data.list, wf, J, "dwt", "zero", "periodic", "auto")
  } else if(mode=="MODWT") {
    dwt.list<- modwt.vt(data.list, wf, J, "periodic", "auto")
  } else {
    dwt.list<- at.vt(data.list, wf, J, "zero", "periodic", "auto")
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
#' Calculate Partial Weight
#'
#' @param x       A vector of response.
#' @param z       A matrix containing identified predictors of x.
#' @param cpyPIC  Partial informational correlation (cpyPIC).
#'
#' @return A vector of partial weights(pw) of the same length of z.
#' @export
#'
#' @references Sharma, A., Mehrotra, R., 2014. An information theoretic alternative to model a natural system using observational information alone. Water Resources Research, 50(1): 650-660.

#' @examples
#'
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


