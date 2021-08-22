library(WASP)
library(NPRED)

#-------------------------------------------------------------------------------------
#load response and predictor variables
data(SPI.12); data(data.CI); data(Ind_AWAP.2.5)
#study grids and period
Grid <- sample(Ind_AWAP.2.5,1)
Grid = 149 #Grid used in the SI
SPI.12.ts <- window(SPI.12, start=c(1910,1),end=c(2009,12))
data.CI.ts <- window(data.CI, start=c(1910,1),end=c(2009,12))
#partition into two folds
folds <- cut(seq(1,nrow(SPI.12.ts)),breaks=2,labels=FALSE)
sub.cali <- which(folds==1, arr.ind=TRUE); sub.vali <- which(folds==2, arr.ind=TRUE)
#-------------------------------------------------------------------------------------
###calibration and selection
data <- list(x=SPI.12.ts[sub.cali,Grid],dp=data.CI.ts[sub.cali,])

#variance transformation - calibration
modwt <- modwt.vt(data, wf="haar", J=8, boundary="periodic")

#stepwise PIC selection
sel <- NPRED::stepwise.PIC(modwt$x, modwt$dp.n)
#-------------------------------------------------------------------------------------
###validation and prediction
data.val <- list(x=SPI.12.ts[sub.vali,Grid],dp=data.CI.ts[sub.vali,])

#variance transformation - validation
modwt.val <- modwt.vt.val(data.val, J=8, modwt)

#knn prediction
cpy <- sel$cpy; wt <- sel$wt
x=data$x; z=modwt$dp.n[,cpy]; zout=modwt.val$dp.n[,cpy]
mod <- knn(x, z, zout, k=5, pw=wt, extrap=T)

#-------------------------------------------------------------------------------------
###plot
start.cal <- c(1910,1); start.val <- c(1960,1)
ndim = ncol(data.CI.ts); CI.names = colnames(data.CI.ts)
op <- par(mfcol=c(ndim+1,2),mar=c(2,4,2,2),mgp=c(1.5,0.5,0),
    bg = "white",pty="m", cex.lab=1.5, ps=8)
#----------------------------------------------
#plot before and after vt - calibration
if(TRUE){

  x <- ts(modwt$x, start=start.cal, frequency = 12)
  dp <- ts(modwt$dp, start=start.cal, frequency = 12)
  dp.n <- ts(modwt$dp.n, start=start.cal, frequency = 12)

  ts.plot(x,xlab=NA, main=paste0("Sampled Grid: ", Grid), ylab=paste0("SPI",12), col=c("black"),lwd=c(1))
  for(nc in 1:ndim)
    ts.plot(dp[,nc],dp.n[,nc],xlab=NA,ylab=paste0(CI.names[nc]),
            col=c("red","blue"),lwd=c(1,1),lty=c(1,2))

}
#----------------------------------------------
#plot before and after vt - validation
if(TRUE){

  x <- ts(modwt.val$x, start=start.val, frequency = 12)
  dp <- ts(modwt.val$dp, start=start.val, frequency = 12)
  dp.n <- ts(modwt.val$dp.n, start=start.val, frequency = 12)

  ts.plot(x, xlab=NA, main=paste0("Sampled Grid: ",Grid), ylab=paste0("SPI",12), col=c("black"),lwd=c(1))
  for(nc in 1:ndim)
    ts.plot(dp[,nc],dp.n[,nc],xlab=NA,ylab=paste0(CI.names[nc]),
            col=c("red","blue"),lwd=c(1,1),lty=c(1,2))

}

par(op)
