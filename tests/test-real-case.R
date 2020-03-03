#-------------------------------------------------------------------------------------
rm(list=ls()) # remove all variables
graphics.off() # remove all figures
op <- par() # store original par setup

#variance transformation with different wavelet transforms
library(zoo)
library(fitdistrplus)
library(waveslim)
library(wavethresh)

library(WASP)
library(NPRED)

flag.save = 1
#setwd("~/OneDrive - UNSW/R_Package/WASP/tests/realcase")
setwd("C:/Users/Ze/OneDrive - UNSW/R_Package/WASP/tests/realcase")

for(sc in c(12,36)){
#-------------------------------------------------------------------------------------
###Study Index and Grid
#sc=switch(2, 12,36) #c(6,12,24,36,48) #SPI moving window
k.folds = 2 #cross-validation partition
#Grid = sample(Ind_AWAP.2.5,4)
Grid = c(130, 87, 149, 28) #177 249

###knn predictive model
k=5 # The number of nearest neighbours used
start.yr = c(1910,1); end.yr = c(2009,12)  #study period

### wavelet method selection - for proposed method
mode <- switch(3,"MRA", "MODWT","a trous")
cov.opt <- switch(3,"pos","neg","auto")

if(mode=="MRA") method <- switch(1,"dwt","modwt")
if(mode=="MODWT") vt.opt <- switch(2,"Sxx","Cov")

# wavelet family, extension mode and package
wf <- "d4" # wavelet family D8 or db4
boundary <- "periodic"
pad <-  "zero"
if(wf!="haar") v <- as.integer(as.numeric(substr(wf,2,3))/2) else v <- 1

#------------------------------------------------------
###load data
data("data.CI")
data("Ind_AWAP.2.5")
data("data.AWAP.2.5")
summary(data.CI); summary(Ind_AWAP.2.5); summary(data.AWAP.2.5)

do.call("data", list(paste0("SPI.",sc)))

###subset data
data.CI.ts <- window(data.CI, start=start.yr, end=end.yr)
start(data.CI.ts); end(data.CI.ts)
data.SPI.ts <- eval(parse(text=paste0("SPI.",sc)))
start(data.SPI.ts); end(data.SPI.ts)

#-------------------------------------------------------------------------------------
###Wavelet-based variance transformation
###output
data.SPI.mod <- matrix(NA,ncol=ncol(data.SPI.ts),nrow = nrow(data.SPI.ts))
data.SPI.ref <- matrix(NA,ncol=ncol(data.SPI.ts),nrow = nrow(data.SPI.ts))
sel.cv <- vector("list",k.folds)
wt.cv <- vector("list",k.folds)

#------------------------------------------------------
###Perform k fold cross validation
#create k.folds folds
folds <- cut(seq(1,nrow(data.SPI.ts)),breaks=k.folds,labels=FALSE)

for(kf in 1:k.folds){
# kf =1
#Segement your data by fold using the which() function
trainIndexes <- which(folds==kf, arr.ind=TRUE)

#create paired response and predictors dataset for each grid
data.list <- list()
for(id in seq_along(Grid)){
  x <- data.SPI.ts[trainIndexes,Grid[id]]
  dp <- data.CI.ts[trainIndexes,]
  data.list[[id]] <- list(x=as.numeric(x), dp=dp)
}

#Maximum decomposition level J
n <- nrow(data.CI.ts[trainIndexes,])
J <- ceiling(log(n/(2*v-1))/log(2)) #(Kaiser, 1994)
if(wf=="haar"&&mode=="MODWT") J = J-1 #since modwt no need a dyadic number size
print(paste0("Calibration: Decomposition Levels J= ",J))

#variance transfrom
if(mode=="MRA"){
  dwt.list<- lapply(data.list, function(x) dwt.vt(x, wf, J, method, pad, boundary, cov.opt))
} else if(mode=="MODWT") {
  dwt.list<- lapply(data.list, function(x) modwt.vt(x, wf, J, pad, boundary,vt.opt, cov.opt))
} else {
  dwt.list<- lapply(data.list, function(x) at.vt(x, wf, pad, boundary, cov.opt))
}

#------------------------------------------------------------------------------
#calibration
sel <- vector("list"); wt <- vector("list")
RMSE <- NULL; df <- NULL
for(i in 1:length(dwt.list)){
  # selection
  dwt <- dwt.list[[i]]
  x <- dwt$x; dp <- dwt$dp; dp.n <- dwt$dp.n
  #plot.ts(cbind(x,dp)); plot.ts(cbind(x,dp.n))

  s1 <- NPRED::stepwise.PIC(x, dp, method=F)
  s2 <- NPRED::stepwise.PIC(x, dp.n, method=F)

  # if(is.null(s1)) {s1$cpy <- 1:ncol(data.CI.ts);s1$wt <- rep(1,ncol(data.CI.ts))}
  # if(is.null(s2)) {s2$cpy <- 1:ncol(data.CI.ts);s2$wt <- rep(1,ncol(data.CI.ts))}

  if(is.null(s1)) {s1$cpy <- 1;s1$wt <- 1}
  if(is.null(s2)) {s2$cpy <- 1;s2$wt <- 1}

  sel[[i]] <- list(origin=s1$cpy, vt=s2$cpy)
  wt[[i]]  <- list(origin=s1$wt, vt=s2$wt)

  if(TRUE){

    m1 <- knn(x, z=as.matrix(dp[,s1$cpy]), zout=as.matrix(dp[,s1$cpy]), k=k, pw=s1$wt, extrap=F)
    m2 <- knn(x, z=dp.n[,s2$cpy], zout=dp.n[,s2$cpy], pw=s2$wt, extrap=T)
    #ts.plot(ts(cbind(x,m1,m2)), col=c("black","red","blue"),lwd=c(2,1,1))
  }

  data.RMSE <- round(sqrt(mean((m1-x)^2)),3)
  dwt.RMSE <- round(sqrt(mean((m2-x)^2)),3)
  RMSE <- rbind(RMSE,c(data.RMSE,dwt.RMSE))

  df1 <- data.frame(Group=1, No=Grid[i],N=1:n,Pred=m1, Obs=x)
  df2 <- data.frame(Group=2, No=Grid[i],N=1:n,Pred=m2, Obs=x)
  df <- rbind(df, rbind(df1,df2))

}
summary(df)

#------------------------------------------------------------------------------
#prediction - cali
###bar plot: improved RMSE
df.RMSE <- data.frame(RMSE, dif=RMSE[,1]-RMSE[,2])
print(df.RMSE)
par(op); barplot(df.RMSE$dif,names.arg = Grid)

###qq plot: sample quantiles
par(mfrow=c(2,length(Grid)/2), mar=c(0.5, 1, 2, 1),# margin of the plot
    oma = c(2, 1, 1, 2), # move plot to the right and up
    mgp=c(2, 1, 0), # move axis labels closer to axis
    bg = "transparent", pty="s", # maximal plotting region
    ps=8)
for(i in Grid){
  tmp <- df[df$No==i,]
  plot(qnorm(ppoints(n)),sort(tmp[tmp$Group==1,"Obs"]),
       xlim=c(-3,3),ylim=c(-3,3), xlab="",ylab="",
       # xlab="Theoretical Quantiles",
       # ylab="Sample Quantiles",
       main=paste0("Sampled Grid: ",i),
       pch=16,
       col="black")
  points(qnorm(ppoints(n)),sort(tmp[tmp$Group==1,"Pred"]),pch=16,col="red")
  points(qnorm(ppoints(n)),sort(tmp[tmp$Group==2,"Pred"]),pch=16,col="blue")

  legend('topleft', box.lty = 0,
         legend=c("Observed","Predicted","Predicted(VT)"),  # text in the legend
         col=c("black","red","blue"),  # point colors
         pch=16)  # specify the point type to be a square

}

###xy plot: observed vs predicted
par(mfrow=c(2,length(Grid)/2), mar=c(0.5, 1, 2, 1),# margin of the plot
    oma = c(2, 1, 1, 2), # move plot to the right and up
    mgp=c(2, 1, 0), # move axis labels closer to axis
    bg = "transparent", pty="s", # maximal plotting region
    ps=8)
for(i in Grid){
  tmp <- df[df$No==i,]

  plot(tmp[tmp$Group==1,"Obs"],tmp[tmp$Group==1,"Pred"],
       xlim=c(-3,3),ylim=c(-3,3), xlab="",ylab="",
       # xlab="Theoretical Quantiles",
       # ylab="Sample Quantiles",
       main=paste0("Station No. ",i),
       pch=16,
       col="red")
  points(tmp[tmp$Group==2,"Obs"],tmp[tmp$Group==2,"Pred"],pch=16,col="blue")
  abline(a=0,b=1)

  legend('topleft', box.lty = 0,
         legend=c("Predicted","Predicted(VT)"),  # text in the legend
         col=c("black","red","blue"),  # point colors
         pch=16)  # specify the point type to be a square

}

#------------------------------------------------------------------------------
# create paired response and predictors dataset for each station
data.list.val <- list()
for(id in seq_along(Grid)){
  x <- data.SPI.ts[-trainIndexes,Grid[id]]
  dp <- data.CI.ts[-trainIndexes,]
  data.list.val[[id]] <- list(x=as.numeric(x), dp=dp)
}

#Maximum decomposition level J
n <- nrow(data.CI.ts[-trainIndexes,])
J <- ceiling(log(n/(2*v-1))/log(2)) #(Kaiser, 1994)
if(wf=="haar"&&mode=="MODWT") J = J-1 #since modwt no need a dyadic number size
print(paste0("Validation: Decomposition Levels J= ",J))

#------------------------------------------------------------------------------
#variance transform
if(mode=="MRA"){
    dwt.list.val<- lapply(1:length(data.list.val), function(i) dwt.vt.val(data.list.val[[i]], J, dwt.list[[i]]))
} else if(mode=="MODWT"){
    dwt.list.val<- lapply(1:length(data.list.val), function(i) modwt.vt.val(data.list.val[[i]], J, dwt.list[[i]]))
} else {
    dwt.list.val<- lapply(1:length(data.list.val), function(i) at.vt.val(data.list.val[[i]], dwt.list[[i]]))
}

#------------------------------------------------------------------------------
#validation
RMSE.val <- NULL; df.val <- NULL
for(i in 1:length(dwt.list.val)){
  # cali
  dwt <- dwt.list[[i]]
  x.train <- dwt$x
  dp <- dwt$dp
  dp.n <- dwt$dp.n

  # vali
  dwt <- dwt.list.val[[i]]
  x <- dwt$x
  dp.v <- dwt$dp
  dp.n.v <- dwt$dp.n

  s1 <- sel[[i]]$origin; s2 <- sel[[i]]$vt
  w1 <- wt[[i]]$origin; w2 <- wt[[i]]$vt

  if(TRUE){

    m1 <- knn(x.train, z=as.matrix(dp[,s1]), zout=as.matrix(dp.v[,s1]),k=k, tailcorrection = F, pw=w1, extrap=F)
    m2 <- knn(x.train, z=dp.n[,s2], zout=dp.n.v[,s2], k=k, tailcorrection=F, pw=w2, extrap=T)
    #ts.plot(ts(cbind(x,m1,m2)), col=c("black","red","blue"),lwd=c(2,1,1))

  }
  data.RMSE <- round(sqrt(mean((m1-x)^2)),3)
  dwt.RMSE <- round(sqrt(mean((m2-x)^2)),3)
  RMSE.val <- rbind(RMSE.val,c(data.RMSE,dwt.RMSE))

  df1 <- data.frame(Group=1, No=Grid[i],N=1:n,Pred=m1, Obs=x)
  df2 <- data.frame(Group=2, No=Grid[i],N=1:n,Pred=m2, Obs=x)

  df.val <- rbind(df.val, rbind(df1,df2))

  data.SPI.ref[-trainIndexes,Grid[i]] <- m1
  data.SPI.mod[-trainIndexes,Grid[i]] <- m2

}
summary(df.val)

#------------------------------------------------------------------------------
#prediction - vali
###bar plot: improved RMSE
df.RMSE <- data.frame(RMSE.val, dif=RMSE.val[,1]-RMSE.val[,2])
print(df.RMSE)
par(op);barplot(df.RMSE$dif,names.arg = Grid)

###qq plot: sample quantiles
par(mfrow=c(2,length(Grid)/2), mar=c(0.5, 1, 2, 1),# margin of the plot
    oma = c(2, 1, 1, 2), # move plot to the right and up
    mgp=c(2, 1, 0), # move axis labels closer to axis
    bg = "transparent", pty="s", # maximal plotting region
    ps=8)
for(i in Grid){
  tmp <- df.val[df.val$No==i,]
  plot(qnorm(ppoints(n)),sort(tmp[tmp$Group==1,"Obs"]),
       xlim=c(-3,3),ylim=c(-3,3), xlab="",ylab="",
       # xlab="Theoretical Quantiles",
       # ylab="Sample Quantiles",
       main=paste0("Sampled Grid: ",i),
       pch=16,
       col="black")
  points(qnorm(ppoints(n)),sort(tmp[tmp$Group==1,"Pred"]),pch=16,col="red")
  points(qnorm(ppoints(n)),sort(tmp[tmp$Group==2,"Pred"]),pch=16,col="blue")

  legend('topleft', box.lty = 0,
         legend=c("Observed","Predicted","Predicted(VT)"),  # text in the legend
         col=c("black","red","blue"),  # point colors
         pch=16)  # specify the point type to be a square

}

###xy plot: observed vs predicted
par(mfrow=c(2,length(Grid)/2), mar=c(0.5, 1, 2, 1),# margin of the plot
    oma = c(2, 1, 1, 2), # move plot to the right and up
    mgp=c(2, 1, 0), # move axis labels closer to axis
    bg = "transparent", pty="s", # maximal plotting region
    ps=8)
for(i in Grid){
  tmp <- df.val[df.val$No==i,]
  plot(tmp[tmp$Group==1,"Obs"],tmp[tmp$Group==1,"Pred"],
       xlim=c(-3,3),ylim=c(-3,3), xlab="",ylab="",
       # xlab="Theoretical Quantiles",
       # ylab="Sample Quantiles",
       main=paste0("Station No. ",i),
       pch=16,
       col="red")
  points(tmp[tmp$Group==2,"Obs"],tmp[tmp$Group==2,"Pred"],pch=16,col="blue")
  abline(a=0,b=1)

  legend('topleft', box.lty = 0,
         legend=c("Predicted","Predicted(VT)"),  # text in the legend
         col=c("black","red","blue"),  # point colors
         pch=16)  # specify the point type to be a square

}



#------------------------------------------------------------------------------
#selection
t1 <- table(factor(unlist(sapply(sel, function(x) x[[1]])),levels=1:ncol(data.CI.ts)))
t2 <- table(factor(unlist(sapply(sel, function(x) x[[2]])),levels=1:ncol(data.CI.ts)))

output <- rbind(t1,t2)
rownames(output) <- c("Origin","VT")
print(output)

sel.cv[[kf]] <- sel
wt.cv[[kf]] <- wt

}

#------------------------------------------------------------------------------
#save file
summary(sel.cv)
summary(wt.cv)
# summary(data.SPI.mod)
# summary(data.SPI.ref)
if(flag.save){
if(mode!="MODWT") {
  save(sel.cv, file=paste0("data.SPI.",sc,".selection_",mode,"_",wf,"_",cov.opt,"_",k.folds,"folds.noEns.Rdata"))
  save(wt.cv, file=paste0("data.SPI.",sc,".weights_",mode,"_",wf,"_",cov.opt,"_",k.folds,"folds.noEns.Rdata"))
  save(data.SPI.mod, file=paste0("data.SPI.",sc,".mod.vt_",mode,"_",wf,"_",cov.opt,"_",k.folds,"folds.noEns.Rdata"))
  save(data.SPI.ref, file=paste0("data.SPI.",sc,".ref_",mode,"_",wf,"_",cov.opt,"_",k.folds,"folds.noEns.Rdata"))
} else {
  save(sel.cv, file=paste0("data.SPI.",sc,".selection_",mode,"_",vt.opt,"_",wf,"_",cov.opt,"_",k.folds,"folds.noEns.Rdata"))
  save(wt.cv, file=paste0("data.SPI.",sc,".weights_",mode,"_",vt.opt,"_",wf,"_",cov.opt,"_",k.folds,"folds.noEns.Rdata"))
  save(data.SPI.mod, file=paste0("data.SPI.",sc,".mod.vt_",mode,"_",vt.opt,"_",wf,"_",cov.opt,"_",k.folds,"folds.noEns.Rdata"))
  save(data.SPI.ref, file=paste0("data.SPI.",sc,".ref_",mode,"_",vt.opt,"_",wf,"_",cov.opt,"_",k.folds,"folds.noEns.Rdata"))
}
}


}
