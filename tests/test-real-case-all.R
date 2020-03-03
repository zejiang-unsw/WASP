#-------------------------------------------------------------------------------------
rm(list=ls()) # remove all variables
graphics.off() # remove all figures

#variance transformation with different wavelet transforms
library(NPRED)

library(zoo)
library(fitdistrplus)
library(waveslim)
library(wavethresh)
# library(WASP)

flag.save = 0
label = c("(a)","(b)")
#setwd("~/OneDrive - UNSW/R_Package/WASP/tests/realcase")
setwd("C:/Users/Ze/OneDrive - UNSW/R_Package/WASP/tests/realcase")
#-------------------------------------------------------------------------------------
###Study Index and Grid
sc=switch(4, 6,12,24,36,48) #SPI moving window
k.folds = 2 #cross-validation partition
Grid = sample(Ind_AWAP.2.5,4)
Grid = c(130, 87, 149, 28)  #177 249 #28 87

###knn predictive model
k=5 # The number of nearest neighbours used
start.yr = c(1910,1); end.yr = c(2009,12)  #study period

### wavelet method selection - for proposed method
mode <- switch(2,"MRA", "MODWT","a trous")
cov.opt <- switch(3,"pos","neg","auto")

if(mode=="MRA") method <- switch(1,"dwt","modwt")
if(mode=="MODWT") vt.opt <- switch(1,"Sxx","Cov")

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

}

df.RMSE <- data.frame(RMSE, dif=RMSE[,1]-RMSE[,2])
par(op); boxplot(df.RMSE$dif)

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

    m1 <- knn(x.train, z=as.matrix(dp[,s1]), zout=as.matrix(dp.v[,s1]),k=k, pw=w1, extrap=F)
    m2 <- knn(x.train, z=dp.n[,s2], zout=dp.n.v[,s2], pw=w2, extrap=T)
    #ts.plot(ts(cbind(x,m1,m2)), col=c("black","red","blue"),lwd=c(2,1,1))

  }
  data.RMSE <- round(sqrt(mean((m1-x)^2)),3)
  dwt.RMSE <- round(sqrt(mean((m2-x)^2)),3)
  RMSE.val <- rbind(RMSE.val,c(data.RMSE,dwt.RMSE))

  data.SPI.ref[-trainIndexes,Grid[i]] <- m1
  data.SPI.mod[-trainIndexes,Grid[i]] <- m2

}

df.RMSE.val <- data.frame(RMSE.val, dif=RMSE.val[,1]-RMSE.val[,2])
par(op);boxplot(df.RMSE.val$dif)

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



