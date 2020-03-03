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

library(parallel)

flag.save = 0
#setwd("~/OneDrive - UNSW/R_Package/WASP/tests/realcase")
setwd("C:/Users/Ze/OneDrive - UNSW/R_Package/WASP/tests/realcase")

#---------------------------------------
cores=detectCores()
cl <- makeCluster(cores[1]/2)
clusterEvalQ(cl,{library(NPRED)})

for(sc in c(12,36)){
#-------------------------------------------------------------------------------------
###Study Index and Grid
#sc=switch(2, 12,36) #c(6,12,24,36,48) #SPI moving window
k.folds = 2 #cross-validation partition
#Grid = sample(Ind_AWAP.2.5,4)
Grid = c(45,117,142,149) #c(107, 142, 153, 190)

###knn predictive model
k=5 # The number of nearest neighbours used
alpha=0.1 # selection stop criteria
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
#kf =1
#Segement your data by fold using the which() function
testIndexes <- which(folds==kf, arr.ind=TRUE)

#-----------------------------------------------------------------------------------------------------
###calibration and selection
S.k <- vector("list",k.folds-1)
dwt.list.k <- vector("list",k.folds-1)

cal.folds <- seq(1:k.folds)[-kf] # calibration folds
for(i in 1:(k.folds-1)){

  trainIndexes <- which(folds==cal.folds[i], arr.ind=TRUE)
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
  if(mode=="a trous") J=ceiling(log2(n))
  if(wf=="haar"&&mode=="MODWT") J = J-1 #since modwt no need a dyadic number size
  print(paste0("Calibration: Decomposition Levels J= ",J))

  #variance transfrom
  if(mode=="MRA"){
    dwt.list<- mclapply(data.list, function(x) dwt.vt(x, wf, J, method, pad, boundary, cov.opt))
  } else if(mode=="MODWT") {
    dwt.list<- mclapply(data.list, function(x) modwt.vt(x, wf, J, pad, boundary,vt.opt, cov.opt))
  } else {
    dwt.list<- mclapply(data.list, function(x) at.vt(x, wf, pad, boundary, cov.opt))
  }

  S.k[[i]] <- lapply(dwt.list, function(ls) ls$S)
  dwt.list.k[[i]] <- dwt.list
}

#-----------------------------------------
###prepare for validation
if(TRUE){
  S.avg <- vector("list",length=length(Grid)) #gather then compute averaged S
  dwt.list.n <- dwt.list                      #gather all the cal folds together

  for(i in 1:length(Grid)){
    #averaged S
    S.avg[[i]] <- matrix(NA,nrow=J+1,ncol=ncol(data.CI.ts))
    for(nc in 1:ncol(data.CI.ts)){
      tmp <- NULL
      for(j in 1:length(S.k)) tmp <- rbind(tmp, S.k[[j]][[i]][,nc])
      S.avg[[i]][,nc] <-  colMeans(tmp)
    }
    #save S in dwt.list.n
    dwt.list.n[[i]]$S <- S.avg[[i]]

    #create new dwt.list.n, append them together
    dwt.x <- NULL; dwt.dp <- NULL; dwt.dp.n <- NULL
    for(j in 1:length(dwt.list.k)){
      dwt.x <- c(dwt.x, dwt.list.k[[j]][[i]]$x)
      dwt.dp <- rbind(dwt.dp, dwt.list.k[[j]][[i]]$dp)
      dwt.dp.n <- rbind(dwt.dp.n, dwt.list.k[[j]][[i]]$dp.n)
    }

    dwt.list.n[[i]]$x <- dwt.x
    dwt.list.n[[i]]$dp <- dwt.dp
    dwt.list.n[[i]]$dp.n <- dwt.dp.n

  }

}

#-----------------------------------------------------------------------------------------------------
#selection and weights
print(paste0("Selection: Stepwise PIC alpha=",alpha))
sel <- list(); wt <- list()
data <- lapply(dwt.list.n, function(ls) list(x=ls$x, dp=ls$dp))
data.n <- lapply(dwt.list.n, function(ls) list(x=ls$x, dp=ls$dp.n))

clusterExport(cl,varlist=c("alpha"),envir=environment())
tmp1 <- parLapply(cl, data, function(df) stepwise.PIC(df$x, df$dp, alpha, method=F))
tmp2 <- parLapply(cl, data.n, function(df) stepwise.PIC(df$x, df$dp, alpha, method=F))

sel$origin <- lapply(tmp1, function(ls) ls$cpy); wt$origin <- lapply(tmp1, function(ls) ls$wt)
sel$vt <- lapply(tmp2, function(ls) ls$cpy); wt$vt <- lapply(tmp2, function(ls) ls$wt)

sel.cv[[kf]] <- sel
wt.cv[[kf]] <- wt

#-----------------------------------------------------------------------------------------------------
###validation and prediction
# create paired response and predictors dataset for each station
data.list.val <- list()
for(id in seq_along(Grid)){
  x <- data.SPI.ts[testIndexes,Grid[id]]
  dp <- data.CI.ts[testIndexes,]
  data.list.val[[id]] <- list(x=as.numeric(x), dp=dp)
}

#Maximum decomposition level J
n <- nrow(data.CI.ts[testIndexes,])
J <- ceiling(log(n/(2*v-1))/log(2)) #(Kaiser, 1994)
if(mode=="a trous") J=ceiling(log2(n))
if(wf=="haar"&&mode=="MODWT") J = J-1 #since modwt no need a dyadic number size
print(paste0("Validation: Decomposition Levels J= ",J))

#------------------------------------------------------------------------------
#variance transform
if(mode=="MRA"){
    dwt.list.val<- mclapply(1:length(data.list.val), function(i) dwt.vt.val(data.list.val[[i]], J, dwt.list.n[[i]]))
} else if(mode=="MODWT"){
    dwt.list.val<- mclapply(1:length(data.list.val), function(i) modwt.vt.val(data.list.val[[i]], J, dwt.list.n[[i]]))
} else {
    dwt.list.val<- mclapply(1:length(data.list.val), function(i) at.vt.val(data.list.val[[i]], dwt.list.n[[i]]))
}

#------------------------------------------------------------------------------
#validation
RMSE.val <- NULL; df.val <- NULL
for(i in 1:length(dwt.list.val)){
  # cali
  dwt <- dwt.list.n[[i]]
  x.train <- dwt$x
  dp <- dwt$dp
  dp.n <- dwt$dp.n

  # vali
  dwt <- dwt.list.val[[i]]
  x <- dwt$x
  dp.v <- dwt$dp
  dp.n.v <- dwt$dp.n

  s1 <- sel$origin[[i]]; s2 <- sel$vt[[i]]
  w1 <- wt$origin[[i]]; w2 <- wt$vt[[i]]

  if(is.null(s1)) {s1 <- 1;w1 <- 1}
  if(is.null(s2)) {s2 <- 1;w2 <- 1}

  if(TRUE){
    m1 <- knn(x.train, z=as.matrix(dp[,s1]), zout=as.matrix(dp.v[,s1]),k=k,tailcorrection=F,extrap=F)
    m2 <- knn(x.train, z=dp.n[,s2], zout=dp.n.v[,s2], k=k, pw=w2)
    #ts.plot(ts(cbind(x,m1,m2)), col=c("black","red","blue"),lwd=c(2,1,1))
  }
  data.RMSE <- round(sqrt(mean((m1-x)^2)),3)
  dwt.RMSE <- round(sqrt(mean((m2-x)^2)),3)
  RMSE.val <- rbind(RMSE.val,c(data.RMSE,dwt.RMSE))

  df1 <- data.frame(Group=1, No=Grid[i],N=1:n,Pred=m1, Obs=x)
  df2 <- data.frame(Group=2, No=Grid[i],N=1:n,Pred=m2, Obs=x)

  df.val <- rbind(df.val, rbind(df1,df2))

  data.SPI.ref[testIndexes,Grid[i]] <- m1
  data.SPI.mod[testIndexes,Grid[i]] <- m2

}


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
t1 <- table(factor(unlist(sel$origin),levels=1:ncol(data.CI.ts)))
t2 <- table(factor(unlist(sel$vt),levels=1:ncol(data.CI.ts)))

output <- rbind(t1,t2)
rownames(output) <- c("Origin","VT")
print(output)

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


stopCluster(cl)
