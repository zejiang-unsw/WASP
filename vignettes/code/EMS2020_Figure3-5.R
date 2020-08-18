# devtools::install_github("zejiang-unsw/WASP")
# devtools::install_github("zejiang-unsw/NPRED@v1.0.2")

library(WASP)
library(waveslim)
library(parallel)

#setwd to a folder with /result/

#---------------------------------------
cores=detectCores()
cl <- makeCluster(cores[1]-4)
clusterEvalQ(cl,{library(NPRED)})

#-------------------------------------------------------------------------------------
###load data
data("data.CI")
data("Ind_AWAP.2.5")
data("data.AWAP.2.5")
summary(data.CI); summary(Ind_AWAP.2.5)
summary(data.AWAP.2.5[,Ind_AWAP.2.5])

###Study grid and period
k.folds = 4 #cross-validation partition
Grids = Ind_AWAP.2.5 #run over Australia
#Grids = c(45,117,142,149) # run with sampled grids
start.yr = c(1910,1); end.yr = c(2009,12)  #study period

###knn predictive model
k=0 # The number of nearest neighbours used
alpha=0.1 # selection stop criteria

### wavelet method selection - for proposed method
mode <- switch(2,"MRA", "MODWT","AT")
cov.opt <- switch(1,"auto","pos","neg")
flag = switch(1, "biased", "unbiased")
detrend = F

if(mode=="MRA") method <- switch(1,"dwt","modwt")

# wavelet family, extension mode and package
wf <- "haar" # wavelet family d4 or db2
pad <-  "zero" # padding for non dyadic sample size
boundary <- "periodic"
if(wf!="haar") v <- as.integer(readr::parse_number(wf)/2) else v <- 1 #vanish moments

#------------------------------------------------------
for(sc in c(12,36)){
do.call("data", list(paste0("SPI.",sc)))

###subset data to study period of interest
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
  for(id in seq_along(Grids)){
    x <- data.SPI.ts[trainIndexes,Grids[id]]
    dp <- data.CI.ts[trainIndexes,]
    data.list[[id]] <- list(x=as.numeric(x), dp=dp)
  }

  #Maximum decomposition level J
  n <- nrow(data.CI.ts[trainIndexes,])
  J <- ceiling(log(n/(2*v-1))/log(2)) - 1  #(Kaiser, 1994)
  print(paste0("Calibration: Decomposition Levels J= ",J))

  #variance transfrom
  if(mode=="MRA"){
    dwt.list<- mclapply(data.list, function(x) dwt.vt(x, wf, J, method, pad, boundary, cov.opt, flag, detrend))
  } else if(mode=="MODWT") {
    dwt.list<- mclapply(data.list, function(x) modwt.vt(x, wf, J, boundary, cov.opt, flag, detrend))
  } else {
    dwt.list<- mclapply(data.list, function(x) at.vt(x, wf, J, boundary, cov.opt, flag, detrend))
  }

  S.k[[i]] <- lapply(dwt.list, function(ls) ls$S)
  dwt.list.k[[i]] <- dwt.list
}

#-----------------------------------------
###prepare for validation
if(TRUE){
  S.avg <- vector("list",length=length(Grids)) #gather then compute averaged S
  dwt.list.n <- dwt.list                      #gather all the cal folds together

  for(i in 1:length(Grids)){
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
tmp1 <- parLapply(cl, data, function(df) stepwise.PIC(df$x, df$dp, alpha))
tmp2 <- parLapply(cl, data.n, function(df) stepwise.PIC(df$x, df$dp, alpha))

sel$origin <- lapply(tmp1, function(ls) ls$cpy); wt$origin <- lapply(tmp1, function(ls) ls$wt)
sel$vt <- lapply(tmp2, function(ls) ls$cpy); wt$vt <- lapply(tmp2, function(ls) ls$wt)

sel.cv[[kf]] <- sel
wt.cv[[kf]] <- wt

#-----------------------------------------------------------------------------------------------------
###validation and prediction
# create paired response and predictors dataset for each station
data.list.val <- list()
for(id in seq_along(Grids)){
  x <- data.SPI.ts[testIndexes,Grids[id]]
  dp <- data.CI.ts[testIndexes,]
  data.list.val[[id]] <- list(x=as.numeric(x), dp=dp)
}

#Maximum decomposition level J
n <- nrow(data.CI.ts[testIndexes,])
J <- ceiling(log(n/(2*v-1))/log(2)) - 1 #(Kaiser, 1994)
print(paste0("Validation: Decomposition Levels J= ",J))

#------------------------------------------------------------------------------
#variance transform
if(mode=="MRA"){
    dwt.list.val<- mclapply(1:length(data.list.val), function(i) dwt.vt.val(data.list.val[[i]], J, dwt.list.n[[i]], detrend))
} else if(mode=="MODWT"){
    dwt.list.val<- mclapply(1:length(data.list.val), function(i) modwt.vt.val(data.list.val[[i]], J, dwt.list.n[[i]], detrend))
} else {
    dwt.list.val<- mclapply(1:length(data.list.val), function(i) at.vt.val(data.list.val[[i]], J, dwt.list.n[[i]], detrend))
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
    m1 <- NPRED::knn(x.train, z=as.matrix(dp[,s1]), zout=as.matrix(dp.v[,s1]), k=k, pw=w1, extrap=F)#tailcorrection=F,extrap=F)
    m2 <- NPRED::knn(x.train, z=dp.n[,s2], zout=dp.n.v[,s2], k=k, pw=w2)
    #ts.plot(ts(cbind(x,m1,m2)), col=c("black","red","blue"),lwd=c(2,1,1))
  }
  data.RMSE <- round(sqrt(mean((m1-x)^2)),3)
  dwt.RMSE <- round(sqrt(mean((m2-x)^2)),3)
  RMSE.val <- rbind(RMSE.val,c(data.RMSE,dwt.RMSE))

  df1 <- data.frame(Group=1, No=Grids[i],N=1:n,Pred=m1, Obs=x)
  df2 <- data.frame(Group=2, No=Grids[i],N=1:n,Pred=m2, Obs=x)

  df.val <- rbind(df.val, rbind(df1,df2))

  data.SPI.ref[testIndexes,Grids[i]] <- m1
  data.SPI.mod[testIndexes,Grids[i]] <- m2

}

}

#------------------------------------------------------------------------------
#save file
save(sel.cv, file=paste0("./result/data.SPI.",sc,".selection_",mode,"_",wf,"_",k.folds,"folds.Rdata"))
save(wt.cv, file=paste0("./result/data.SPI.",sc,".weights_",mode,"_",wf,"_",k.folds,"folds.Rdata"))
save(data.SPI.mod, file=paste0("./result/data.SPI.",sc,".mod_",mode,"_",wf,"_",k.folds,"folds.Rdata"))
save(data.SPI.ref, file=paste0("./result/data.SPI.",sc,".ref_",mode,"_",wf,"_",k.folds,"folds.Rdata"))

}


stopCluster(cl)
