######################################################################
rm(list=ls()) # remove all variables
graphics.off() # remove all figures

setwd("C:/Users/z5154277/OneDrive - UNSW/PhD Research Work/4th Paper/JoH_2020")
source("radarchart.R")

library(WASP)
library(synthesis)

#library(plot3D)
#library(rgl) # the best tool to work in 3D from R
library(ggplot2)
library(cowplot)

library(fmsb) # radar chart in R
library(scales)

#--------------------------------------------------
flag.plt = 0
flag.save= 0

sample= 100000
sample.cal=sample/2

k=5
s=c(0.1,0.5,1.0)
ar=0.6
set.seed(2020)

flag.noise = switch(1,"white","red")
model = switch(1,"Rossler","Lorenz")

### wavelet method selection - for proposed method
mode <- switch(1,"MRA", "MODWT","a trous")
cov.opt <- switch(1,"auto","pos","neg")

if(mode=="MRA") method <- switch(1,"dwt","modwt")
if(mode=="MODWT") vt.opt <- switch(2,"Sxx","Cov")

# wavelet family, extension mode and package
wf <- "haar" # wavelet family D8 or db4
method <- switch(1,"dwt","modwt")
boundary <- "periodic"
pad <-  "zero"

#Maximum decomposition level J
# v is number of vanishing moment of dbv
if(wf!="haar") v <- as.integer(as.numeric(substr(wf,2,3))/2) else v <- 1

#########################################################################################################
###synthetic example - Rossler
ts.list <- list()
for(i in seq_along(s)){
  if(model=="Rossler")
    ts.r <- data.gen.Rossler(a = 0.2, b = 0.2, w = 5.7, start = c(-2, -10, 0.2),
                             time = seq(0, 50, length.out = sample))
  else if(model=="Lorenz")
    ts.r <- data.gen.Lorenz(sigma = 10, beta = 8/3,rho = 28,start = c(-13, -14, 47),
                            time = seq(0, 50, length.out = sample))

  #add noise
  if(flag.noise=="white"){

    noise=rnorm(n = sample,mean=0, sd=s[i])
    ts.r$x <- ts(ts.r$x + rnorm(n = sample, mean=0, sd=s[i]))
    ts.r$y <- ts(ts.r$y + rnorm(n = sample, mean=0, sd=s[i]))
    ts.r$z <- ts(ts.r$z + rnorm(n = sample, mean=0, sd=s[i]))
  } else {

    noise=arima.sim(n = sample, list(order = c(1,0,0), ar = ar),sd=s[i])
    ts.r$x <- ts(ts.r$x + arima.sim(n = sample, list(order = c(1,0,0), ar = ar),sd=s[i]))
    ts.r$y <- ts(ts.r$y + arima.sim(n = sample, list(order = c(1,0,0), ar = ar),sd=s[i]))
    ts.r$z <- ts(ts.r$z + arima.sim(n = sample, list(order = c(1,0,0), ar = ar),sd=s[i]))
  }

  ts.list[[i]]<- ts.r
}

#--------------------------------------------------
#calibration dataset
data.list <- lapply(ts.list, function(ts) list(x=ts$z[1:sample.cal],
                                               dp=cbind(ts$x[1:sample.cal],ts$y[1:sample.cal])))

n <- sample.cal
J <- ceiling(log(n/(2*v-1))/log(2))
if(wf=="haar"&&mode=="MODWT") J = J-1 #since modwt no need a dyadic number size
print(paste0("Calibration: Decomposition Levels J= ",J))

#variance transfrom
if(mode=="MRA"){
  dwt.list<- lapply(data.list, function(x) dwt.vt(x, wf, J, method, pad, boundary, cov.opt))
} else if(mode=="MODWT") {
  dwt.list<- lapply(data.list, function(x) modwt.vt(x, wf, J, pad, boundary,vt.opt, cov.opt))
} else {
  dwt.list<- lapply(data.list, function(x) at.vt(x, wf, J, pad, boundary, cov.opt))
}

#--------------------------------------------------
# calibration
df <- NULL;data.RMSE<-NULL;dwt.RMSE<-NULL
sd.cal<-NULL; cor.cal<-NULL
for(i in 1:length(dwt.list)){

  dwt <- dwt.list[[i]]
  dp <- dwt$dp; dp.n <- dwt$dp.n; x <- dwt$x

  m1 <- FNN::knn.reg(dp, y=x, k=k)$pred
  m2 <- FNN::knn.reg(dp.n, y=x, k=k)$pred

  data.RMSE <-c(data.RMSE, round(sqrt(mean((x-m1)^2)),3))
  dwt.RMSE <- c(dwt.RMSE, round(sqrt(mean((x-m2)^2)),3))

  sd.cal <- cbind(sd.cal, as.vector(c(sd(x),sd(m1),sd(m2))))
  cor.cal <-cbind(cor.cal, cor(cbind(x,m1,m2))[,1])

  df1 <- data.frame(Group=1, s=s[i], No=1:sample.cal,Pred=m1, Obs=x)
  df2 <- data.frame(Group=2, s=s[i], No=1:sample.cal,Pred=m2, Obs=x)

  df <- rbind(df, rbind(df1,df2))

}
#summary(df)
print(rbind(data.RMSE,dwt.RMSE))

t1 <- rbind(0,data.RMSE,dwt.RMSE)
sd.cal;cor.cal
#########################################################################################################
#validataion dataset
data.list.val <- lapply(ts.list, function(ts) list(x=ts$z[(sample.cal+1):sample],
                                               dp=cbind(ts$x[(sample.cal+1):sample],
                                                        ts$y[(sample.cal+1):sample])))
sample.val <- sample-sample.cal
n <- sample.val
J <- ceiling(log(n/(2*v-1))/log(2))
if(wf=="haar"&&mode=="MODWT") J = J-1 #since modwt no need a dyadic number size
print(paste0("Validation: Decomposition Levels J= ",J))

#--------------------------------------------------
#variance transform
if(mode=="MRA"){
  dwt.list.val<- lapply(1:length(data.list.val), function(i) dwt.vt.val(data.list.val[[i]], J, dwt.list[[i]]))
} else if(mode=="MODWT"){
  dwt.list.val<- lapply(1:length(data.list.val), function(i) modwt.vt.val(data.list.val[[i]], J, dwt.list[[i]]))
} else {
  dwt.list.val<- lapply(1:length(data.list.val), function(i) at.vt.val(data.list.val[[i]], J, dwt.list[[i]]))
}

#--------------------------------------------------
# validation
df.val <- NULL;data.RMSE <-NULL;dwt.RMSE<-NULL
sd.val<-NULL; cor.val<-NULL
for(i in 1:length(dwt.list.val)){

  dwt <- dwt.list[[i]]
  dp <- dwt$dp; dp.n <- dwt$dp.n; x.train <- dwt$x

  dwt <- dwt.list.val[[i]]
  dp.v <- dwt$dp; dp.n.v <- dwt$dp.n; x <- dwt$x

  m1 <- FNN::knn.reg(train=dp, test=dp.v, y=x.train, k=k)$pred
  m2 <- FNN::knn.reg(train=dp.n, test=dp.n.v, y=x.train, k=k)$pred

  data.RMSE <-c(data.RMSE, round(sqrt(mean((m1-x)^2)),3))
  dwt.RMSE <- c(dwt.RMSE, round(sqrt(mean((m2-x)^2)),3))

  sd.val <- cbind(sd.val, as.vector(c(sd(x),sd(m1),sd(m2))))
  cor.val <- cbind(cor.val, cor(cbind(x,m1,m2))[,1])

  df1 <- data.frame(Group=1, s=s[i], No=1:sample.val,Pred=m1, Obs=x)
  df2 <- data.frame(Group=2, s=s[i], No=1:sample.val,Pred=m2, Obs=x)

  df.val <- rbind(df.val, rbind(df1,df2))

}

#summary(df.val)
print(rbind(data.RMSE,dwt.RMSE))

t2 <- rbind(0,data.RMSE,dwt.RMSE)
sd.val;cor.val
#--------------------------------------------------
print(cbind(t1,t2)) #slightly different from the one presented in the WRR paper
print(cbind(sd.cal, sd.val))
print(cbind(cor.cal,cor.val))

max.rmse = ceiling(max(unlist(cbind(t1,t2)))); min.rmse=floor(min(unlist(cbind(t1,t2))))
max.sd = ceiling(max(unlist(cbind(sd.cal, sd.val)))); min.sd=floor(min(unlist(cbind(sd.cal, sd.val))))
max.r = ceiling(max(unlist(cbind(cor.cal,cor.val)))); min.r=floor(min(unlist(cbind(cor.cal,cor.val))))


#--------------------------------------------------
if(T){
  jpeg(paste0("Figure6_JoH2020_",model,".jpg"),
       width = 230, height = 115, units = "mm", res = 600)

if(TRUE){
  layout(matrix(1:6,nrow=2))
  par(mar=c(0, 0, 0, 0),# margin of the plot
      oma = c(0, 0, 0, 0), # move plot to the right and up
      mgp=c(0, 0, 0), # move axis labels closer to axis
      bg = "transparent", pty="m", # maximal plotting region
      cex.lab=2, adj=0.15,
      ps=12)

  # Color vector
  color=c("red", "blue")
  colors_border=c("black",alpha(color,1))
  colors_in=c("transparent",alpha(color,0.3))

for(i in seq_along(s)){
  # reshape data
  #data <- as.data.frame(cbind(t1[,i], sd.cal[2:3,i], cor.cal[,i]))
  data <- as.data.frame(cbind(t1[,i], sd.cal[,i], cor.cal[,i]))

  colnames(data) <- c("RMSE", "Standard Deviation", "Correlation")
  rownames(data) <- c("obs", paste0("predict" , c("","(vt)")))

  # the max and min of each variable to show on the plot!
  data <- rbind(c(max.rmse,max.sd,max.r), rep(0,3), data)

  data
  # plot with default options:
  radarchart1(data, axistype=2, #title=paste0("s= ", s[i]),
             #custom polygon
             pcol=colors_border , pfcol=colors_in , plwd=1, plty=1,
             pty=0:2,
             #custom the grid
             #cglcol="grey", cglty=1, #axislabcol="grey", #caxislabels=seq(0,20,5), cglwd=0.8,
             #custom labels
             vlcex=0.8)
  title(paste0("s= ", as.character(format(s[i],nsmall=1))),line=-2)
  if(i==3)
  # Add a legend
  legend(x=0.6, y=1.2, title= "Calibration", legend = rownames(data[-c(1,2),]), bty = "n",
         pch=0:2, col= colors_border, cex=1.2)

  #--------------------------------------------------
  #data <- as.data.frame(cbind(t2[,i], sd.val[2:3,i], cor.val[,i]))
  data <- as.data.frame(cbind(t2[,i], sd.val[,i], cor.val[,i]))

  colnames(data) <- c("RMSE", "Standard Deviation", "Correlation")
  rownames(data) <- c("obs", paste0("predict" , c("","(vt)")))

  # the max and min of each variable to show on the plot!
  data <- rbind(c(max.rmse,max.sd,max.r), c(min.rmse,min.sd,min.r), data)

  # plot with default options:
  radarchart1(data, axistype=2, #title=paste0("s= ", s[i]),
             #custom polygon
             pcol=colors_border, pfcol=colors_in, plwd=1, plty=1,
             pty=0:2,
             #custom the grid
             #cglcol="grey", cglty=1, #axislabcol="grey", #caxislabels=seq(0,20,5), cglwd=0.8,
             #custom labels
             vlcex=0.8)
  title(paste0("s= ", as.character(format(s[i],nsmall=1))),line=-2)
  if(i==3)
  # Add a legend
  legend(x=0.6, y=1.2, title= "Validation", legend = rownames(data[-c(1,2),]), bty = "n",
         pch=0:2, col= colors_border, cex=1.2)

}
}
  dev.off()

}
