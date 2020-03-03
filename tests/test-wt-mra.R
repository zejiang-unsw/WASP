#-------------------------------------------------------------------------------------
rm(list=ls()) # remove all variables
graphics.off() # remove all figures

#wavelet transforms and multiresolution analysis
library(waveslim)
library(wavethresh)
library(WASP)
library(ggplot2)
#--------------------------------------------------------------------------------------------
# wavelet family, extension mode and package
wf <- "d4" # wavelet family D8 or db4
boundary <- "periodic"
pad <-  "zero"
if(wf!="haar") v <- as.integer(as.numeric(substr(wf,2,3))/2) else v <- 1

flag.plt = 1
label = c("(a)","(b)","(c)")
path <- "~/OneDrive - UNSW/R_Package/WASP/tests/"
#path <- "C:/Users/Ze/OneDrive - UNSW/R_Package/WASP/tests/"
#-------------------------------------------------------------------------------------
###synthetic example
sample=500 #sample size

#frequency, sampled from a given range
fd <- c(3,5,10,15,25,30,55,70,95)

set.seed(101)
#data.SW <- data.gen.SW(nobs=sample,fp=95,fd=fd)
#data.SW <- data.gen.SW(nobs=sample,fp=c(15,25,30),fd=fd)
data.SW <- data.gen.Rossler(time=seq(0,50,length.out = sample))
#----------------------------------------------------------------------------
###variable
# data <- data.SW
# x <- data$x - mean(data$x)
x= scale(data.SW$x,scale=F)
xx <- padding(x,pad="zero")

#Maximum decomposition level J
n <- length(x)
J <- ceiling(log(n/(2*v-1))/log(2)) #(Kaiser, 1994)
#if(wf=="haar"&&mode=="MODWT") J = J-1 #since modwt no need a dyadic number size
print(paste0("Calibration: Decomposition Levels J= ",J))

range(x)
limits.y <- c(-25,25)
#----------------------------------------------------------------------------
###decomposition
#DWT
x.dwt <- waveslim::dwt(xx, wf=wf, n.levels = J, boundary="periodic")
x.dwt.m <- unlist(x.dwt)

var(xx); sum(unlist(lapply(x.dwt,var)))                  #variance check
#Variance decomposition
1/(length(xx)-1)*sum(unlist(lapply(x.dwt,function(z) norm(as.matrix(z),"2")))^2) - mean(xx)^2

#----------------------------------------------------------------------------
#DWT-MRA
x.mra <- waveslim::mra(xx,wf=wf, J=J, method="dwt", boundary="periodic")
x.mra.m <- matrix(unlist(x.mra), ncol=J+1)

print(sum(abs(x-rowSums(x.mra.m[1:n,]))))               #additive check
var(x);sum(apply(x.mra.m[1:n,],2,var))                  #variance check
var(xx);sum(apply(x.mra.m,2,var))                       #variance check

if(TRUE){
  limits.x <- c(0,n); #limits.y <- c(-3,3)
  mra.plot(x, x.mra.m, limits.x, limits.y)
  p1 <- recordPlot()
}
#----------------------------------------------------------------------------
#MODWT
x.modwt <- waveslim::modwt(x, wf=wf, n.levels = J, boundary = "periodic")
x.modwt.m <- matrix(unlist(x.modwt), ncol=J+1)

print(sum(abs(x-rowSums(x.modwt.m))))                   #additive check
var(x);sum(apply(x.modwt.m,2,var))                      #variance check

if(TRUE){
  limits.x <- c(0,n); #limits.y <- c(-3,3)
  par(oma = c(2, 1, 1, 1))
  mra.plot(x, x.modwt.m, limits.x, limits.y)
  p2 <- recordPlot()
}


#----------------------------------------------------------------------------
#a trous
x.at <- at.wd(xx, v=v, nthresh=J, boundary="periodic")

print(sum(abs(x-rowSums(x.at[1:n,]))))                  #additive check
var(x);sum(apply(x.at[1:n,],2,var))                     #variance check
var(xx);sum(apply(x.at,2,var))                          #variance check

if(TRUE){
  limits.x <- c(0,n); #limits.y <- c(-5,5)
  mra.plot(x, x.at, limits.x, limits.y)
  p3 <- recordPlot()
}
#----------------------------------------------------------------------------
#plot and save
fig <- cowplot::plot_grid(p1,p2,p3, ncol=3, labels = label, label_size = 12, vjust=1, hjust=0)
fig

if(flag.plt){
  ggsave(paste0(path,"Figure_WT_Rossler_",wf,".jpg"),
         width = 200, height = 115, units = "mm", dpi = 600, fig)
}
