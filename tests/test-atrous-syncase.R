#-------------------------------------------------------------------------------------
#a trous algorithm
rm(list=ls()) # remove all variables
graphics.off() # remove all figures

#wavelet transforms and multiresolution analysis
library(R.matlab)
library(matlabr)

library(wavethresh)
library(WASP)
library(ggplot2)
#--------------------------------------------------------------------------------------------
# wavelet family, extension mode and package
wf <- "d4" # wavelet family D8 or db4
boundary <- "periodic"
pad <-  "zero"
if(wf!="haar") v <- as.integer(as.numeric(substr(wf,2,3))/2) else v <- 1

flag.plt = 0
label = c("(a)","(b)")
path <- "~/OneDrive - UNSW/R_Package/WASP/tests/"
#path <- "C:/Users/Ze/OneDrive - UNSW/R_Package/WASP/tests/"
#-------------------------------------------------------------------------------------
###synthetic example
sample=500 #sample size

#frequency, sampled from a given range
fd <- c(3,5,10,15,25,30,55,70,95)

set.seed(101)
data.SW <- data.gen.SW(nobs=sample,fp=95,fd=fd)
#data.SW <- data.gen.SW(nobs=sample,fp=c(15,25,30),fd=fd)

#----------------------------------------------------------------------------
###variable
data <- data.SW
x <- data$x - mean(data$x)
xx <- padding(x,pad="zero")

#Maximum decomposition level J
n <- length(x)
J <- ceiling(log(n/(2*v-1))/log(2)) #(Kaiser, 1994)
if(wf=="haar"&&mode=="MODWT") J = J-1 #since modwt no need a dyadic number size
print(paste0("Calibration: Decomposition Levels J= ",J))


#----------------------------------------------------------------------------
#a trous - WASP
x.at <- at.wd(xx, v=v, nthresh=J, boundary="periodic")

print(sum(abs(x-rowSums(x.at[1:n,]))))                  #additive check
var(x);sum(apply(x.at[1:n,],2,var))                     #variance check
var(xx);sum(apply(x.at,2,var))                          #variance check

if(TRUE){
  limits.x <- c(0,n); limits.y <- c(-2,2)
  mra.plot(x, x.at, limits.x, limits.y)
  p1 <- recordPlot()
}


#----------------------------------------------------------------------------
#a trous - wavethresh
x.at.n <- wavethresh::wd(xx, type="station", filter.number=v, family="DaubExPhase")

x.at.D <- NULL
x.at.D <- accessC(x.at.n,0)
for (i in (J:0)) {

  x.at.D <- cbind(x.at.D, accessD(x.at.n, i))

}

plot.ts(x.at.D)

print(sum(abs(x-rowSums(x.at.D[1:n,]))))                  #additive check
var(x);sum(apply(x.at.D[1:n,],2,var))                     #variance check
var(xx);sum(apply(x.at.D,2,var))                          #variance check

x.at.C <- NULL
for (i in (J+1):0) {

  x.at.C <- cbind(x.at.C, accessC(x.at.n, i))

}

plot.ts(x.at.C)

J=9
x.mra <- matrix(NA, nrow=nrow(x.at), ncol = J+1)

x.mra[,1] <- x.at.C[,J+1]
for (j in 2:(J+1)) {
  x.mra[,j] <- x.at.C[,j-1] - x.at.C[,j]
}

x.mra0<-x.mra
print(sum(abs(x-rowSums(x.mra0[1:n,]))))                  #additive check

var(x);sum(apply(x.mra0[1:n,],2,var))                     #variance check
var(xx);sum(apply(x.mra0,2,var))                          #variance check

plot.ts(x.mra0)

#----------------------------------------------------------------------------
#a trous - Matlab
#set a variable in R and save in a Mat file
#x.mat=rvec_to_matlab(x)
writeMat('./data-raw/x.mat', x=xx)

#run our MATLAB script
run_matlab_script(fname="C:/Users/z5154277/OneDrive - UNSW/R_Package/WASP/tests/Atrous.m")

#read in the output Mat we created in MATLAB
x.atrous = readMat("./data-raw/x-atrous.mat")$dwt

print(sum(abs(x-rowSums(x.atrous[1:n,]))))                  #additive check
var(x);sum(apply(x.atrous[1:n,],2,var))                     #variance check
var(xx);sum(apply(x.atrous,2,var))                          #variance check

plot.ts(x.atrous)

#read in the output Mat we created in MATLAB
x.atrous.mra = readMat("./data-raw/x-atrous-mra.mat")$idwt

print(sum(abs(x-rowSums(x.atrous.mra[1:n,]))))                  #additive check
var(x);sum(apply(x.atrous.mra[1:n,],2,var))                     #variance check
var(xx);sum(apply(x.atrous.mra,2,var))                          #variance check

plot.ts(x.atrous.mra)


#----------------------------------------------------------------------------
# plot.ts(x.atrous.mra); plot.ts(x.at.C[,10:1])
#
# plot.ts(x.at.D); plot.ts(x.atrous)
x.mra <- matrix(NA, nrow=nrow(x.at), ncol = ncol(x.at))

x.mra[,1] <- x.atrous.mra[,1]
for (j in 2:ncol(x.at)) {
  x.mra[,j] <- x.atrous.mra[,j] - x.atrous.mra[,j-1]
}

print(sum(abs(x-rowSums(x.mra[1:n,]))))                  #additive check
var(x);sum(apply(x.mra[1:n,],2,var))                     #variance check
var(xx);sum(apply(x.mra,2,var))                          #variance check

plot.ts(x.mra)
plot.ts(x.at)


#----------------------------------------------------------------------------
#DWT-MRA
x.mra <- waveslim::mra(xx,wf=wf, J=J, method="dwt", boundary="periodic")
x.mra.m <- matrix(unlist(x.mra), ncol=J+1)

print(sum(abs(x-rowSums(x.mra.m[1:n,]))))               #additive check
var(x);sum(apply(x.mra.m[1:n,],2,var))                  #variance check
var(xx);sum(apply(x.mra.m,2,var))                       #variance check

if(TRUE){
  limits.x <- c(0,n); limits.y <- c(-2,2)
  mra.plot(xx, x.mra.m, limits.x, limits.y)
  p2 <- recordPlot()
}

#----------------------------------------------------------------------------
#plot and save
fig <- cowplot::plot_grid(p1,p2, ncol=2, labels = label, label_size = 12, vjust=1, hjust=0)
fig

if(flag.plt){
  ggsave(paste0(path,"Figure_atrous_SW_",wf,".jpg"),
         width = 200, height = 115, units = "mm", dpi = 600, fig)
}
