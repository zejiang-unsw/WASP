#--------------------------------------------------------------------------------------------
# rm(list=ls()) # remove all variables
# graphics.off() # remove all figures

setwd("C:/Users/Ze/OneDrive - UNSW/R_Package/WASP/tests")

#variance transformation with different wavelet transforms
library(zoo)
library(fitdistrplus)
library(waveslim)
library(wavethresh)
library(WASP)

library(rgdal)
library(ggplot2)
#--------------------------------------------------------------------------------------------
### wavelet method selection
mode <- switch(2,"MRA", "MODWT","a trous")
cov.opt <- switch(1,"auto", "pos","neg")
if(mode=="MRA") method <- switch(1,"dwt","modwt")
if(mode=="MODWT") vt.opt <- switch(2,"Sxx","Cov")

#study Index and Period
sc=c(12,36)
Grid = 117 #177 249
end.calc = c(2009,12)
start.yr = c(1910,1); end.yr = c(2009,12)

#plot save
flag.plt = 1
label = c("(a)","(b)")

#--------------------------------------------------------------------------------------------
###Real-world example
#import basemap
map <- "AUS_boundary_simplified.shp"
Aus_map <- readOGR(map)

#load response and predictors
data("data.CI")
data("data.AWAP.2.5")

#sampled grid in Southern Australia
data("Ind_AWAP.2.5")
data("lat_lon.2.5")
Grid.xy=lat_lon.2.5[Grid,]
Grid.xy=lat_lon.2.5[c(45,117,142,149),]

#overview
par(mfrow=c(1,1),pty="m")
plot(Aus_map)
points(Grid.xy, col="red", pch=16)

#------------------------------------------------------
p.list <- list()
for(i in sc){
start.calc = switch(as.character(i),"12"=c(1909,1),"36"=c(1907,1))

#SPI calculation
AWAP.mon.ts <- window(ts(data.AWAP.2.5, start=c(1900,1), frequency = 12),
                      start=start.calc, end=end.calc)

SPI  <- matrix(NA, nrow=nrow(AWAP.mon.ts), ncol=ncol(AWAP.mon.ts))
SPI[,Ind_AWAP.2.5] <- SPI.calc(AWAP.mon.ts[,Ind_AWAP.2.5], sc=i, method="mme")

SPI.ts <- ts(SPI, start=start.calc, frequency = 12)
start(SPI.ts); end(SPI.ts)

# #-------------------
# SPI1[,Ind_AWAP.2.5] <- SPI.calc(AWAP.mon.ts[,Ind_AWAP.2.5], sc=i, method="mle")
# SPI1 <- matrix(NA, nrow=nrow(AWAP.mon.ts), ncol=ncol(AWAP.mon.ts))
#
# SPI.ts1 <- ts(SPI1, start=start.calc, frequency = 12)
# start(SPI.ts1); end(SPI.ts1)
#
# SPI.12 <- get(load("~/OneDrive - UNSW/PhD Research Work/EMS Software/EMS_2019/data-raw/data.SPI.12.obs.2.5.Rdata"))
# SPI.36 <- get(load("~/OneDrive - UNSW/PhD Research Work/EMS Software/EMS_2019/data-raw/data.SPI.36.obs.2.5.Rdata"))
#
#
# if(TRUE){
#   #Grid = sample(Ind_AWAP.2.5,1)
#
#   par(mfcol=c(1,1), mar=c(0,3,2,1),# margin of the plot
#       oma = c(2, 1, 1, 2), # move plot to the right and up
#       mgp=c(1, 0.5, 0), # move axis labels closer to axis
#       bg = "transparent", pty="m", # maximal plotting region
#       ps=8)
#   SPI.12.ts <- ts(SPI.12, start=c(1910,1), frequency = 12)
#   SPI.36.ts <- ts(SPI.36, start=c(1910,1), frequency = 12)
#
#   SPI.ref.ts <- switch(as.character(i),"12"=SPI.12.ts,"36"=SPI.36.ts)
#   SPI.ts.sub <- window(SPI.ts,start=start.yr,end=end.yr)
#
#   ts.plot(SPI.ts.sub[,Grid], SPI.ref.ts[,Grid],
#           main=paste0("Sample Grid: ",Grid),
#           col=c("red","blue"),lwd=c(2,1))
#
# }

#------------------------------------------------------
#Wavelet-based variance transformation
#data set preparation for vt
x <- window(SPI.ts[,Grid],start=start.yr,end=end.yr)
dp <- window(data.CI[,3:ncol(data.CI)],start=start.yr,end=end.yr)
CI.names <- colnames(dp)

data.list <- lapply(seq_along(Grid), function(id) list(x=as.matrix(x)[,id],dp=dp))

###wavelet transform
# wavelet family, extension mode and package
wf <- "d4" # wavelet family D8 or db4
boundary <- "periodic"
pad <-  "zero"
if(wf!="haar") v <- as.integer(as.numeric(substr(wf,2,3))/2) else v <- 1

# Maximum decomposition level J
n <- length(x)
J <- ceiling(log(n/(2*v-1))/log(2)) #(Kaiser, 1994)
if(wf=="haar"&&mode=="MODWT") J = J-1 #since modwt no need a dyadic number size

#variance transfrom
if(mode=="MRA"){
  dwt.list<- lapply(data.list, function(x) dwt.vt(x, wf, J, method, pad, boundary, cov.opt))
} else if(mode=="MODWT") {
  dwt.list<- lapply(data.list, function(x) modwt.vt(x, wf, J, pad, boundary,vt.opt, cov.opt))
} else {
  dwt.list<- lapply(data.list, function(x) at.vt(x, wf, pad, boundary, cov.opt))
}


#------------------------------------------------------
#plot before and after vt
for(j in 1:length(dwt.list)){

  dwt <- dwt.list[[j]]

  par(mfrow=c(ncol(dp)+1,1),mar=c(2,4,2,2),bg = "transparent",pty="m",ps=8)
  #ts.plot(x,xlab=NA, main=paste0("Sampled Grid: ",Grid), ylab=paste0("SPI",i), col=c("black"),lwd=c(1))
  ts.plot(x,xlab=NA, ylab=paste0("SPI",i), col=c("black"),lwd=c(1))
  for(nc in 1:ncol(dp))
    ts.plot(dwt$dp[,nc],dwt$dp.n[,nc],xlab=NA,ylab=paste0(CI.names[nc]),
            col=c("red","blue"),lwd=c(1,1),lty=c(1,2))

  p.list[[length(p.list)+1]] <- recordPlot()

  # plot.ts(cbind(dwt$x,dwt$dp))
  # plot.ts(cbind(dwt$x,dwt$dp.n))

}

}
#------------------------------------------------------
fig <- cowplot::plot_grid(plotlist = p.list, ncol=2, labels = label, label_size = 12, hjust=0)
fig

if(flag.plt){
  if(mode=="MRA") ggsave(paste0("Figure_VT_",mode,"_",cov.opt,"_",method,"_",wf,"_",Grid,".jpg"),
                         width = 200, height = 230, units = "mm", dpi = 600, fig)
  if(mode=="MODWT") ggsave(paste0("Figure_VT_",mode,"_",cov.opt,"_",vt.opt,"_",wf,"_",Grid,".jpg"),
                           width = 200, height = 230, units = "mm", dpi = 600, fig)
  if(mode=="a trous") ggsave(paste0("Figure_VT_",mode,"_",cov.opt,"_",wf,"_",Grid,".jpg"),
                             width = 200, height = 230, units = "mm", dpi = 600, fig)
}




