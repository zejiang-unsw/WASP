#--------------------------------------------------------------------------------------------
#sensitivity ananlysis of options of covariance matrix with Rossler model
rm(list=ls()) # remove all variables
graphics.off()

library(WASP)
library(biwavelet)
library(FNN)
library(ggplot2)
library(tidyr)
library(dplyr)

#--------------------------------------------------------------------------------------------
### wavelet method selection
mode <- "MODWT"
cov.opt <- switch(1,"auto","pos","neg")
#if(mode=="MODWT") vt.opt <- ifelse(1,"Sxx","Cov")

#data sample
s=0
sample=5000
dt=NULL
#set.seed(101)
ts.list <- list(data.gen.Rossler(a = 0.2, b = 0.2, w = 5.7, start = c(-2, -10, 0.2),
                                 time = seq(0, 50, length.out = sample)))

#knn prediction
k=5

# wavelet family, extension mode and package
wf <- "d4" # wavelet family D8 or db4
boundary <- "periodic"
pad <-  "zero"

#Maximum decomposition level J
if(wf!="haar") v <- as.integer(as.numeric(substr(wf,2,3))/2) else v <- 1
n <- sample
J <- ceiling(log(n/(2*v-1))/log(2)) #(Kaiser, 1994)
if(wf=="haar"&&mode=="MODWT") J = J-1 #since modwt no need a dyadic number size

flag.plt = 0
label = c("(a)","(b)")
path <- "~/OneDrive - UNSW/R_Package/WASP/tests/"
#--------------------------------------------------------------------------------------------
if(TRUE){
  p1.list <- list(); p2.list <- list()
  max=3.5

  par(mfrow=c(1,1),
      mar=c(2,4,2,3),
      oma = c(0,0,0,3),
      mgp=c(1.5, 0.8, 0),
      bg = "transparent",pty="m",
      ps=8)

  x=ts.list[[1]]$z
  x.ts <- cbind(1:length(x), x)
  wt.t1 <- biwavelet::wt(x.ts, pad=TRUE, mother="morlet",dt=dt)

  plot(wt.t1, plot.cb = TRUE, plot.phase = FALSE, xlab=NA, xlim=c(0,sample),
       legend.loc=c(.9, .92, .2, .8))

  p1.list[[length(p1.list)+1]] <- recordPlot()

  modwt.x <- waveslim::modwt(x, wf = wf, n.levels = J, boundary = boundary)
  Sxx <- unlist(lapply(modwt.x, var))


  x=ts.list[[1]]$x
  x.ts <- cbind(1:length(x), x)
  wt.t1 <- biwavelet::wt(x.ts, pad=TRUE, mother="morlet",dt=dt)

  plot(wt.t1, plot.cb = TRUE, plot.phase = FALSE, xlab=NA, xlim=c(0,sample),
       legend.loc=c(.9, .92, .2, .8))

  p1.list[[length(p1.list)+1]] <- recordPlot()

  modwt.x <- waveslim::modwt(x, wf = wf, n.levels = J, boundary = boundary)
  Sxx.dp <- unlist(lapply(modwt.x, var))

  print(sum(Sxx.dp)-var(x))

}

barplot(cbind(Sxx,Sxx.dp), beside=TRUE)

#--------------------------------------------------------------------------------------------
opts <- c("Sxx","Cov")
RMSE <- NULL; p.list <- list();Sxx.dp.n<- NULL
for(vt.opt in opts){
#vt.opt = "Sxx"

data.list <- lapply(ts.list, function(ts) list(x=ts$z,dp=cbind(ts$x,ts$y)))

#variance transfrom - calibration
if(mode=="MRA"){
  dwt.list<- lapply(data.list, function(x) dwt.vt(x, wf, J, method, pad, boundary, cov.opt))
} else if(mode=="MODWT") {
  dwt.list<- lapply(data.list, function(x) modwt.vt(x, wf, J, pad, boundary,vt.opt, cov.opt))
} else {
  dwt.list<- lapply(data.list, function(x) at.vt(x, wf, pad, boundary, cov.opt))
}

#-------------------------------------------------------------------------------
for(j in 1:length(dwt.list)){

  dwt <- dwt.list[[j]]

  #layout(rbind(matrix(1:4,nrow=2),c(5,5)))
  par(mfrow=c(2,1),
      mar=c(2,3,2,2),
      oma = c(1, 1, 1, 0),
      mgp=c(2, 1, 0),
      bg = "transparent",pty="m",
      ps=8)
  len=5
  for(i in 1:2){
    ts.plot(ts(dwt$dp[,i]),ts(dwt$dp.n[,i]),ylim=c(-15,15), xlab=NA,ylab=paste0("X",i),
            col=c("black","red"),lwd=c(2,2))
    legend("topleft", legend=c("Original", "Variance-transformed"),box.lty=0,bg='transparent',
           col=c("black", "red"), lty=c(1,1), seg.len=len)
  }

  p.list[[length(p.list)+1]] <- recordPlot()
}


#-------------------------------------------------------------------------------
#######after variance transform
x.n= dwt.list[[1]]$dp.n[,1]
x.n.ts <- cbind(1:length(x.n), x.n)
wt.t2 <- biwavelet::wt(x.n.ts, pad=TRUE, mother="morlet",dt=dt)

plot.ts(x.n)

modwt.x <- waveslim::modwt(x.n, wf = wf, n.levels = J, boundary = boundary)
Sxx.dp.n <- cbind(Sxx.dp.n, unlist(lapply(modwt.x, var)))

print(sum(unlist(lapply(modwt.x, var)))-var(x.n))

## Plot power
## Make room to the right for the color bar
#jpeg(file = "predictor_vt.jpeg",width=320,height=160, units="mm",res = 500)

par(mfrow=c(1,1),
    mar=c(2,4,2,3),
    oma = c(0,0,0,3),
    mgp=c(1.5, 0.8, 0),
    bg = "transparent",pty="m",
    ps=8)
# plot(wt.t1, plot.cb = TRUE, plot.phase = FALSE, xlab=NA,
#      legend.loc=c(.9, .92, .2, .8))
plot(wt.t2, plot.cb = TRUE, plot.phase = FALSE, xlab=NA, xlim=c(0,sample),
     legend.loc=c(.9, .92, .2, .8))

p1.list[[length(p1.list)+1]] <- recordPlot()

#dev.off()

print(sum(unlist(lapply(modwt.x, var)))-sum(Sxx.dp))

}


#-----------------------------
barplot(cbind(Sxx, Sxx.dp, Sxx.dp.n), beside=TRUE)

fig1 <- cowplot::plot_grid(plotlist = p.list, labels=label)
fig1

fig2 <- cowplot::plot_grid(plotlist = p1.list, ncol=2, labels=c("(a)","(b)","(c)","(d)"))
#fig2

fig3 <- cowplot::plot_grid(plotlist = p1.list[c(1,2,4)], ncol=1, labels=c("(a)","(b)","(c)"))
#fig3

if(flag.plt){
  if(mode=="MRA") ggsave(paste0(path,"Figure_SynCase_",mode,"_",method,"_",wf,"_",J,".jpg"),
                         width = 230, height = 115, units = "mm", dpi = 500, fig1)
  if(mode=="MODWT") ggsave(paste0(path,"Figure_SynCase_",mode,"_vt.opt_",wf,"_",J,".jpg"),
                           width = 230, height = 115, units = "mm", dpi = 500, fig1)
  if(mode=="a trous") ggsave(paste0(path,"Figure_SynCase_",mode,"_",wf,"_",J,".jpg"),
                             width = 230, height = 115, units = "mm", dpi = 500, fig1)


  if(mode=="MRA") ggsave(paste0(path,"Figure_SynCase_power_",mode,"_",method,"_",wf,"_",J,".jpg"),
                         width = 240, height = 120, units = "mm", dpi = 500, fig2)
  if(mode=="MODWT") ggsave(paste0(path,"Figure_SynCase_power_",mode,"_vt.opt_",wf,"_",J,".jpg"),
                           width = 240, height = 120, units = "mm", dpi = 500, fig2)
  if(mode=="a trous") ggsave(paste0(path,"Figure_SynCase_power_",mode,"_",wf,"_",J,".jpg"),
                             width = 240, height = 120, units = "mm", dpi = 500, fig2)

  if(mode=="MRA") ggsave(paste0(path,"Figure_SynCase_ps_",mode,"_",method,"_",wf,"_",J,".jpg"),
                         width = 120, height = 160, units = "mm", dpi = 500, fig3)
  if(mode=="MODWT") ggsave(paste0(path,"Figure_SynCase_ps_",mode,"_vt.opt_",wf,"_",J,".jpg"),
                           width = 120, height = 160, units = "mm", dpi = 500, fig3)
  if(mode=="a trous") ggsave(paste0(path,"Figure_SynCase_power_",mode,"_",wf,"_",J,".jpg"),
                             width = 120, height = 160, units = "mm", dpi = 500, fig3)
}






