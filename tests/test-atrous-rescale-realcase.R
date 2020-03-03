#--------------------------------------------------------------------------------------------
#sensitivity ananlysis of sign of covariance matrix with real-world data
rm(list=ls()) # remove all variables
graphics.off()

library(WASP)
library(ggplot2)

at.vt1 <- function(data, wf, pad, boundary, cov.opt=c("auto","pos","neg")){

  # initialization
  x= data$x; dp= data$dp
  mu.dp <- apply(dp,2,mean)
  if(wf!="haar") v <- as.integer(as.numeric(substr(wf,2,3))/2) else v <- 1

  # variance transfrom
  ndim=ncol(dp); n=nrow(dp); J <- ceiling(log2(n))
  S <- matrix(nrow=J+1, ncol=ndim)
  dp.n <- matrix(nrow=n,ncol=ndim)

  for(i in 1:ndim){

    # Padding
    dp.c <- scale(dp[,i],scale=F)
    pp <- padding(dp.c, pad=pad)

    # Multiresolution Analysis
    coefs <- at.wd(pp, v, bc=boundary)
    B <- coefs[1:n,]

    Bn <- scale(B)
    Bn[is.na(Bn)] <- 0 # for haar
    V <- as.numeric(apply(B,2,sd))

    dif <- sum(abs(Bn%*%V+mu.dp[i]-dp[,i]));
    if(dif>10^-10) warning(paste0("Difference between Reconstructed and original:",dif))

    cov <- cov(x, Bn)
    if(cov.opt=="pos") cov <- cov else if(cov.opt=="neg") cov <- -cov
    S[,i] <- as.vector(cov)

    #Vr <- as.numeric(cov/norm(cov,type="2")*sd(dp.c)) # For variance decompoistion
    Vr <- as.numeric(cov/norm(cov,type="2")*norm(V,type="2")) #For a trous!

    #dp.n[,i] <- scale(Bn%*%Vr)*sd(dp[,i]) + mu.dp[i]
    dp.n[,i] <- Bn%*%Vr + mu.dp[i]

    #check the correlation after vt then decide the direction of C
    if(cov.opt=="auto"){
      #if(cor(dp.n[,i],dp[,i])<0) cov <- -cov
      if(cor.test(dp.n[,i],dp[,i])$estimate<0&cor.test(dp.n[,i],dp[,i])$p.value<0.05) cov <- -cov

      Vr <- as.numeric(cov/norm(cov,type="2")*norm(V,type="2")) #For a trous!
      #dp.n[,i] <- scale(Bn%*%Vr)*sd(dp[,i]) + mu.dp[i]
      dp.n[,i] <- Bn%*%Vr + mu.dp[i]

    }

    dif.var <- (var(dp[,i])-var(dp.n[,i]))/var(dp[,i])
    if(dif.var>0.15) warning(paste0("Variance difference between Reconstructed and original(percentage):",dif.var*100))

  }

  dwt <- list(wavelet = wf,
              J = J,
              boundary = boundary,
              pad = pad,

              x=x,
              dp=dp,
              dp.n=dp.n,

              S=S
  )
  class(dwt) <- "atrous"

  return(dwt)

}
#--------------------------------------------------------------------------------------------
### wavelet method selection
mode <- switch(3,"MRA", "MODWT","a trous")
if(mode=="MRA") method <- switch(1,"dwt","modwt")
if(mode=="MODWT") vt.opt <- switch(2,"Sxx","Cov")

# wavelet family, extension mode and package
wf <- "d4" # wavelet family D8 or db4
pad <-  "zero"
boundary <- "periodic"
if(wf!="haar") v <- as.integer(as.numeric(substr(wf,2,3))/2) else v <- 1

flag.plt = 1
label = c("(a)","(b)","(c)")
path <- "~/OneDrive - UNSW/R_Package/WASP/tests/"
#--------------------------------------------------------------------------------------------
opts <- c("auto","pos","neg")[1]
p1.list <- list()
for(cov.opt in opts){
#cov.opt = "auto"
###Real-world example
SPI.12 <- SPI.calc(rain.mon,sc=12)
x <- window(SPI.12,start=c(1950,1),end=c(2009,12))
dp <- window(obs.mon,start=c(1950,1),end=c(2009,12))
lab.names <- colnames(obs.mon)

station.id = c(5) # station to investigate
data.list <- lapply(station.id, function(id) list(x=x[,id],dp=dp))

#Maximum decomposition level J
n <- nrow(x)
J <- ceiling(log(n/(2*v-1))/log(2)) #(Kaiser, 1994)
if(wf=="haar"&&mode=="MODWT") J = J-1 #since modwt no need a dyadic number size

#variance transfrom - calibration
if(mode=="MRA"){
  dwt.list<- lapply(data.list, function(x) dwt.vt(x, wf, J, method, pad, boundary, cov.opt))
} else if(mode=="MODWT") {
  dwt.list<- lapply(data.list, function(x) modwt.vt(x, wf, J, pad, boundary,vt.opt, cov.opt))
} else {
  dwt.list<- lapply(data.list, function(x) at.vt(x, wf, pad, boundary, cov.opt))
}

for(j in 1:length(dwt.list)){

  dwt <- dwt.list[[j]]

  par(mfrow=c(ncol(dp),1),mar=c(0,3,2,2),
      oma = c(2, 1, 0, 0), # move plot to the right and up
      mgp=c(2, 1, 0), # move axis labels closer to axis
      pty="m",bg = "transparent",
      ps=8)
  for(i in 1:ncol(dp))
    ts.plot(dwt$dp[,i],dwt$dp.n[,i],xlab=NA,ylab=paste0(lab.names[i]),
            col=c("black","red"),lwd=c(1,2))

  p1.list[[length(p1.list)+1]] <- recordPlot()

}

}
#--------------------------------------------------------------------------------------------
for(cov.opt in opts){
  #cov.opt = "auto"
  ###Real-world example
  SPI.12 <- SPI.calc(rain.mon,sc=12)
  x <- window(SPI.12,start=c(1950,1),end=c(2009,12))
  dp <- window(obs.mon,start=c(1950,1),end=c(2009,12))
  lab.names <- colnames(obs.mon)

  station.id = c(5) # station to investigate
  data.list <- lapply(station.id, function(id) list(x=x[,id],dp=dp))

  #Maximum decomposition level J
  n <- nrow(x)
  J <- ceiling(log(n/(2*v-1))/log(2)) #(Kaiser, 1994)
  if(wf=="haar"&&mode=="MODWT") J = J-1 #since modwt no need a dyadic number size

  #variance transfrom - calibration
  if(mode=="MRA"){
    dwt.list<- lapply(data.list, function(x) dwt.vt(x, wf, J, method, pad, boundary, cov.opt))
  } else if(mode=="MODWT") {
    dwt.list<- lapply(data.list, function(x) modwt.vt(x, wf, J, pad, boundary,vt.opt, cov.opt))
  } else {
    dwt.list<- lapply(data.list, function(x) at.vt1(x, wf, pad, boundary, cov.opt))
  }

  for(j in 1:length(dwt.list)){

    dwt <- dwt.list[[j]]

    par(mfrow=c(ncol(dp),1),mar=c(0,3,2,2),
        oma = c(2, 1, 0, 0), # move plot to the right and up
        mgp=c(2, 1, 0), # move axis labels closer to axis
        pty="m",bg = "transparent",
        ps=8)
    for(i in 1:ncol(dp))
      ts.plot(dwt$dp[,i],dwt$dp.n[,i],xlab=NA,ylab=paste0(lab.names[i]),
              col=c("black","red"),lwd=c(1,2))

    p1.list[[length(p1.list)+1]] <- recordPlot()

  }

}
#-----------------------------
fig1 <- cowplot::plot_grid(plotlist =p1.list[1:2], ncol=2, labels=label[1:2])
fig1

# fig1 <- cowplot::plot_grid(plotlist =p1.list[1:3], ncol=3, labels=label)
# fig1

if(flag.plt){

  if(mode=="a trous") ggsave(paste0(path,"Figure_RealCase_rescale_",mode,"_",wf,"_",J,".jpg"),
                             width = 240, height = 160, units = "mm", dpi = 500, fig1)
}
