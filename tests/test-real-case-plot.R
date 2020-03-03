################################################################
rm(list=ls()) # remove all variables
graphics.off() # remove all figures

setwd("C:/Users/Ze/OneDrive - UNSW/R_Package/WASP/tests/realcase")

library(WASP)
library(zoo)
library(waveslim)
library(wavethresh)
library(fitdistrplus)

library(cowplot)
library(tidyverse)

library(verification) #model performance metrics
library(overlapping) #PDF skill scores

library(sp)
library(raster)
library(rgdal)
library(RColorBrewer)
#-----------------------------------------------------------------------------------------------------
#predictor - obs
data("data.CI");Ind_CI <- colnames(data.CI)
data("SPI.12");data("SPI.36")
data("Ind_AWAP.2.5")
data("lat_lon.2.5"); lon = lat_lon.2.5$lon; lat=lat_lon.2.5$lat

sc.list = c(12,36)
k.folds = 2 #cross-validation partition

flag.plt = 1
Grids = Ind_AWAP.2.5
Grids = c(130, 87, 149, 28)
start.yr <- c(1910,1); end.yr <- c(2009,12)

### wavelet method selection - for proposed method
mode <- switch(2,"MRA", "MODWT","a trous")
cov.opt <- switch(3,"pos","neg","auto")

if(mode=="MRA") method <- switch(1,"dwt","modwt")
if(mode=="MODWT") vt.opt <- switch(2,"Sxx","Cov")

# wavelet family, extension mode and package
wf <- "d4" # wavelet family D8 or db4
boundary <- "periodic"
pad <-  "zero"
if(wf!="haar") v <- as.integer(as.numeric(substr(wf,2,3))/2) else v <- 1

#-----------------------------------------------------------------------------------------------------
#import basemap
map <- "../AUS_boundary_simplified.shp"
Aus_map <- readOGR(map)

p0.list <- list() #selection map
p1.list <- list() #PDF skill score map
p2.list <- list() #Reduced RMSE map
for(sc in sc.list){
#response - obs
data.SPI.obs =  eval(parse(text=paste0("SPI.",sc)))

#model simulated response
if(mode!="MODWT"){
  data.SPI.mod= get(load(paste0("data.SPI.",sc,".mod.vt_",mode,"_",wf,"_",cov.opt,"_",k.folds,"folds.noEns.Rdata")))
  load(paste0("data.SPI.",sc,".ref_",mode,"_",wf,"_",cov.opt,"_",k.folds,"folds.noEns.Rdata"))
  load(paste0("data.SPI.",sc,".selection_",mode,"_",wf,"_",cov.opt,"_",k.folds,"folds.noEns.Rdata"))
  load(paste0("data.SPI.",sc,".weights_",mode,"_",wf,"_",cov.opt,"_",k.folds,"folds.noEns.Rdata"))

} else {
  data.SPI.mod= get(load(paste0("data.SPI.",sc,".mod.vt_",mode,"_",vt.opt,"_",wf,"_",cov.opt,"_",k.folds,"folds.noEns.Rdata")))
  load(paste0("data.SPI.",sc,".ref_",mode,"_",vt.opt,"_",wf,"_",cov.opt,"_",k.folds,"folds.noEns.Rdata"))
  load(paste0("data.SPI.",sc,".selection_",mode,"_",vt.opt,"_",wf,"_",cov.opt,"_",k.folds,"folds.noEns.Rdata"))
  load(paste0("data.SPI.",sc,".weights_",mode,"_",vt.opt,"_",wf,"_",cov.opt,"_",k.folds,"folds.noEns.Rdata"))

}

#-----------------------------------------------------------------------------------------------------
###selection map
#summary(sel.cv); head(sel.cv)
sel <- vector("list",length=length(Grids))
for(i in 1:length(Grids)){
  tmp1 <- NULL; tmp2 <- NULL
  for(j in 1:length(sel.cv)){
    tmp1 <- c(tmp1, as.numeric(as.character(sel.cv[[j]][[i]]$origin)))
    tmp2 <- c(tmp2, as.numeric(as.character(sel.cv[[j]][[i]]$vt)))
  }

  ifelse(!is.null(tmp1), sel[[i]]$origin <- data.frame(table(tmp1))[order(data.frame(table(tmp1))$Freq, decreasing = T),]$tmp1,
         sel[[i]]$origin <- 0)
  ifelse(!is.null(tmp2), sel[[i]]$vt <- data.frame(table(tmp2))[order(data.frame(table(tmp2))$Freq, decreasing = T),]$tmp2,
         sel[[i]]$vt <- 0)
}

summary(sel); head(sel)

#---------------------------------------------------------
#selection map for the entire period
for(i in 1){
  for(j in c("origin","vt")){
  data <- matrix(NA,nrow=1,ncol=ncol(data.SPI.obs))
  data[Grids] <- as.numeric(unlist(lapply(sel,function(ls) as.character(ls[[j]])[i])))

  ###matrix to spatial points
  SPI.obs.mat <- data.frame(t(data))

  SPI.obs.p <- data.frame(lon=lon,lat=lat, SPI.obs.mat)

  coordinates(SPI.obs.p) = ~lon+lat
  proj4string(SPI.obs.p) = "+proj=longlat +datum=WGS84"
  gridded(SPI.obs.p) <- TRUE

  SPI.obs.ras <- brick(SPI.obs.p) #raster() covert only one layer/time step
  SPI.obs.ras
  #---------------------------------------
  ## labels and color
  label <- seq(0.5,6.5,1);labelat = seq(1,length(Ind_CI),1); labeltext = Ind_CI
  #my.palette <- brewer.pal(n = length(label), name = "Blues")
  #my.palette <- topo.colors(length(label))
  #my.palette <- heat.colors(length(label))

  my.palette <- c("red","blue","green","purple","pink","black")

  df.ras.mask <- mask(SPI.obs.ras, Aus_map)
  print(data.frame(table(data[Grids])))
  p <- spplot(df.ras.mask, xlim=c(112,156), ylim=c(-44,-10),
              col.regions = my.palette,
              #cuts=3,
              main = list("Period: 1910-2009",cex=2),
              at = label, # colour breaks
              sp.layout = list("sp.polygons",Aus_map, first=FALSE, col="grey"),
              colorkey = list(space="bottom",
                              height = 0.5,
                              width = 1,
                              labels = list(at = labelat, labels = labeltext, cex=1.2)
              )
  )
  p
  p0.list[[length(p0.list)+1]] <- p
  }
}

#-----------------------------------------------------------------------------------------------------
###Improved skill scores spatial plot
SPI.ref.PDF <- matrix(NA, nrow=1,ncol=ncol(data.SPI.obs))
SPI.mod.PDF <- matrix(NA, nrow=1,ncol=ncol(data.SPI.obs))

SPI.ref.PDF[,Grids] <- sapply(Grids, function(i) overlap(list(data.SPI.obs[,i],data.SPI.ref[,i]),na.rm=T,nbins = 1024)$OV)
SPI.mod.PDF[,Grids] <- sapply(Grids, function(i) overlap(list(data.SPI.obs[,i],data.SPI.mod[,i]),na.rm=T,nbins = 1024)$OV)
SPI.mod.PDF.OL <- (SPI.mod.PDF - SPI.ref.PDF)/SPI.ref.PDF*100

summary(SPI.ref.PDF[,Grids])
summary(SPI.mod.PDF[,Grids])
summary(SPI.mod.PDF.OL[,Grids])


#---------------------------------------
for(data in c("SPI.mod.PDF.OL")){
  ###matrix to spatial points
  SPI.obs.mat <- data.frame(t(eval(parse(text = data))))

  SPI.obs.p <- data.frame(lon=lon,lat=lat, SPI.obs.mat)

  coordinates(SPI.obs.p) = ~lon+lat
  proj4string(SPI.obs.p) = "+proj=longlat +datum=WGS84"
  gridded(SPI.obs.p) <- TRUE

  SPI.obs.ras <- brick(SPI.obs.p) #raster() covert only one layer/time step
  SPI.obs.ras

  #---------------------------------------
  #plot(SPI.obs.ras, main="SPI 12M", ylim=c(-44,-10),xlim=c(112,156))
  #image(SPI.obs.ras, ylim=c(-44,-10),xlim=c(112,156))

  ## labels and color
  label <- seq(0,50,10);labelat = round(label,digits=2); labeltext = round(label,digits=2)
  #my.palette <- brewer.pal(n = length(label), name = "Blues")
  #my.palette <- topo.colors(length(label))
  #my.palette <- heat.colors(length(label))

  pal <-  colorRampPalette(c("lightblue", "blue"))
  my.palette <- pal(length(label))

  df.ras.mask <- mask(SPI.obs.ras, Aus_map)
  p <- spplot(df.ras.mask[[1]], xlim=c(112,156), ylim=c(-44,-10),
              col.regions = my.palette,
              cuts=5,
              at = labelat, # colour breaks
              sp.layout = list("sp.polygons",Aus_map, first=FALSE, col="grey"),
              auto.key = list(cex = 1.5, title="percent"),
              colorkey = list(space="bottom",
                              title= "percent", #'kg/m2',
                              height = 0.5,
                              width = 1,
                              labels = list(at = labelat, labels = labeltext, cex=1.5)
              )
  )

  p1.list[[length(p1.list)+1]] <- p

}

#-----------------------------------------------------------------------------------------------------
###Reduced RMSE spatial plot - package verification
if(TRUE){
  data.SPI.ref.RMSE  <- matrix(NA, nrow=1,ncol=ncol(data.SPI.obs))
  data.SPI.mod.RMSE <- matrix(NA, nrow=1,ncol=ncol(data.SPI.obs))

  data.SPI.ref.RMSE[,Grids] <- sapply(Grids, function(i)
    sqrt(verify(data.SPI.obs[,i],data.SPI.ref[,i], frcst.type = "cont", obs.type = "cont")$MSE))

  data.SPI.mod.RMSE[,Grids] <- sapply(Grids, function(i)
    sqrt(verify(data.SPI.obs[,i],data.SPI.mod[,i], frcst.type = "cont", obs.type = "cont")$MSE))

  data.SPI.mod.RMSE.dif <- (data.SPI.ref.RMSE - data.SPI.mod.RMSE)/data.SPI.ref.RMSE*100
}
summary(data.SPI.ref.RMSE[,Grids])
summary(data.SPI.mod.RMSE[,Grids])
summary(data.SPI.mod.RMSE.dif[,Grids])

#---------------------------------------
for(data in c("data.SPI.mod.RMSE.dif")){
###matrix to spatial points
SPI.obs.mat <- data.frame(t(eval(parse(text = data))))

SPI.obs.p <- data.frame(lon=lon,lat=lat, SPI.obs.mat)

coordinates(SPI.obs.p) = ~lon+lat
proj4string(SPI.obs.p) = "+proj=longlat +datum=WGS84"
gridded(SPI.obs.p) <- TRUE

SPI.obs.ras <- brick(SPI.obs.p) #raster() covert only one layer/time step
SPI.obs.ras

#---------------------------------------
#plot(SPI.obs.ras, main="SPI 12M", ylim=c(-44,-10),xlim=c(112,156))
#image(SPI.obs.ras, ylim=c(-44,-10),xlim=c(112,156))

## labels and color
label <- seq(0,50,10);labelat = round(label,digits=2); labeltext = round(label,digits=2)
#my.palette <- brewer.pal(n = length(label), name = "Blues")
#my.palette <- topo.colors(length(label))
#my.palette <- heat.colors(length(label))

pal <-  colorRampPalette(c("lightblue", "blue"))
my.palette <- pal(length(label))

df.ras.mask <- mask(SPI.obs.ras, Aus_map)
p <- spplot(df.ras.mask[[1]], xlim=c(112,156), ylim=c(-44,-10),
            col.regions = my.palette,
            #cuts=5,
            at = labelat, # colour breaks
            sp.layout = list("sp.polygons",Aus_map, first=FALSE, col="grey"),
            #auto.key = list(cex = 1.5, title="percent"),
            colorkey = list(space="bottom",
                            title= "percent", #'kg/m2',
                            height = 0.5,
                            width = 1,
                            labels = list(at = labelat, labels = labeltext, cex=1.5)
            )
)

p2.list[[length(p2.list)+1]] <- p

}

#------------------------------------------------------------------------------
###density plot
df.ref <- data.frame(Group=1, N=1:nrow(data.SPI.obs),
                     rbind(cbind("Obs",data.SPI.obs[,Grids]),
                           cbind("Pred",data.SPI.ref[,Grids])))
colnames(df.ref) <- c("Group","N","Type",Grids)
df.ref.n <- gather(df.ref,"No","Value",4:(length(Grids)+3)) %>% spread("Type","Value")
df.ref.n$Obs <- as.numeric(df.ref.n$Obs)
df.ref.n$Pred <- as.numeric(df.ref.n$Pred)
df.ref.n$No <- as.numeric(df.ref.n$No)
summary(df.ref.n)

df.mod <- data.frame(Group=2, N=1:nrow(data.SPI.obs),
                     rbind(cbind("Obs",data.SPI.obs[,Grids]),
                           cbind("Pred",data.SPI.mod[,Grids])))
colnames(df.mod) <- c("Group","N","Type",Grids)
df.mod.n <- gather(df.mod,"No","Value",4:(length(Grids)+3)) %>% spread("Type","Value")
df.mod.n$Obs <- as.numeric(df.mod.n$Obs)
df.mod.n$Pred <- as.numeric(df.mod.n$Pred)
df.mod.n$No <- as.numeric(df.mod.n$No)
summary(df.mod.n)

data <- rbind(df.ref.n, df.mod.n)
limits.x <- c(-3.55,3.55); breaks.x <- seq(-3,3,1)
limits.y <- c(0,1); breaks.y <- seq(0,1,0.2)
Predictor.labs <- paste0("Sampled Grid: ",Grids)
names(Predictor.labs) <- Grids

fig4 <- ggplot(data = data) +

  geom_density(aes(x=Obs, fill="Observed"),col="pink") +
  stat_density(aes(x=Pred, color=factor(Group)), geom="line",position="identity", lwd=1) +

  facet_wrap(No~., labeller = labeller(No = Predictor.labs)) +

  scale_y_continuous(breaks=breaks.y, limits=limits.y, expand = c(0,0)) +
  scale_x_continuous(breaks=breaks.x, limits=limits.x, expand = c(0,0)) +

  scale_color_manual(values=c("red","blue"), labels=c("Predicted","Predicted(VT)"))+
  scale_fill_manual(values=c("pink"))+

  theme_bw() +
  theme(text = element_text(size = 8),
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),

        legend.title = element_blank(),
        #legend.position= c(0.9,0.1),
        legend.position="bottom",
        legend.key.width = unit(1,"cm"),

        #x and y axis
        axis.title.x = element_blank()
        #axis.title.y = element_blank()

  )

fig4

#-----------------------------------------------------------------------------------------------------
###Taylor Diagram
par(mfrow=c(1,1),pty="m", ps=8)
if(TRUE){
  op=0.5
  df.obs <-  data.SPI.obs
  df.ref <-  data.SPI.ref
  df.mod <- data.SPI.mod

  taylor.diag(df.obs[,Grids[1]],df.ref[,Grids[1]],pos.cor=FALSE, cex.max = 1.2,
              main="",xlab = "",col=alpha("red",op), normalize=T,pch=16,mar=c(1,4,3,2))
  taylor.diag(df.obs[,Grids[1]],df.mod[,Grids[1]],add=TRUE,col=alpha("blue",op),normalize=T,pch=16)
  for(i in Grids[-1]){
    obs <- df.obs[,i]
    ref <- df.ref[,i]
    model <- df.mod[,i]

    taylor.diag(obs,ref,add=TRUE,col=alpha("red",op),normalize=T,pch=16)
    taylor.diag(obs,model,add=TRUE,col=alpha("blue",op),normalize=T,pch=16)

  }

  legend("topright",inset = c(0.05,0.25), box.lty=0, bg=NA,
         legend=c("Observed","Predicted","Predicted(VT)"),pch=c(15,16,16),col=c("darkgreen","red","blue"))
}
p1 <- plot_grid(recordPlot())

#-----------------------------------------------------------------------------------------------------
###Boxplot
SPI.ref.cor <- matrix(NA, nrow=1,ncol=ncol(data.SPI.obs))
SPI.mod.cor <- matrix(NA, nrow=1,ncol=ncol(data.SPI.obs))

SPI.ref.cor[,Grids] <- sapply(Grids, function(i) cor(data.SPI.obs[,i], data.SPI.ref[,i]))
SPI.mod.cor[,Grids] <- sapply(Grids, function(i) cor(data.SPI.obs[,i], data.SPI.mod[,i]))


if(TRUE){
  par(mfrow=c(1,2),mar=c(3,4,4,2),ps=8, bg="transparent")
  boxplot(cbind(t(SPI.ref.PDF),t(SPI.mod.PDF)), ylim=c(0,1),xaxt="n", cex.main=0.8,
          main="PDF Skill Score", col=c("red","blue"))
  axis(1, at= c(1,2), labels=c("Predicted","Predicted(VT)"))

  boxplot(cbind(t(SPI.ref.cor),t(SPI.mod.cor)), ylim=c(-0.5,0.5),xaxt="n", cex.main=0.8,
          main="Correlation Coefficient",col=c("red","blue"))
  axis(1, at= c(1,2), labels=c("Predicted","Predicted(VT)"))
}

p2 <- plot_grid(recordPlot())

fig <- plot_grid(p1, p2, ncol=1, labels=c("(a)","(b)"),label_size=8,
                  scale = c(1,0.9), rel_widths = c(1, 1))
fig

}
#-----------------------------------------------------------------------------------------------------
fig0 <- cowplot::plot_grid(p0.list[[1]], p0.list[[2]], ncol=2, labels=c("A","B"),
                           scale = c(1,1),label_size=24, rel_widths = c(1, 1))

fig0

fig1 <- cowplot::plot_grid(p0.list[[3]], p0.list[[4]], ncol=2, labels=c("A","B"),
                           scale = c(1,1),label_size=24, rel_widths = c(1, 1))

fig1

fig2 <- cowplot::plot_grid(p1.list[[1]], p1.list[[2]], ncol=2, labels=c("A","B"),
                           scale = c(1,1),label_size=24, rel_widths = c(1, 1))
fig2

fig3 <- cowplot::plot_grid(p2.list[[1]], p2.list[[2]], ncol=2, labels=c("A","B"),
                           scale = c(1,1),label_size=24, rel_widths = c(1, 1))
fig3

if(flag.plt){
  sc = "12-36"
  if(mode!="MODWT"){
    ggsave(paste0("SPI.12_Selection_CI_",mode,"_",wf,"_",cov.opt,"_",k.folds,"folds.r1.jpg"),
           width = 300, height = 180, units = "mm", dpi = 500, fig0)

    ggsave(paste0("SPI.36_Selection_CI_",mode,"_",wf,"_",cov.opt,"_",k.folds,"folds.r1.jpg"),
           width = 300, height = 180, units = "mm", dpi = 500, fig1)

    ggsave(paste0("SPI.",sc,"_Improved_PDF_",mode,"_",wf,"_",cov.opt,"_",k.folds,"folds.2.5.jpg"),
           width = 360, height = 180, units = "mm", dpi = 500, fig2)

    ggsave(paste0("SPI.",sc,"_Reduced_RMSE_",mode,"_",wf,"_",cov.opt,"_",k.folds,"folds.2.5.jpg"),
           width = 360, height = 180, units = "mm", dpi = 500, fig3)

  } else {
    ggsave(paste0("SPI.12_Selection_CI_",mode,"_",vt.opt,"_",wf,"_",cov.opt,"_",k.folds,"folds.r1.jpg"),
           width = 400, height = 300, units = "mm", dpi = 500, fig0)

    ggsave(paste0("SPI.36_Selection_CI_",mode,"_",vt.opt,"_",wf,"_",cov.opt,"_",k.folds,"folds.r1.jpg"),
           width = 400, height = 300, units = "mm", dpi = 500, fig1)

    ggsave(paste0("SPI.",sc,"_Improved_PDF_",mode,"_",vt.opt,"_",wf,"_",cov.opt,"_",k.folds,"folds.2.5.jpg"),
           width = 360, height = 180, units = "mm", dpi = 500, fig2)

    ggsave(paste0("SPI.",sc,"_Reduced_RMSE_",mode,"_",vt.opt,"_",wf,"_",cov.opt,"_",k.folds,"folds.2.5.jpg"),
           width = 360, height = 180, units = "mm", dpi = 500, fig3)

  }

}
