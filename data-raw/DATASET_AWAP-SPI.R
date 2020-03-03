#---------------------------------------------------------------------------------------------
rm(list=ls()) # remove all variables
graphics.off() # remove all figures

library(zoo)
library(R.matlab)
library(devtools)

library(rgdal)
library(pracma)
library(raster)

flag.save=0
#---------------------------------------------------------------------------------------------
###import basemap
map <- "./data-raw/AUS_boundary_simplified.shp"
Aus_map <- readOGR(map)

#--------------------------------------------------
###AWAP rainfall 0.5 deg: 1900 - 2009 monthly data
AWAP.mon.rain.mat <- readMat("./data-raw/AWAP_MonthlyRainfall_89X70_0.50.mat")

AWAP.mon.rain <- AWAP.mon.rain.mat$precp
AWAP.mon.rain.lat <- AWAP.mon.rain.mat$lat
AWAP.mon.rain.lon <- AWAP.mon.rain.mat$lon
AWAP.mon.rain.time <- AWAP.mon.rain.mat$time
head(AWAP.mon.rain.time);tail(AWAP.mon.rain.time)

dim.AWAP = dim(AWAP.mon.rain)
data.AWAP = t(pracma::Reshape(AWAP.mon.rain, dim.AWAP[1]*dim.AWAP[2], dim.AWAP[3]))

#lon and lat
lat_lon <- readMat("./data-raw/lat_lon.mat")
lat_mat <- lat_lon$latitude
lon_mat <- lat_lon$longitude

as.vector(AWAP.mon.rain.lat-unique(lat_mat))
as.vector(AWAP.mon.rain.lon-unique(lon_mat))

#NaN grids check
Ind_AWAP <- which(apply(data.AWAP, 2, function(m) sum(is.na(m)))==0)
#Ind_AWAP_NaN <- which(apply(data.AWAP, 2, function(m) sum(is.na(m)))==nrow(data.AWAP))
#length(Ind_AWAP)+length(Ind_AWAP_NaN)

Ind_AWAP_NaN1 <- which(apply(data.AWAP, 2, function(m) sum(is.na(m)))==nrow(data.AWAP))
Ind_AWAP_NaN2 <- which(apply(data.AWAP, 2, function(m) sum(m==0))>0.3*nrow(data.AWAP))
Ind_AWAP_NaN <- c(Ind_AWAP_NaN1,Ind_AWAP_NaN2)

Ind_NaN.mat <- readMat("./data-raw/Ind_NaN.mat")
Ind_NaN <- as.integer(Ind_NaN.mat$Ind.NaN)

as.vector(Ind_NaN-Ind_AWAP_NaN)

data.AWAP[,Ind_AWAP_NaN] <- NaN
#---------------------------------------------------------------------------------------------
###Aggregate to  2.5 deg
if(TRUE){

  #matrix to spatial points
  AWAP.obs.p <- data.frame(lon=lon_mat,lat=lat_mat, t(data.AWAP))

  coordinates(AWAP.obs.p) = ~lon+lat
  proj4string(AWAP.obs.p) = "+proj=longlat +datum=WGS84"
  gridded(AWAP.obs.p) <- TRUE

  AWAP.obs.ras <- brick(AWAP.obs.p) #raster() covert only one layer/time step
  AWAP.obs.ras

  plot(AWAP.obs.ras[[1]], main="AWAP Rainfall", ylim=c(-45,-10),xlim=c(110,160))
  plot(Aus_map, col=NA, add=TRUE)

  #aggregate to 2.5 degree
  AWAP.ras.2.5 <- raster::aggregate(AWAP.obs.ras, fact=5, fun=mean,na.rm=T)
  AWAP.ras.2.5

  plot(AWAP.ras.2.5[[1]], main="AWAP Rainfall", ylim=c(-45,-10),xlim=c(110,160))
  plot(Aus_map, col=NA, add=TRUE)

  tmp <- as.data.frame(AWAP.ras.2.5, xy=TRUE)
  lat_lon.2.5 <- data.frame(lon= tmp$x, lat= tmp$y)
  data.AWAP.2.5 <- matrix(t(subset(tmp, select=-c(x,y))), ncol=nrow(tmp))

}
summary(data.AWAP.2.5)
#---------------------------------------------------------------------------------------------
###cross check and save
#identify NaN grid
# Ind_AWAP.2.5 <- which(apply(data.AWAP.2.5, 2, function(m) sum(is.na(m)))==0)
Ind_AWAP.2.5_NaN1 <- which(apply(data.AWAP.2.5, 2, function(m) sum(is.na(m)))==nrow(data.AWAP.2.5))
Ind_AWAP.2.5_NaN2 <- which(apply(data.AWAP.2.5, 2, function(m) sum(m==0))>=0.25*nrow(data.AWAP.2.5))
Ind_AWAP.2.5_NaN <- c(Ind_AWAP.2.5_NaN1,Ind_AWAP.2.5_NaN2)

Ind_AWAP.2.5 <- seq(1, ncol(data.AWAP.2.5))[-Ind_AWAP.2.5_NaN]
length(Ind_AWAP.2.5_NaN)+length(Ind_AWAP.2.5)

summary(Ind_AWAP.2.5)
#---------------------------------------------------------------------------------------------
for(sc in c(12,36)){
  start.calc = switch(as.character(sc),"12"=c(1909,1),"36"=c(1907,1))
  end.calc = c(2009,12)
  start.yr = c(1910,1); end.yr = c(2009,12)

  #SPI calculation
  AWAP.mon.ts <- window(ts(data.AWAP.2.5, start=c(1900,1), frequency = 12),
                        start=start.calc, end=end.calc)

  SPI  <- matrix(NA, nrow=nrow(AWAP.mon.ts), ncol=ncol(AWAP.mon.ts))
  SPI[,Ind_AWAP.2.5] <- SPI.calc(AWAP.mon.ts[,Ind_AWAP.2.5], sc=sc, method="mle")

  SPI.ts <- ts(SPI, start=start.calc, frequency = 12)
  start(SPI.ts); end(SPI.ts)

  assign(paste0("SPI.",sc), window(SPI.ts,start=start.yr, end=end.yr))

}
summary(SPI.12[,Ind_AWAP.2.5])

#---------------------------------------------------------------------------------------------
#save
if(flag.save){
  #AWAP
  use_data(lat_lon.2.5, overwrite = TRUE)
  use_data(Ind_AWAP.2.5, overwrite = TRUE)
  use_data(data.AWAP.2.5, overwrite = TRUE)

  #SPI
  for(sc in c(12,36)){
    do.call("use_data", list(as.name(paste0("SPI.",sc)), overwrite=TRUE))
  }

}
