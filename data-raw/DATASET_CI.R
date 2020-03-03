## code to prepare Climate Index dataset goes here
# rm(list=ls()) # remove all variables
# graphics.off() # remove all figures

library(zoo)
library(R.matlab)
library(tidyverse)
#-----------------------------------------------------------------
###Nino3.4, SOI, PDO, and DMI: 1900-2010 monthly data, 1332 data points
CI.obs <- readMat("./data-raw/Climate_Indices_1900_2010.mat")

CI.obs.df <- data.frame(CI.obs$CI)
names(CI.obs.df) <- c("Year", "Month", "SOI",  "Nino",  "PDO",  "DMI")
CI.obs.ts <-  ts(CI.obs.df[,3:ncol(CI.obs.df)],start=c(1900,1), frequency = 12)
start(CI.obs.ts); end(CI.obs.ts)
#-----------------------------------------------------------------
###SAM: 1851-2011 monthly data, 1932 data points
filename <- "./data-raw/sam.20crv2c.long.data.txt"

sam <- read.table(filename,skip=1,nrows=161)
names(sam) <- c("Year",1:12)
summary(sam)

sam.df <- gather(sam,"Month","SST",2:13)
sam.df <- mutate(sam.df, Date=as.Date(paste0(as.character(Year),"-",Month,"-15"),format="%Y-%m-%d"))

sam.zoo <- read.zoo(sam.df[c("Date","SST")])
start(sam.zoo);end(sam.zoo)
summary(sam.zoo);plot(sam.zoo)

sam.ts <- ts(sam.zoo, start=c(1851,1),freq=12)
SAM <- window(sam.ts, start=c(1900,1),end=c(2010,12))
#-----------------------------------------------------------------
###TPI: 1877-2012 monthly data, 1632 data points
filename1 <- "./data-raw/TPI.hadisst11.filt.data.txt"

tpi <- read.table(filename1,skip=8,nrows=136)
names(tpi) <- c("Year",1:12)
summary(tpi)

tpi.df <- gather(tpi,"Month","SST",2:13)
tpi.df <- mutate(tpi.df, Date=as.Date(paste0(as.character(Year),"-",Month,"-15"),format="%Y-%m-%d"))

tpi.zoo <- read.zoo(tpi.df[c("Date","SST")])
start(tpi.zoo);end(tpi.zoo)
summary(tpi.zoo);plot(tpi.zoo)

tpi.ts <- ts(tpi.zoo, start=c(1877,1),freq=12)
TPI <- window(tpi.ts, start=c(1900,1),end=c(2010,12))

#-----------------------------------------------------------------
#combine and save data
data.CI <- do.call("cbind", data.frame(CI.obs.ts,SAM,TPI))
summary(data.CI)
start(data.CI); end(data.CI)
use_data(data.CI, overwrite = TRUE)

