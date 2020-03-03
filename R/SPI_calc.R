#-------------------------------------------------------------------------------
#' Calculate Standardized Precipitation Index, SPI
#' @param prec.zoo		  A zoo series contain date and rainfall vector/matrix.
#' @param sc		        The accumulation period in months. Commonly 6, 12, 24, 36, and 48 months.
#' @param method        A character string coding for the fitting method: "mle" for 'maximum likelihood estimation', "mme" for 'moment matching estimation', "qme" for 'quantile matching estimation' and "mge" for 'maximum goodness-of-fit estimation'.

#' @return A matrix of time sereis.
#' @import fitdistrplus
#' @import zoo
#' @export
#'
#' @examples
#' data(rain.mon)
#' data(SPI.12)
#' ##compute SPI
#' SPI <- SPI.calc(window(rain.mon, start=c(1949,1), end=c(2009,12)),sc=12)
#'
#' ##compare with sample data set SPI.12
#' par(mfrow=c(3,5))
#' for(i in 1:ncol(SPI)) plot(SPI[,i],SPI.12[,i])

SPI.calc <- function(prec.zoo,sc=24, method="mme"){
###spi calculation

### Input
# prec.zoo is the zoo series contain date and rainfall vector/matrix and
# sc is the accumulation period in months

### Computation
# Calculate running sum of different duration like 6,12 24,48 etc)
# rollsumr(..., align = "right"), rank is applied to zoo series, align which is the date you took
#run_sum <- apply(prec,2,rollsum,k=sc) # if prec is a rainfall matrix without date
run_sum <- rollsumr(prec.zoo,k=sc) # if prec is a zoo objective
# number of data after rollsum, length(date) -sc + 1 = nrow(run_sum)

### Estimate shape and rate parameter of gamma distribution from data
alpha<- NA ;beta<- NA
ncol<- ifelse(is.null(ncol(run_sum)),1,ncol(run_sum))
nrow <- ifelse(is.null(ncol(run_sum)),length(run_sum),nrow(run_sum))

for (i in 1:ncol){
  # "mle" for 'maximum likelihood estimation'
  # lower: lower limit of the parameters: shape>0 and rate>0 by definition
  pd <-  fitdist(as.numeric(run_sum[ ,i]),"gamma", method = method, lower = c(0, 0))
  alpha[i]= pd$estimate[1] # shape
  beta[i] = pd$estimate[2] # rate
}
### fit gamma
gammavals<- matrix(NA, nrow = nrow, ncol = ncol)

for(i in 1: ncol){
  # pgamma gives the cumulative distribution function
  gammavals[ ,i] <- pgamma(as.numeric(run_sum[ ,i]),shape=alpha[i], rate=beta[i])
}

### Inverse Normal Distribution: qnorm gives the quantile(value) function
# cumulative probability function: pnorm(1.64)=0.95 --> qnorm(0.95) = 1.64;
# dnorm():probability density/mass function
spi<- qnorm(gammavals)

### create ts spi
spi.ts <- ts(rbind(matrix(NA,nrow=sc-1,ncol=ncol(prec.zoo)),spi),start=start(prec.zoo), end=end(prec.zoo), freq=12)
return(spi.ts)

}
