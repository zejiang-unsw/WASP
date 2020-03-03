## code to prepare Synthetic dataset goes here
# rm(list=ls()) # remove all variables
# graphics.off() # remove all figures

#------------------------------------------------------------------------------
####synthetic example - Rossler
sample=50000*2
s=c(0.1,0.5,1.0) # noise level
ts.list <- list()
for(i in 1:length(s)){

  ts.r <- data.gen.Rossler(a = 0.2, b = 0.2, w = 5.7, start = c(-2, -10, 0.2),
                           time = seq(0, 50, length.out = sample))

  #add noise
  ts.r$x <- ts.r$x + rnorm(sample,mean=0, sd=s[i])
  ts.r$y <- ts.r$y + rnorm(sample,mean=0, sd=s[i])
  ts.r$z <- ts.r$z + rnorm(sample,mean=0, sd=s[i])

  ts.list[[i]] <- ts.r

  #ts.plot(ts(ts.r$x),ts(ts.r$y),ts(ts.r$z), col=c("black","red","blue"))

}
names(ts.list) <- paste0("s=",s)
assign("data.Rossler",ts.list)
use_data(data.Rossler, overwrite = TRUE)

#------------------------------------------------------------------------------
###synthetic example - Sinewave model
#frequency, sampled from a given range
fd <- c(3,5,10,15,25,30,55,70,95)
s=c(0.1,0.5,1.0)
data.SW1.list <- list()
for(i in 1:length(s)){

  data.SW1 <- data.gen.SW(nobs=512,fp=25,fd=fd,sd.x=s[i],sd.y=s[i])
  data.SW1.list[[i]] <- data.SW1

  # plot.ts(cbind(data.SW1$x, data.SW1$dp))
  # plot(data.SW1$dp[,data.SW1$true.cpy],data.SW1$x)

}
names(data.SW1.list) <- paste0("s=",s)
assign("data.SW1",data.SW1.list)
use_data(data.SW1, overwrite = TRUE)

#frequency, sampled from a given range
fd <- c(3,5,10,15,25,30,55,70,95)
s=c(1.0,2.0,3.0)
data.SW3.list <- list()
for(i in 1:length(s)){

  data.SW3 <- data.gen.SW(nobs=512,fp=c(15,25,30),fd=fd,sd.x=s[i],sd.y=s[i])
  data.SW3.list[[i]] <- data.SW3

  # plot.ts(cbind(data.SW3$x, data.SW3$dp))
  # plot(data.SW3$dp[,data.SW3$true.cpy[1]],data.SW3$x)

}
names(data.SW3.list) <- paste0("s=",s)
assign("data.SW3",data.SW3.list)
use_data(data.SW3, overwrite = TRUE)

#------------------------------------------------------------------------------
###synthetic example - Hysteresis loop
#frequency, sampled from a given range
fd <- c(3,5,10,15,25,30,55,70,95)
s=c(0.1,0.5,1.0)
data.HL.list <- list()
for(i in 1:length(s)){

  data.HL <- data.gen.HL(n=3,m=5,nobs=512,fp=25,fd=fd,sd.x=s[i],sd.y=s[i])
  data.HL.list[[i]] <- data.HL

  # plot.ts(cbind(data.HL$x, data.HL$dp))
  # plot(data.HL$dp[,data.HL$true.cpy],data.HL$x)

}
names(data.HL.list) <- paste0("s=",s)
assign("data.HL",data.HL.list)
use_data(data.HL, overwrite = TRUE)
