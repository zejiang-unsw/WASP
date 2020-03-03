######################################################################
rm(list=ls()) # remove all variables
#--------------------------------------------------------------------------------------------
#sensitivity ananlysis of decomposition level with Rossler model
library(WASP)
library(FNN)

library(tidyr)
library(dplyr)
library(plyr)

#--------------------------------------------------------------------------------------------
### wavelet method selection - for proposed method
mode <- switch(1,"MRA", "MODWT","a trous")
sample.data <- TRUE #use sample data in the package or not

if(mode=="MRA") method <- ifelse(1,"dwt","modwt")
if(mode=="MODWT") vt.opt <- ifelse(1,"Sxx","Cov")
#so far the case is that Sxx is good for real example, while Cov is good for Rossler

# wavelet family, extension mode and package
wf <- "haar" # wavelet family D8 or db4
boundary <- "periodic"
pad <-  "zero"
if(wf!="haar") v <- as.integer(as.numeric(substr(wf,2,3))/2) else v <- 1

#--------------------------------------------------------------------------------------------
###synthetic example Rossler
#data generation
sample=50000*2
#sample=2^16*2
sample.cal=sample/2;

k=5
s=c(0.1,0.5,1) # scaling factor for noise level
#s=c(0.1,0.5)
#s=0.1

#--------------------------------------------------------------------------------------------
### data generation or load
ts.list <- list()
for(i in 1:length(s)){
    #i <- 1
    #rossler
    ts.r <- data.gen.Rossler(a = 0.2, b = 0.2, w = 5.7, start = c(-2, -10, 0.2),
                             time = seq(0, 50, length.out = sample))

    #add noise
    ts.r$x <- ts.r$x + rnorm(length(ts.r$time),mean=0, sd=s[i])
    ts.r$y <- ts.r$y + rnorm(length(ts.r$time),mean=0, sd=s[i])
    ts.r$z <- ts.r$z + rnorm(length(ts.r$time),mean=0, sd=s[i])

    ts.list[[i]] <- ts.r

    #ts.plot(ts(ts.r$x),ts(ts.r$y),ts(ts.r$z), col=c("black","red","blue"))

}

if(sample.data){
    data("data.Rossler")
    ts.list <- data.Rossler[1:length(s)]
}


#--------------------------------------------------------------------------------------------
###proposed method
#calibration
#Maximum decomposition level J
n <- sample.cal
J <- ceiling(log(n/(2*v-1))/log(2)) #(Kaiser, 1994)
#J <- floor(log10(n)) # (Nourani et al., 2008) for standard approach
if(wf=="haar"&&mode=="MODWT") J = J-1 #since modwt no need a dyadic number size
print(paste0("Calibration: Decomposition Levels J= ",J))

#--------------------------------------------------------------------------------------------
if(TRUE){

    data.list <- lapply(ts.list, function(ts) list(x=ts$z[1:sample.cal],
                                                   dp=cbind(ts$x[1:sample.cal],ts$y[1:sample.cal])))

    #########################################################################################################
    #calibration
    #variance transfrom
    if(mode=="MRA"){
        dwt.list<- lapply(data.list, function(x) dwt.vt(x, wf, J, method, pad, boundary))
    } else if(mode=="MODWT") {
        dwt.list<- lapply(data.list, function(x) modwt.vt(x, wf, J, pad, boundary, vt.opt))
    } else {
        dwt.list<- lapply(data.list, function(x) at.vt(x, wf, pad, boundary))
    }

    #-------------------------------------------------------------------------------
    df <- NULL;data.RMSE<-NULL;dwt.RMSE<-NULL
    for(i in 1:length(dwt.list)){

        dwt <- dwt.list[[i]]
        x <- dwt$x; dp <- dwt$dp; dp.n <- dwt$dp.n

        if(TRUE){
            data.RMSE <-c(data.RMSE, round(sqrt(mean(FNN::knn.reg(dp, y=x, k=k)$residuals^2)),3))
            dwt.RMSE <- c(dwt.RMSE, round(sqrt(mean(FNN::knn.reg(dp.n, y=x, k=k)$residuals^2)),3))

            df1 <- data.frame(Group=1, s=s[i], No=1:sample.cal,Pred=FNN::knn.reg(dp, y=x, k=k)$pred, Obs=x)
            df2 <- data.frame(Group=2, s=s[i], No=1:sample.cal,Pred=FNN::knn.reg(dp.n, y=x, k=k)$pred, Obs=x)
        }

        df <- rbind(df, rbind(df1,df2))

    }
    summary(df)
    #------------------------------------------------------------------------------
    #prediction - cali
    ###bar plot: reduced RMSE
    #print(rbind(data.RMSE,dwt.RMSE))
    #expect_gt(data.RMSE,dwt.RMSE)

    dwt.RMSE.cal <- rbind(data.RMSE,dwt.RMSE)

    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------
    #validataion
    sample.val <- sample-sample.cal
    n <- sample.val
    J <- ceiling(log(n/(2*v-1))/log(2))
    if(wf=="haar"&&mode=="MODWT") J = J-1 #since modwt no need a dyadic number size
    print(paste0("Validation: Decomposition Levels J= ",J))

    #data for training
    dp <- lapply(dwt.list, function(ts) ts$dp)
    dp.n <- lapply(dwt.list, function(ts) ts$dp.n)
    x <- lapply(dwt.list, function(ts) ts$x)

    #new input data
    data.list <- lapply(ts.list, function(ts) list(x=ts$z[(sample.cal+1):sample],
                                                   dp=cbind(ts$x[(sample.cal+1):sample],
                                                            ts$y[(sample.cal+1):sample])))

    #------------------------------------------------------------------------------
    #variance transform
    if(mode=="MRA"){
        dwt.list.val<- lapply(1:length(data.list), function(i) dwt.vt.val(data.list[[i]], J, dwt.list[[i]]))

    } else if(mode=="MODWT"){
        dwt.list.val<- lapply(1:length(data.list), function(i) modwt.vt.val(data.list[[i]], J, dwt.list[[i]]))

    } else {
        dwt.list.val<- lapply(1:length(data.list), function(i) at.vt.val(data.list[[i]], dwt.list[[i]]))

    }

    #-------------------------------------------------------------------------------
    df.val <- NULL;data.RMSE <-NULL;dwt.RMSE<-NULL
    for(i in 1:length(dwt.list.val)){

        dwt <- dwt.list.val[[i]]
        dp.v <- dwt$dp; dp.n.v <- dwt$dp.n;

        if(TRUE){
            m1 <- FNN::knn.reg(train=dp[[i]], test=dp.v, y=unlist(x[[i]]), k=k)
            m2 <- FNN::knn.reg(train=dp.n[[i]], test=dp.n.v, y=unlist(x[[i]]), k=k)

            data.RMSE <-c(data.RMSE, round(sqrt(mean((m1$pred-dwt$x)^2)),3))
            dwt.RMSE <- c(dwt.RMSE, round(sqrt(mean((m2$pred-dwt$x)^2)),3))

            df1 <- data.frame(Group=1, s=s[i], No=1:sample.val,Pred=m1$pred, Obs=dwt$x)
            df2 <- data.frame(Group=2, s=s[i], No=1:sample.val,Pred=m2$pred, Obs=dwt$x)

        }

        df.val <- rbind(df.val, rbind(df1,df2))

    }
    summary(df.val)

    #------------------------------------------------------------------------------
    #prediction - vali
    ###bar plot: reduced RMSE
    #print(rbind(data.RMSE,dwt.RMSE))
    #expect_gt(data.RMSE,dwt.RMSE)

    dwt.RMSE.val <- rbind(data.RMSE,dwt.RMSE)

}

#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
###standard method
#calibration
n <- sample.cal
J <- floor(log10(n)) # (Nourani et al., 2008)

#--------------------------------------------------------------------------------------------
if(TRUE){
    # form new response and predictors dataset with same generated data ts.list
    data.list <- list()
    for(i in 1:length(ts.list)){
        #i <- 1
        x <- ts.list[[i]]$x[1:sample.cal]
        y <- ts.list[[i]]$y[1:sample.cal]
        z <- ts.list[[i]]$z[1:sample.cal]

        xx <- padding(x, pad); yy <- padding(y, pad)

        if(mode=="MRA"){
            mra.x <- matrix(unlist(lapply(mra(xx,wf,J,method,boundary), function(z) z[1:n])), ncol=J+1)
            mra.y <- matrix(unlist(lapply(mra(yy,wf,J,method,boundary), function(z) z[1:n])), ncol=J+1)

            # sum(abs(x-rowSums(mra.x)))
            # sum(abs(y-rowSums(mra.y)))

            data.list[[i]] <- list(x=as.numeric(z), dp=cbind(mra.x, mra.y))
        } else if(mode=="MODWT"){
            modwt.x <- modwt(x,wf,J,boundary)
            modwt.y <- modwt(y,wf,J,boundary)

            data.list[[i]] <- list(x=as.numeric(z), dp=cbind(modwt.x, modwt.y))
        } else {

            at.x <- at.wd(xx, v=2)
            at.y <- at.wd(yy, v=2)

            # sum(abs(x-rowSums(at.x[1:n,])))
            # sum(abs(y-rowSums(at.y[1:n,])))

            data.list[[i]] <- list(x=as.numeric(z), dp=cbind(at.x, at.y))

        }
    }

    #-------------------------------------------------------------------------------
    #prediction
    df <- NULL;data.RMSE<-NULL;dwt.RMSE<-NULL
    for(i in 1:length(data.list)){
        # i <- 1
        dwt <- data.list[[i]]
        x <- dwt$x; dp <- dwt$dp; #dp.n <- dwt$dp.n

        if(TRUE){
            data.RMSE <-c(data.RMSE, round(sqrt(mean(FNN::knn.reg(dp, y=x, k=k)$residuals^2)),3))
            #dwt.RMSE <- c(dwt.RMSE, round(sqrt(mean(FNN::knn.reg(dp.n, y=x, k=k)$residuals^2)),3))

            df1 <- data.frame(Group=1, s=s[i], No=1:sample.cal,Pred=FNN::knn.reg(dp, y=x, k=k)$pred, Obs=x)
            #df2 <- data.frame(Group=2, s=s[i], No=1:sample.cal,Pred=FNN::knn.reg(dp.n, y=x, k=k)$pred, Obs=x)
        }

        #df <- rbind(df, rbind(df1,df2))

    }
    #summary(df)
    #------------------------------------------------------------------------------
    #prediction - cali
    ###bar plot: reduced RMSE
    #print(rbind(data.RMSE,dwt.RMSE))
    #expect_gt(data.RMSE,dwt.RMSE)

    data.RMSE.Std.cal <- data.RMSE
    #dwt.RMSE.Std.cal <- dwt.RMSE

    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------
    #-------------------------------------------------------------------------------
    #validataion
    n <- sample.val
    J <- floor(log10(n)) # (Nourani et al., 2008)

    # form new response and predictors dataset
    data.list.val <- list()
    for(i in 1:length(ts.list)){
        #i <- 1
        x <- ts.list[[i]]$x[(sample.cal+1):sample]
        y <- ts.list[[i]]$y[(sample.cal+1):sample]
        z <- ts.list[[i]]$z[(sample.cal+1):sample]

        xx <- padding(x, pad); yy <- padding(y, pad)

        if(mode=="MRA"){
            mra.x <- matrix(unlist(lapply(mra(xx,wf,J,method,boundary), function(z) z[1:n])), ncol=J+1, byrow=FALSE)
            mra.y <- matrix(unlist(lapply(mra(yy,wf,J,method,boundary), function(z) z[1:n])), ncol=J+1, byrow=FALSE)

            data.list.val[[i]] <- list(x=as.numeric(z), dp=cbind(mra.x, mra.y))
        } else if(mode=="MODWT"){
            modwt.x <- modwt(x,wf,J,boundary)
            modwt.y <- modwt(y,wf,J,boundary)

            data.list.val[[i]] <- list(x=as.numeric(z), dp=cbind(modwt.x, modwt.y))
        } else {

            at.x <- at.wd(xx, v=2)
            at.y <- at.wd(yy, v=2)

            # sum(abs(x-rowSums(at.x[1:n,])))
            # sum(abs(y-rowSums(at.y[1:n,])))

            data.list.val[[i]] <- list(x=as.numeric(z), dp=cbind(at.x, at.y))
        }

    }

    #data for training
    dp <- lapply(data.list, function(ts) ts$dp)
    x <- lapply(data.list, function(ts) ts$x)

    #-------------------------------------------------------------------------------
    #prediction
    df.val <- NULL;data.RMSE <-NULL;dwt.RMSE<-NULL
    for(i in 1:length(data.list.val)){

        dwt <- data.list.val[[i]]
        dp.v <- dwt$dp; #dp.n.v <- dwt$dp.n;


        if(TRUE){
            m1 <- FNN::knn.reg(train=dp[[i]], test=dp.v, y=unlist(x[[i]]), k=k)
            #m2 <- FNN::knn.reg(train=dp.n, test=dp.n.v, y=unlist(x[[i]]), k=k)

            data.RMSE <-c(data.RMSE, round(sqrt(mean((m1$pred-dwt$x)^2)),3))
            #dwt.RMSE <- c(dwt.RMSE, round(sqrt(mean((m2$pred-dwt$x)^2)),3))

            df1 <- data.frame(Group=1, s=s[i], No=1:sample.val,Pred=m1$pred, Obs=dwt$x)
            #df2 <- data.frame(Group=2, s=s[i], No=1:sample.val,Pred=m2$pred, Obs=dwt$x)
        }

        #df.val <- rbind(df.val, rbind(df1,df2))

    }
    #summary(df.val)
    #------------------------------------------------------------------------------
    #prediction - vali
    ###bar plot: reduced RMSE
    #print(rbind(data.RMSE,dwt.RMSE))
    #expect_gt(data.RMSE,dwt.RMSE)

    data.RMSE.Std.val <- data.RMSE
    #dwt.RMSE.Std.val <- dwt.RMSE

    data.RMSE.Std <- rbind(data.RMSE.Std.cal,data.RMSE.Std.val)
    #dwt.RMSE.Std <- rbind(dwt.RMSE.Std.cal,dwt.RMSE.Std.val)
    print(data.RMSE.Std)
}


#------------------------------------------------------------------------------
#comparison
df.RMSE <- rbind(rbind(dwt.RMSE.cal,dwt.RMSE.val),data.RMSE.Std)
rownames(df.RMSE) <- NULL
df.RMSE.n <- data.frame(df.RMSE,Group=c("cal", "cal", "val", "val", "cal", "val"),
                        Model = c("Ori", "VT", "Ori", "VT", "Std", "Std"))

df.RMSE.avg <- aggregate(.~Group+Model,df.RMSE.n,mean)

t1 <- gather(df.RMSE.avg, s, Value,3:ncol(df.RMSE.avg)) %>% spread(Group, Value)
t1.n <- t1[order(t1$s),c("s","Model","cal","val")]

print(paste0(mode,"----------",method))
print(t1.n)

