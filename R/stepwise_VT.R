#' Calculate stepwise high order VT in calibration
#'
#' @param data    A list of data, including response and predictors
#' @param alpha   The significance level used to judge whether the sample estimate is significant. A default alpha value is 0.1.
#' @param nvarmax The maximum number of variables to be selected.
#' @param mode    A mode of variance transformation, i.e., MRA, MODWT, or AT
#' @param wf      Wavelet family
#' @param J      	Specifies the depth of the decomposition. This must be a number less than or equal to log(length(x),2).
#' @param method  Either "dwt" or "modwt" of MRA.
#' @param pad      The method used for extend data to dyadic size. Use "per", "zero", or "sym".
#' @param boundary Character string specifying the boundary condition. If boundary=="periodic" the default, then the vector you decompose is assumed to be periodic on its defined interval, if boundary=="reflection", the vector beyond its boundaries is assumed to be a symmetric reflection of itself.
#' @param cov.opt Options of Covariance matrix sign. Use "pos", "neg", or "auto".
#' @param flag    Biased or Unbiased variance transformation.
#' @param detrend Detrend the input time series or just center, default (F).
#'
#'
#' @return A list of 2 elements: the column numbers of the meaningful predictors (cpy), and partial informational correlation (cpyPIC).
#' @export
#'
#' @references Sharma, A., Mehrotra, R., 2014. An information theoretic alternative to model a natural system using observational information alone. Water Resources Research, 50(1): 650-660.
#'
#' Jiang, Z., Sharma, A., & Johnson, F. (2021). Variable transformations in the spectral domain – Implications for hydrologic forecasting. Journal of Hydrology, 126816.
#'
#' @examples
#' ### Real-world example
#' data("rain.mon")
#' data("obs.mon")
#' mode <- switch(1,
#'   "MRA",
#'   "MODWT",
#'   "AT"
#' )
#' wf <- "d4"
#' station.id <- 5 # station to investigate
#' SPI.12 <- SPEI::spi(rain.mon, scale = 12)$fitted
#' lab.names <- colnames(obs.mon)
#' # plot.ts(SPI.12[,1:10])
#'
#' x <- window(SPI.12[, station.id], start = c(1950, 1), end = c(1979, 12))
#' dp <- window(obs.mon[, lab.names], start = c(1950, 1), end = c(1979, 12))
#'
#' data <- list(x = x, dp = matrix(dp, ncol = ncol(dp)))
#'
#' dwt <- stepwise.VT(data, mode = mode, wf = wf, flag = "biased")
#'
#' ### plot transformed predictor before and after
#' cpy <- dwt$cpy
#' op <- par(mfrow = c(length(cpy), 1), mar = c(2, 3, 2, 1))
#' for (i in seq_along(cpy)) {
#'   ts.plot(cbind(dwt$dp[, i], dwt$dp.n[, i]), xlab = "NA", col = 1:2)
#' }
#' par(op)
stepwise.VT <- function(data, alpha = 0.1, nvarmax = 4, mode = c("MRA", "MODWT", "AT"), wf, J, method = "dwt", pad = "zero",
                        boundary = "periodic", cov.opt = "auto", flag = "biased", detrend = FALSE) {
  x <- as.matrix(data$x)
  py <- as.matrix(data$dp)

  n <- nrow(x)
  npy <- ncol(py)
  cpy <- cpyPIC <- NULL
  icpy <- 0
  z <- NULL
  z.vt <- NULL
  S <- NULL
  r2 <- NULL
  isig <- T
  icoloutz <- 1:npy

  # cat("calc.PIC-----------","\n")
  while (isig) {
    npicmax <- npy - icpy
    pictemp <- rep(0, npicmax)
    y <- py[, icoloutz]

    temp <- pic.calc(x, y, z, mode, wf, J, method, pad, boundary, cov.opt, flag, detrend)
    pictemp <- temp$pic

    pytmp <- temp$py
    Stmp <- temp$S

    ctmp <- order(-pictemp)[1]
    cpytmp <- icoloutz[ctmp]
    picmaxtmp <- pictemp[ctmp]

    cpy <- c(cpy, cpytmp)
    cpyPIC <- c(cpyPIC, picmaxtmp)
    z <- cbind(z, py[, cpytmp])

    S <- cbind(S, Stmp[, ctmp])
    z.vt <- cbind(z.vt, pytmp[, ctmp])

    z <- z.vt # mathematically this is more valid

    if (!is.null(z)) {
      z <- as.matrix(z)
      df <- n - ncol(z)
    } else {
      df <- n
    }
    t <- qt(1 - alpha, df = df)
    picthres <- sqrt(t^2 / (t^2 + df))
    # picthres <- qt((0.5 + alpha/2), df)

    # method 1 and 2
    u <- x - knnregl1cv(x, z.vt[1:n, ])
    # u <- lm.fit(z.vt, x)$residuals
    r2 <- c(r2, 1 - sum(u^2) / sum((x - mean(x))^2))

    # method 3
    # r2 <- c(r2, FNN::knn.reg(z.vt, y=x, k=ceiling(sqrt(length(x)/2)))$R2Pred)

    icoloutz <- icoloutz[-ctmp]
    icpy <- icpy + 1

    if (icpy > 1) {
      # r2thres <- r2.boot(z.vt, x, prob=1-alpha)
      # cat("r2thres: ",r2thres,"\n")

      # if(r2[icpy]<r2[icpy-1]|r2[icpy]<r2thres) {
      # cat(r2[icpy]<r2[icpy-1],"|",r2[icpy]<r2thres,"\n")
      if (r2[icpy] < r2[icpy - 1] | picmaxtmp < picthres) {
        # cat(r2[icpy]<r2[icpy-1],"|",picmaxtmp<picthres,"\n")
        # if(r2[icpy]<r2[icpy-1]) {
        isig <- F

        cpy <- cpy[-icpy]
        cpyPIC <- cpyPIC[-icpy]

        z <- as.matrix(z[, -icpy])
        z.vt <- as.matrix(z.vt[, -icpy])
        S <- as.matrix(S[, -icpy])
      }
    }
    if ((npy - icpy) == 0 | icpy >= nvarmax) isig <- F
  }
  # cat("R2: ",r2,"\n")
  # cat("calc.PW------------","\n")
  if (!is.null(z)) {
    out <- pw.calc(x, z.vt, cpyPIC)

    outwt <- out$pw
    lstwt <- abs(lsfit(z.vt[seq_along(x), ], x)$coef)

    z.n <- py[, cpy]
    ncpy <- length(cpy)
    if (ncpy > 1) {
      for (i in 2:ncpy) {
        tmp <- z[, 1:(i - 1)]
        z.n[, i] <- z.n[, i] - knnregl1cv(z.n[, i], tmp)
        # z.n[,i] <- lm.fit(as.matrix(tmp), z.n[,i])$residuals
      }
    }

    return(list(
      cpy = cpy, cpyPIC = cpyPIC, wt = outwt, lstwet = lstwt,
      x = x, py = py, r2 = r2,
      dp = z.n, dp.n = z.vt, S = S,
      wavelet = wf, method = method, pad = pad, boundary = boundary
    ))
  } else {
    message("None of the provided predictors is related to the response variable")
  }
}

#' Calculate stepwise high order VT in validation
#'
#' @param data    A list of data, including response and predictors
#' @param J      	Specifies the depth of the decomposition. This must be a number less than or equal to log(length(x),2).
#' @param dwt     Output from dwt.vt(), including the transformation covariance
#' @param mode    A mode of variance transformation, i.e., MRA, MODWT, or AT
#' @param detrend Detrend the input time series or just center, default (F)
#'
#' @return        A list of objects, including transformed predictors
#' @export
#'
#' @examples
#' ### Real-world example
#' data("rain.mon")
#' data("obs.mon")
#' mode <- switch(1,
#'   "MRA",
#'   "MODWT",
#'   "a trous"
#' )
#' wf <- "d4"
#' station.id <- 5 # station to investigate
#' SPI.12 <- SPEI::spi(rain.mon, scale = 12)$fitted
#' lab.names <- colnames(obs.mon)
#' # plot.ts(SPI.12[,1:10])
#'
#' #--------------------------------------
#' ### calibration
#' x <- window(SPI.12[, station.id], start = c(1950, 1), end = c(1979, 12))
#' dp <- window(obs.mon[, lab.names], start = c(1950, 1), end = c(1979, 12))
#'
#' data <- list(x = x, dp = matrix(dp, ncol = ncol(dp)))
#' dwt <- stepwise.VT(data, mode = mode, wf = wf, flag = "biased")
#' cpy <- dwt$cpy
#' #--------------------------------------
#' ### validation
#' x <- window(SPI.12[, station.id], start = c(1980, 1), end = c(2009, 12))
#' dp <- window(obs.mon[, lab.names], start = c(1980, 1), end = c(2009, 12))
#'
#' data.n <- list(x = x, dp = matrix(dp, ncol = ncol(dp)))
#' dwt.val <- stepwise.VT.val(data = data.n, dwt = dwt, mode = mode)
#'
#' ### plot transformed predictor before and after
#' op <- par(mfrow = c(length(cpy), 1), mar = c(0, 3, 2, 1))
#' for (i in seq_along(cpy))
#' {
#'   ts.plot(cbind(dwt.val$dp[, i], dwt.val$dp.n[, i]), xlab = "NA", col = 1:2)
#' }
#' par(op)
stepwise.VT.val <- function(data, J, dwt, mode = c("MRA", "MODWT", "AT"), detrend = FALSE) {
  # initialization
  x <- data$x
  py <- as.matrix(data$dp)
  cpy <- dwt$cpy
  ncpy <- length(cpy)
  wf <- dwt$wavelet
  boundary <- dwt$boundary
  if (mode == "MRA") {
    method <- dwt$method
    pad <- dwt$pad
  }

  if (wf != "haar") v <- as.integer(readr::parse_number(wf) / 2) else v <- 1
  # Maximum decomposition level J
  n <- length(x)
  # if(wf=="haar") J <- ceiling(log(n/(2*v-1))/log(2))-1 else J <- ceiling(log(n/(2*v-1))/log(2))#(Kaiser, 1994)
  if (missing(J)) J <- ceiling(log(n / (2 * v - 1)) / log(2)) - 1

  dwt.n <- c(dwt, method = method, boundary = boundary, pad = pad)
  dp <- dp.n <- as.matrix(py[, cpy])

  for (i in 1:ncpy) {
    data.n <- list(x = x, dp = as.matrix(dp[, 1:i]))

    # variance transform
    if (mode == "MRA") {
      dwt.val <- dwt.vt.val(data.n, J, dwt.n, detrend)
    } else if (mode == "MODWT") {
      dwt.val <- modwt.vt.val(data.n, J, dwt.n, detrend)
    } else {
      dwt.val <- at.vt.val(data.n, J, dwt.n, detrend)
    }

    dp.n[, 1:i] <- dwt.val$dp.n

    if ((i + 1) <= ncpy) {
      dp[, i + 1] <- dp[, i + 1] - knnregl1cv(dp[, i + 1], dp.n[, 1:i])
      # dp[,i+1] <- lm.fit(as.matrix(dp.n[,1:i]), dp[,i+1])$residuals
    } else {
      break
    }
  }

  return(c(dwt.val, py = list(py)))
}

# Calculate the ratio of conditional error standard deviations
#
# @param x     A vector of response.
# @param zin   A matrix containing the meaningful predictors selected from a large set of possible predictors (z).
# @param zout  A matrix containing the remaining possible predictors after taking out the meaningful predictors (zin).
#
# @return The STD ratio.
#
# @references Sharma, A., Mehrotra, R., 2014. An information theoretic alternative to model a natural system using observational information alone. Water Resources Research, 50(1): 650-660.
calc.scaleSTDratio <- function(x, zin, zout) {
  if (!missing(zout)) {
    zout <- as.matrix(zout)
    xhat <- knnregl1cv(x, zout[seq_along(x), ])
    stdratxzout <- sqrt(var(x - xhat) / var(x))
    zinhat <- knnregl1cv(zin, zout)
    stdratzinzout <- sqrt(var(zin - zinhat) / var(zin))
    return(0.5 * (stdratxzout + stdratzinzout))
  } else {
    return(1)
  }
}
#-------------------------------------------------------------------------------
#' R2 threshold by re-sampling approach
#'
#' @param z.vt  Identified independent variables
#' @param x     Response or dependent variable
#' @param prob  Probability with values in [0,1].
#'
#' @return A quantile assosciated with prob.
#' @export
#'
r2.boot <- function(z.vt, x, prob) {
  z.boot <- vapply(
    1:100, function(i) sample(z.vt[, ncol(z.vt)], replace = FALSE),
    numeric(nrow(z.vt))
  )

  # method 1 and 2
  # u.boot <- apply(z.boot, 2, function(i) x-knnregl1cv(x, cbind(z.vt[,-ncol(z.vt)],i)))
  u.boot <- apply(z.boot, 2, function(i) lm.fit(cbind(z.vt[, -ncol(z.vt)], i), x)$residuals)
  r2.boot <- apply(u.boot, 2, function(i) 1 - sum(i^2) / sum((x - mean(x))^2))

  # method 3
  # r2.boot <- apply(z.boot, 2, function(i) FNN::knn.reg(cbind(z.vt[,-ncol(z.vt)],i), y=x, k=ceiling(sqrt(length(x)/2)))$R2Pred)

  r2thres <- quantile(r2.boot, probs = prob, type = 8)

  return(r2thres)
}
#-------------------------------------------------------------------------------
kernel.est.uvn <- function(Z) {
  N <- length(Z)
  d <- 1
  # compute sigma & constant
  sigma <- 1.5 * bw.nrd0(Z)
  # sigma <- bw.nrd(Z)
  constant <- sqrt(2 * pi) * sigma * N

  # Commence main loop
  dens <- vector()
  for (h in 1:N) {
    dis.Z <- (Z - Z[h])^2
    exp.dis <- exp(-dis.Z / (2 * sigma^2))
    dens[h] <- sum(exp.dis) / constant
  }
  return(dens)
}
#-------------------------------------------------------------------------------
kernel.est.mvn <- function(Z) {

  # Compute covariance and determinant of cov
  N <- nrow(Z)
  d <- ncol(Z)
  Cov <- cov(Z)
  det.Cov <- det(Cov)

  # Compute sigma & constant
  sigma <- 1.5 * (4 / (d + 2))^(1 / (d + 4)) * N^(-1 / (d + 4))
  # sigma <- (4/(d + 2))^(1/(d + 4)) * N^(-1/(d + 4))
  constant <- (sqrt(2 * pi) * sigma)^d * sqrt(det.Cov) * N

  # Commence main loop
  dens <- vector()
  for (h in 1:N) {
    dist.val <- mahalanobis(Z, center = Z[h, ], cov = Cov)
    exp.dis <- exp(-dist.val / (2 * sigma^2))
    dens[h] <- sum(exp.dis) / constant
  }

  return(dens)
}
#-------------------------------------------------------------------------------
pmi.calc <- function(X, Y) {
  N <- length(X)

  pdf.X <- kernel.est.uvn(X)
  pdf.Y <- kernel.est.uvn(Y)
  pdf.XY <- kernel.est.mvn(cbind(X, Y))

  calc <- log(pdf.XY / (pdf.X * pdf.Y))
  return(sum(calc) / N)
}
#-------------------------------------------------------------------------------
#' Calculate PIC
#'
#' @param X       A vector of response.
#' @param Y       A matrix of new predictors.
#' @param Z       A matrix of pre-existing predictors that could be NULL if no prior predictors exist.
#' @param mode    A mode of variance transfomration, i.e., MRA, MODWT, or AT
#' @param J       The maximum decomposition level
#' @param wf      Wavelet family
#' @param method  Either "dwt" or "modwt" of MRA.
#' @param pad      The method used for extend data to dyadic size. Use "per", "zero", or "sym".
#' @param boundary Character string specifying the boundary condition. If boundary=="periodic" the default, then the vector you decompose is assumed to be periodic on its defined interval, if boundary=="reflection", the vector beyond its boundaries is assumed to be a symmetric reflection of itself.
#' @param cov.opt Options of Covariance matrix sign. Use "pos", "neg", or "auto".
#' @param flag    Biased or Unbiased variance transformation.
#' @param detrend Detrend the input time series or just center, default (F).
#'
#' @return A list of 2 elements: the partial mutual information (pmi), and partial informational correlation (pic).
#' @export
#'
#' @references Sharma, A., Mehrotra, R., 2014. An information theoretic alternative to model a natural system using observational information alone. Water Resources Research, 50(1): 650-660.
#' @references Galelli S., Humphrey G.B., Maier H.R., Castelletti A., Dandy G.C. and Gibbs M.S. (2014) An evaluation framework for input variable selection algorithms for environmental data-driven models, Environmental Modelling and Software, 62, 33-51, DOI: 10.1016/j.envsoft.2014.08.015.

pic.calc <- function(X, Y, Z, mode, wf, J, method = "dwt", pad = "zero",
                     boundary = "periodic", cov.opt = "auto", flag = "biased", detrend = F) {
  Y <- as.matrix(Y)

  if (wf != "haar") v <- as.integer(readr::parse_number(wf) / 2) else v <- 1
  # Maximum decomposition level J
  n <- length(X)
  # if(wf=="haar") J <- ceiling(log(n/(2*v-1))/log(2))-1 else J <- ceiling(log(n/(2*v-1))/log(2))#(Kaiser, 1994)
  if (missing(J)) J <- ceiling(log(n / (2 * v - 1)) / log(2)) - 1

  if (is.null(Z)) {
    x.in <- X
    y.in <- Y
  } else {
    Z <- as.matrix(Z)
    x.in <- X - knnregl1cv(X, Z[seq_along(X), ])
    y.in <- apply(Y, 2, function(i) i - knnregl1cv(i, Z))

    # x.in <- lm.fit(as.matrix(Z[seq_along(X),]), X)$residuals
    # y.in <- apply(Y, 2, function(i) lm.fit(Z, i)$residuals)
  }

  data.list <- list(x = x.in, dp = y.in)

  # variance transform
  if (mode == "MRA") {
    dwt.list <- dwt.vt(data.list, wf, J, method, pad, boundary, cov.opt, flag, detrend)
  } else if (mode == "MODWT") {
    dwt.list <- modwt.vt(data.list, wf, J, boundary, cov.opt, flag, detrend)
  } else {
    dwt.list <- at.vt(data.list, wf, J, boundary, cov.opt, flag, detrend)
  }

  y.in <- dwt.list$dp.n

  pmi <- apply(y.in[1:n, ], 2, function(i) pmi.calc(x.in, i))
  pmi[pmi < 0] <- 0
  pic <- sqrt(1 - exp(-2 * pmi))

  return(list(
    pic = as.numeric(pic),
    # x = dwt.list$x, dp = dwt.list$dp,
    py = dwt.list$dp.n, S = dwt.list$S
  ))
}
#-------------------------------------------------------------------------------
# Calculate Partial Weight
#
# @param x       A vector of response.
# @param z       A matrix containing identified predictors of x.
# @param cpyPIC  Partial informational correlation (cpyPIC).
#
# @return A vector of partial weights(pw) of the same length of z.
#
# @references Sharma, A., Mehrotra, R., 2014. An information theoretic alternative to model a natural system using observational information alone. Water Resources Research, 50(1): 650-660.
pw.calc <- function(x, z, cpyPIC) {
  wt <- NA
  Z <- as.matrix(z)
  if (ncol(Z) == 1) {
    wt <- calc.scaleSTDratio(x, Z) * cpyPIC
  } else {
    for (i in seq_len(ncol(Z))) wt[i] <- calc.scaleSTDratio(x, Z[, i], Z[, -i]) * cpyPIC[i]
  }
  return(list(pw = wt))
  #return(list(pw = wt / sum(wt)))
}
