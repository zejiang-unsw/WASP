#--------------------------------------------------------------------------------
#' Variance Transformation Operation - MRA
#' @param data		  A list of response x and dependent variables dp.
#' @param wf		    Name of the wavelet filter to use in the decomposition.
#' @param J      	  Specifies the depth of the decomposition. This must be a number less than or equal to log(length(x),2).
#' @param method	  Either "dwt" or "modwt".
#' @param pad		    The method used for extend data to dyadic size. Use "per", "zero", or "sym".
#' @param boundary  Character string specifying the boundary condition. If boundary=="periodic" the default, then the vector you decompose is assumed to be periodic on its defined interval, if boundary=="reflection", the vector beyond its boundaries is assumed to be a symmetric reflection of itself.
#' @param cov.opt   Options of Covariance matrix sign. Use "pos", "neg", or "auto".
#' @param flag      Biased or Unbiased variance transformation, c("biased","unbiased").
#' @param detrend   Detrend the input time series or just center, default (F)
#'
#' @return A list of 8 elements: wf, method, boundary, pad, x (data), dp (data), dp.n (variance trasnformed dp), and S (covariance matrix).
#' @import waveslim
#' @export
#'
#' @references Jiang, Z., Sharma, A., & Johnson, F. (2020). Refining Predictor Spectral Representation Using Wavelet Theory for Improved Natural System Modeling. Water Resources Research, 56(3), e2019WR026962.
#'
#' @examples
#' data(rain.mon)
#' data(obs.mon)
#'
#' ## response SPI - calibration
#' # SPI.cal <- SPI.calc(window(rain.mon, start=c(1949,1), end=c(1979,12)),sc=12)
#' SPI.cal <- SPEI::spi(window(rain.mon, start = c(1949, 1), end = c(1979, 12)), scale = 12)$fitted
#'
#' ## create paired response and predictors dataset for each station
#' data.list <- list()
#' for (id in seq_len(ncol(SPI.cal))) {
#'   x <- window(SPI.cal[, id], start = c(1950, 1), end = c(1979, 12))
#'   dp <- window(obs.mon, start = c(1950, 1), end = c(1979, 12))
#'   data.list[[id]] <- list(x = as.numeric(x), dp = matrix(dp, nrow = nrow(dp)))
#' }
#'
#' ## variance transformation
#' dwt.list <- lapply(data.list, function(x) {
#'   dwt.vt(x, wf = "d4", J = 7, method = "dwt", pad = "zero", boundary = "periodic", cov.opt = "auto")
#' })
#'
#' ## plot original and reconstrcuted predictors for each station
#' for (i in seq_len(length(dwt.list))) {
#'   # extract data
#'   dwt <- dwt.list[[i]]
#'   x <- dwt$x # response
#'   dp <- dwt$dp # original predictors
#'   dp.n <- dwt$dp.n # variance transformed predictors
#'
#'   plot.ts(cbind(x, dp))
#'   plot.ts(cbind(x, dp.n))
#' }
dwt.vt <- function(data, wf, J, method, pad, boundary, cov.opt = "auto",
                   flag = "biased", detrend = FALSE) {
  # initialization
  x <- data$x
  dp <- as.matrix(data$dp)
  mu.dp <- apply(dp, 2, mean)

  # variance transfrom
  ndim <- ncol(dp)
  n <- nrow(dp)
  S <- matrix(nrow = J + 1, ncol = ndim)
  dp.n <- matrix(nrow = n, ncol = ndim)
  idwt.dp <- vector("list", ndim)

  for (i in 1:ndim) {
    # center or detrend
    if (!detrend) {
      dp.c <- scale(dp[, i], scale = F)
    } else {
      # dp.c <- lm(dp[,i]~c(1:n))$resid
      dp.c <- dp[, i] - smooth.spline(1:n, dp[, i], spar = 1)$y
    }
    dp.p <- padding(dp.c, pad = pad)

    # Multiresolution Analysis
    idwt.dp[[i]] <- waveslim::mra(dp.p, wf = wf, J = J, method = method, boundary = boundary)
    B <- matrix(unlist(lapply(idwt.dp[[i]], function(z) z[1:n])), ncol = J + 1, byrow = FALSE)

    Bn <- scale(B)
    V <- as.numeric(apply(B, 2, sd))

    dif <- sum(abs(Bn %*% V - dp.c))
    if (dif > 10^-10) print(paste0("Difference between reconstructed and
                                   original series: ", dif))

    # variance transformation
    cov <- cov(x, Bn[seq_len(length(x)), ])
    # cat("Biased cov: ", round(cov,3),"\n")
    # cat("Biased cov1: ", round(cor(x, B[seq_len(length(x)),])*sd(x),3),"\n")

    if (flag == "unbiased") { ### unbiased wavelet variance - only change cov
      idwt.dp.n <- non.bdy(idwt.dp[[i]], wf = wf, method = "mra")

      B.n <- matrix(unlist(lapply(idwt.dp.n, function(z) z[1:n])), ncol = J + 1, byrow = FALSE)
      cov <- cov(x, scale(B.n)[seq_len(length(x)), ], use = "pairwise.complete.obs")

      # if unbiased cov is not available then use biased
      cov[is.na(cov)] <- cov(x, Bn[seq_len(length(x)), ])[is.na(cov)]
      # cat("Unbiased: ",round(cov,3),"\n")
    }

    if (cov.opt == "pos") cov <- cov else if (cov.opt == "neg") cov <- -cov
    S[, i] <- as.vector(cov)

    Vr <- as.numeric(cov / norm(cov, type = "2") * sd(dp.c))
    # cat(norm(cov,type="2"),"\n")
    # cat(norm(cor(x, B[seq_len(length(x)),]),type="2"),"\n")
    # cat("Biased alpha: ", round(Vr,3),"\n")
    # cat("Biased alpha1: ", round(cor(x, B[seq_len(length(x)),])/norm(cor(x, B[seq_len(length(x)),])/sd(dp.c),type="2"),3),"\n")

    if (!detrend) {
      dp.n[, i] <- Bn %*% Vr + mu.dp[i]
    } else {
      # dp.n[,i] <- Bn%*%Vr + lm(dp[,i]~c(1:n))$fitted
      dp.n[, i] <- Bn %*% Vr + smooth.spline(1:n, dp[, i], spar = 1)$y
    }

    # check the correlation after vt then decide the direction of C
    if (cov.opt == "auto") {
      if (cor.test(dp.n[, i], dp[, i])$estimate < 0 & cor.test(dp.n[, i], dp[, i])$p.value < 0.05) cov <- -cov
      S[, i] <- as.vector(cov)

      Vr <- as.numeric(cov / norm(cov, type = "2") * sd(dp.c))
      if (!detrend) {
        dp.n[, i] <- Bn %*% Vr + mu.dp[i]
      } else {
        # dp.n[,i] <- Bn%*%Vr + lm(dp[,i]~c(1:n))$fitted
        dp.n[, i] <- Bn %*% Vr + smooth.spline(1:n, dp[, i], spar = 1)$y
      }
    }

    dif.var <- abs(var(dp[, i]) - var(dp.n[, i])) / var(dp[, i])
    if (dif.var > 0.15) print(paste0("Variance difference between transformed
                        and original series by percentage: ", dif.var * 100))
  }

  dwt <- list(
    wavelet = wf,
    method = method,
    boundary = boundary,
    pad = pad,
    x = x,
    dp = dp,
    dp.n = dp.n,
    S = S
  )
  class(dwt) <- "dwt"

  return(dwt)
}

#--------------------------------------------------------------------------------
#' Variance Transformation Operation for Validation
#' @param data		  A list of response x and dependent variables dp.
#' @param J      	  Specifies the depth of the decomposition. This must be a number less than or equal to log(length(x),2).
#' @param dwt       A class of "dwt" data. Output from dwt.vt().
#' @param detrend   Detrend the input time series or just center, default (F)
#'
#' @return          A list of 8 elements: wf, method, boundary, pad, x (data), dp (data), dp.n (variance trasnformed dp), and S (covariance matrix).
#' @export
#' @references Jiang, Z., Sharma, A., & Johnson, F. (2020). Refining Predictor Spectral Representation Using Wavelet Theory for Improved Natural System Modeling. Water Resources Research, 56(3), e2019WR026962. doi:10.1029/2019wr026962
#'
#' @examples
#' data(rain.mon)
#' data(obs.mon)
#'
#' ## response SPI - calibration
#' # SPI.cal <- SPI.calc(window(rain.mon, start=c(1949,1), end=c(1979,12)),sc=12)
#' SPI.cal <- SPEI::spi(window(rain.mon, start = c(1949, 1), end = c(1979, 12)), scale = 12)$fitted
#'
#' ## create paired response and predictors dataset for each station
#' data.list <- list()
#' for (id in seq_len(ncol(SPI.cal))) {
#'   x <- window(SPI.cal[, id], start = c(1950, 1), end = c(1979, 12))
#'   dp <- window(obs.mon, start = c(1950, 1), end = c(1979, 12))
#'   data.list[[id]] <- list(x = as.numeric(x), dp = matrix(dp, nrow = nrow(dp)))
#' }
#'
#' ## variance transformation - calibration
#' dwt.list <- lapply(data.list, function(x) {
#'   dwt.vt(x, wf = "d4", J = 7, method = "dwt", pad = "zero", boundary = "periodic", cov.opt = "auto")
#' })
#'
#' ## response SPI - validation
#' # SPI.val <- SPI.calc(window(rain.mon, start=c(1979,1), end=c(2009,12)),sc=12)
#' SPI.val <- SPEI::spi(window(rain.mon, start = c(1979, 1), end = c(2009, 12)), scale = 12)$fitted
#'
#' ## create paired response and predictors dataset for each station
#' data.list <- list()
#' for (id in seq_len(ncol(SPI.val))) {
#'   x <- window(SPI.val[, id], start = c(1980, 1), end = c(2009, 12))
#'   dp <- window(obs.mon, start = c(1980, 1), end = c(2009, 12))
#'   data.list[[id]] <- list(x = as.numeric(x), dp = matrix(dp, nrow = nrow(dp)))
#' }
#'
#' # variance transformation - validation
#' dwt.list.val <- lapply(
#'   seq_len(length(data.list)),
#'   function(i) dwt.vt.val(data.list[[i]], J = 7, dwt.list[[i]])
#' )
#'
#' ## plot original and reconstrcuted predictors for each station
#' for (i in seq_len(length(dwt.list.val))) {
#'   # extract data
#'   dwt <- dwt.list.val[[i]]
#'   x <- dwt$x # response
#'   dp <- dwt$dp # original predictors
#'   dp.n <- dwt$dp.n # variance transformed predictors
#'
#'   plot.ts(cbind(x, dp))
#'   plot.ts(cbind(x, dp.n))
#' }
dwt.vt.val <- function(data, J, dwt, detrend = FALSE) {

  # initialization
  x <- data$x
  dp <- as.matrix(data$dp)
  wf <- dwt$wavelet
  method <- dwt$method
  boundary <- dwt$boundary
  pad <- dwt$pad
  mu.dp <- apply(dp, 2, mean)

  # variance transfrom
  ndim <- ncol(dp)
  n <- nrow(dp)
  dp.n <- matrix(nrow = n, ncol = ndim)
  idwt.dp <- vector("list", ndim)

  for (i in 1:ndim) {
    # center or detrend
    if (!detrend) {
      dp.c <- scale(dp[, i], scale = FALSE)
    } else {
      # dp.c <- lm(dp[,i]~c(1:n))$resid
      dp.c <- dp[, i] - smooth.spline(1:n, dp[, i], spar = 1)$y
    }
    dp.p <- padding(dp.c, pad = pad)

    # Multiresolution Analysis
    idwt.dp[[i]] <- waveslim::mra(dp.p, wf = wf, J = J, method = method, boundary = boundary)
    B <- matrix(unlist(lapply(idwt.dp[[i]], function(z) z[1:n])), ncol = J + 1, byrow = FALSE)

    Bn <- scale(B)
    V <- as.numeric(apply(B, 2, sd))

    dif <- sum(abs(Bn %*% V - dp.c))
    if (dif > 10^-10) print(paste0("Difference between reconstructed
                                   and original series: ", dif))

    # in case different J
    cov <- rep(0, J + 1)
    if (length(dwt$S[, i]) > (J + 1)) {
      cov <- dwt$S[, i][1:(J + 1)]
    } else {
      cov[seq_len(length(dwt$S[, i]))] <- dwt$S[, i]
    }

    Vr <- as.numeric(cov / norm(cov, type = "2") * sd(dp.c))

    if (!detrend) {
      dp.n[, i] <- Bn %*% Vr + mu.dp[i]
    } else {
      # dp.n[,i] <- Bn%*%Vr + lm(dp[,i]~c(1:n))$fitted
      dp.n[, i] <- Bn %*% Vr + smooth.spline(1:n, dp[, i], spar = 1)$y
    }


    # dif.var <- abs(var(dp[,i])-var(dp.n[,i]))/var(dp[,i])
    # if (dif.var > 0.15) print(paste0("Variance difference between transformed
    #                     and original series by percentage: ", dif.var * 100))
  }

  dwt <- list(
    wavelet = wf,
    method = method,
    boundary = boundary,
    pad = pad,
    x = x,
    dp = dp,
    dp.n = dp.n,
    S = dwt$S
  )
  class(dwt) <- "dwt"

  return(dwt)
}

#-------------------------------------------------------------------------------
#' Padding data to dyadic sample size
#' @param x     A vector or time series containing the data be to decomposed.
#' @param pad   Method for padding, including periodic, zero and symetric padding.
#'
#' @return      A dyadic length (power of 2) vector or time series.
#' @export
#'
#' @examples
#' x <- rnorm(360)
#' x1 <- padding(x, pad = "per")
#' x2 <- padding(x, pad = "zero")
#' x3 <- padding(x, pad = "sym")
#' ts.plot(cbind(x, x1, x2, x3), col = 1:4)
padding <- function(x, pad = c("per", "zero", "sym")) {
  n <- length(x)
  N <- 2^(ceiling(log(n, 2)))
  if (pad == "per") {
    xx <- c(x, x)[1:N]
  } else if (pad == "zero") {
    xx <- c(x, rep(0, N - n))
  } else {
    xx <- c(x, rev(x))[1:N]
  }
}
