#' Modified k-nearest neighbour conditional bootstrap or regression function estimation with extrapolation
#' @param x			  A vector of response.
#' @param z			  A matrix of existing predictors.
#' @param zout    A matrix of predictor values the response is to be estimated at.
#' @param k     	The number of nearest neighbours used. The default value is 0, indicating Lall and Sharma default is used.
#' @param pw   		A vector of partial weights of the same length of z.
#' @param reg		  A logical operator to inform whether a conditional expectation should be output or not nensemble, Used if reg=F and represents the number of realisations that are generated Value.
#' @param nensemble				  An integer the specifies the number of ensembles used. The default is 100.
#' @param tailcorrection	  A logical value, T (default) or F, that denotes whether a reduced value of k (number of nearest neighbours) should be used in the tails of any conditioning plane. Whether one is in the tails or not is determined based on the nearest neighbour response value.
#' @param tailprob				  A scalar that denotes the p-value of the cdf (on either extreme) the tailcorrection takes effect. The default value is 0.25.
#' @param tailfac				    A scalar that specifies the lowest fraction of the default k that can be used in the tails. Depending on the how extreme one is in the tails, the actual k decreases linearly from k (for a p-value greater than tailprob) to tailfac*k proportional to the actual p-value of the nearest neighbour response, divided by tailprob. The default value is 0.2.
#' @param extrap				    A logical value, T (default) or F, that denotes whether a kernel extraplation method is used to predict x.

#' @return A matrix of responses having same rows as zout if reg=T, or having nensemble columns is reg=F.
#' @export
#' @references Sharma, A., Tarboton, D.G. and Lall, U., 1997. Streamflow simulation: A nonparametric approach. Water resources research, 33(2), pp.291-308.
#' @references Sharma, A. and O'Neill, R., 2002. A nonparametric approach for representing interannual dependence in monthly streamflow sequences. Water resources research, 38(7), pp.5-1.
#' @examples
#' \donttest{
#' # AR9 model   x(i)=0.3*x(i-1)-0.6*x(i-4)-0.5*x(i-9)+eps
#' data.ar9 <- data.gen.ar9(500)
#' x <- data.ar9$x # response
#' z <- data.ar9$dp # possible predictors
#'
#' zout <- ts(data.gen.ar9(500, ndim = ncol(z))$dp) # new input
#'
#' xhat1 <- xhat2 <- x
#' xhat1 <- knn(x, z, zout, k = 5, reg = TRUE, extrap = FALSE) # without extrapolation
#' xhat2 <- knn(x, z, zout, k = 5, reg = TRUE, extrap = TRUE) # with extrapolation
#'
#' ts.plot(ts(x), ts(xhat1), ts(xhat2), col = c("black", "red", "blue"),
#' ylim = c(-5, 5), lwd = c(2, 2, 1))
#' }
knn <- function(x, z, zout, k = 0, pw, reg = TRUE, nensemble = 100, tailcorrection = TRUE,
                tailprob = 0.25, tailfac = 0.2, extrap = TRUE) {
  x <- as.matrix(x)
  n <- nrow(x)
  if (k == 0) {
    k <- floor(0.5 + 1 * (sqrt(n)))
  }
  z <- as.matrix(z)
  nz <- ncol(z)
  zout <- as.matrix(zout)
  nczout <- ncol(zout)
  if (nz > 1) {
    if (nczout != nz) {
      zout <- t(zout)
    }
  }
  nzout <- nrow(zout)
  kern <- 1 / (1:k) / sum(1 / (1:k))
  kerncdf <- cumsum(kern)
  if (tailcorrection) {
    mink <- ceiling(tailfac * k)
    empcdf <- rank(x) / (n + 1)
  }
  if (missing(pw)) {
    pw <- rep(1, nz)
  }
  if (extrap) {
    Sxz <- var(x, z)
    Scond <- as.matrix(Sxz)
  }
  if (reg) {
    xhat <- rep(0, nzout)
  } else {
    xhat <- matrix(0, nzout, nensemble)
    randnum <- array(runif(nzout * nensemble), c(nzout, nensemble, k))
    randnum1 <- sweep(randnum, 3, kerncdf)
    randnum1[randnum1 < 0] <- 1
  }
  sd <- sqrt(diag(var(z)))
  for (j in 1:nz) { # standardize
    z[, j] <- z[, j] / (sd[j] / pw[j])
    zout[, j] <- zout[, j] / (sd[j] / pw[j])
  }
  for (i in 1:nzout) {
    z_out <- matrix(rep(zout[i, ], nrow(z)), nrow = nrow(z), byrow = TRUE)
    dtmp <- z_out - z
    d <- apply(dtmp * dtmp, 1, sum)
    ord <- order(d)[1:k]
    if (tailcorrection) {
      pvalue <- empcdf[ord[1]]
      if (pvalue > 0.5) {
        pvalue <- 1 - pvalue
      }
      if (pvalue < tailprob) {
        ktmp <- ceiling(mink + (k - mink) * pvalue / tailprob)
        oldkern <- kern
        oldkerncdf <- kerncdf
        kern <- rep(0, k)
        kern[1:ktmp] <- 1 / (1:ktmp) / sum(1 / (1:ktmp))
        kerncdf <- cumsum(kern)
      }
    } else {
      pvalue <- 0
    } # give some initial value
    if (reg) {
      if (!extrap) xhat[i] <- sum(x[ord] * kern) else xhat[i] <- sum((x[ord] + dtmp[ord, ] %*% t(Scond)) * kern)
    }
    if (!reg) {
      if (tailcorrection) {
        rand1 <- sweep(randnum[i, , ], 2, kerncdf)
        rand1[rand1 < 0] <- 1
      } else {
        rand1 <- randnum1[i, , ]
      }
      ord2 <- apply(-rand1, 1, order, na.last = FALSE)
      ord3 <- ord2[k, ]
      ord3 <- ord3 + 1
      ord3[ord3 == (k + 1)] <- 1
      xhat[i, ] <- x[ord[ord3]]
    }
    if (tailcorrection & (pvalue < tailprob)) {
      kern <- oldkern
      kerncdf <- oldkerncdf
    }
  }
  return(xhat)
}


#' Leave one out cross validation.
#'
#' @param x A vector of response.
#' @param z A matrix of predictors.
#' @param k The number of nearest neighbours used. The default is 0, indicating Lall and Sharma default is used.
#' @param pw A vector of partial weights of the same length of z.
#'
#' @return A vector of L1CV estimates of the response.
#' @export
#'
#' @references Lall, U., Sharma, A., 1996. A Nearest Neighbor Bootstrap For Resampling Hydrologic Time Series. Water Resources Research, 32(3): 679-693.
#' @references Sharma, A., Mehrotra, R., 2014. An information theoretic alternative to model a natural system using observational information alone. Water Resources Research, 50(1): 650-660.
knnregl1cv <- function(x, z, k = 0, pw) {
  x <- as.matrix(x)
  n <- nrow(x)
  if (k == 0) {
    k <- floor(0.5 + 3 * (sqrt(n)))
  }
  z <- as.matrix(z)
  nz <- ncol(z)
  sd <- sqrt(diag(var(z)))
  if (missing(pw)) {
    pw <- rep(1, nz)
  }
  for (j in 1:nz) z[, j] <- z[, j] / (sd[j] / pw[j])
  d <- as.matrix(dist(z))
  ord1 <- apply(d, 2, order)
  ord <- ord1[2:(k + 1), ]
  kern <- 1 / (1:k) / sum(1 / (1:k))
  xhat <- rep(0, n)
  for (j in 1:k) xhat <- xhat + x[ord[j, ]] * kern[j]
  return(xhat)
}
