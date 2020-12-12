#' Replace Boundary Wavelet Coefficients with Missing Values (NA).
#'
#' @param x       DWT/MODWT/AT object
#' @param wf      Character string; name of wavelet filter
#' @param method  Either dwt or modwt or mra
#'
#' @return Same object as x only with some missing values (NA).
#' @export
#'
#' @references Cornish, C. R., Bretherton, C. S., & Percival, D. B. (2006). Maximal overlap wavelet statistical analysis with application to atmospheric turbulence. Boundary-Layer Meteorology, 119(2), 339-374.
#'
non.bdy <- function (x, wf, method = c("dwt","modwt","mra"))
{
  m <- wave.filter(wf)$length
  for (j in 1:(length(x) - 1)) {
    if (method == "dwt")
      n <- ceiling((m - 2) * (1 - 1/2^j))
    else n <- (2^j - 1) * (m - 1)
    n <- min(n, length(x[[j]]))
    x[[j]][1:n] <- NA
    if(method=="mra") x[[j]][(length(x[[j]])-n+1):length(x[[j]])] <- NA
  }
  x[[j + 1]][1:n] <- NA
  if(method=="mra") x[[j+1]][(length(x[[j+1]])-n+1):length(x[[j+1]])] <- NA
  return(x)
}

#' Produces an estimate of the multiscale variance along with approximate confidence intervals.
#'
#' @param x       DWT/MODWT/AT object
#' @param type    character string describing confidence interval calculation; valid methods are gaussian, eta1, eta2, eta3, nongaussian
#' @param p       (one minus the) two-sided p-value for the confidence interval
#'
#' @return        Dataframe with as many rows as levels in the wavelet transform object. The first column provides the point estimate for the wavelet variance followed by the lower and upper bounds from the confidence interval.
#' @export
#'
#' @references Percival, D. B. (1995) Biometrika, 82, No. 3, 619-631.
#'
wave.var <- function (x, type = "eta3", p = 0.025)
{
    ci.gaussian <- function(x, y, p) {
        find.first <- function(v) {
            na.length <- sum(is.na(v))
            v[na.length + 1]
        }
        x.acf <- lapply(x, FUN = my.acf)
        Aj <- unlist(lapply(x.acf, FUN = function(v) sum(v *
            v, na.rm = TRUE))) - unlist(lapply(x.acf, FUN = find.first))^2/2
        wv.var <- 2 * Aj/unlist(lapply(x, FUN = function(v) sum(!is.na(v))))
        return(data.frame(wavevar = y, lower = y - qnorm(1 -
            p) * sqrt(wv.var), upper = y + qnorm(1 - p) * sqrt(wv.var)))
    }
    ci.eta1 <- function(x, y, p) {
        return(0)
    }
    ci.eta2 <- function(x, y, p) {
        return(0)
    }
    ci.eta3 <- function(x, y, p) {
        x.length <- unlist(lapply(x, FUN = function(v) sum(!is.na(v))))
        eta3 <- pmax(x.length/2^(1:length(x)), 1)
        return(data.frame(wavevar = y, lower = eta3 * y/qchisq(1 -
            p, eta3), upper = eta3 * y/qchisq(p, eta3)))
    }
    ci.nongaussian <- function(x, y, p) {
        K <- 5
        J <- length(x)
        x.length <- unlist(lapply(x, FUN = function(v) sum(!is.na(v))))
        x.ss <- unlist(lapply(x, FUN = function(v) v[!is.na(v)]^2))
        mt.var <- numeric(J)
        for (j in 1:J) {
            x.dpss <- dpss.taper(x.length[j], K, 4)
            V <- apply(x.dpss, 2, sum)
            J <- apply(x.dpss * x.ss[[j]], 2, sum)
            mt.var[j] <- sum((J - y[j] * V)^2)/K/x.length[j]
        }
        return(data.frame(wavevar = y, lower = y - qnorm(1 -
            p) * sqrt(mt.var), upper = y + qnorm(1 - p) * sqrt(mt.var)))
    }

    x.ss <- unlist(lapply(x, FUN = function(v) sum(v * v, na.rm = TRUE)))
    x.length <- unlist(lapply(x, FUN = function(v) sum(!is.na(v))))
    y <- x.ss/x.length # use (x.length-1) will result in same as var()

    switch(type, gaussian = ci.gaussian(x, y, p), eta1 = ci.eta1(x,
        y, p), eta2 = ci.eta2(x, y, p), eta3 = ci.eta3(x, y,
        p), nongaussian = ci.nongaussian(x, y, p), stop("Invalid selection of \"type\" for the confidence interval"))
}
