#' Plot function: Plot original time series and decomposed frequency components
#'
#' @param y           Original time series (Y).
#' @param y.mra       Decomposed frequency components (d1,d2,..,aJ).
#' @param limits.x    x limit for plot.
#' @param limits.y    y limit for plot.
#' @param type        type of wavelet coefficients, details or approximations.
#' @param ps          integer; the point size of text (but not symbols).
#' @param ...         arguments for plot().
#'
#' @return A plot with original time series and decomposed frequency components.
#' @export
#'
#' @examples
#' ### synthetic example
#' # frequency, sampled from a given range
#' fd <- c(3, 5, 10, 15, 25, 30, 55, 70, 95)
#' data.SW3 <- data.gen.SW(nobs = 512, fp = c(15, 25, 30), fd = fd)
#'
#' x <- data.SW3$x
#' xx <- padding(x, pad = "zero")
#' ### wavelet transfrom
#' # wavelet family, extension mode and package
#' wf <- "d4" # wavelet family D8 or db4
#' boundary <- "periodic"
#' pad <- "zero"
#' if (wf != "haar") v <- as.integer(as.numeric(substr(wf, 2, 3)) / 2) else v <- 1
#'
#' # Maximum decomposition level J
#' n <- length(x)
#' J <- ceiling(log(n / (2 * v - 1)) / log(2)) # (Kaiser, 1994)
#'
#' ### decomposition
#' x.mra <- waveslim::mra(xx, wf = wf, J = J, method = "dwt", boundary = "periodic")
#' x.mra.m <- matrix(unlist(x.mra), ncol = J + 1)
#'
#' print(sum(abs(x - rowSums(x.mra.m[1:n, ])))) # additive check
#' var(x)
#' sum(apply(x.mra.m[1:n, ], 2, var)) # variance check
#'
#' limits.x <- c(0, n)
#' limits.y <- c(-3, 3)
#' mra.plot(x, x.mra.m, limits.x, limits.y, type = "details")
mra.plot <- function(y, y.mra, limits.x, limits.y, type = c("details", "coefs"), ps=12,...) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  if (type == "details") {
    ylab <- c("d", "a")
  } else {
    ylab <- c("D", "A")
  }

  J <- ncol(y.mra) - 1
  par(
    mfcol = c(J + 2, 1), mar = c(0, 3, 1, 1), # margin of the plot
    oma = c(2, 1, 1, 1), # move plot to the right and up
    mgp = c(1, 0.6, 0), # move axis labels closer to axis
    bg = "transparent", pty = "m", # maximal plotting region
    ps = ps # text size

  )

  plot(y,
    type = "l", axes = FALSE, xlab = NA, xaxs = "i",
    xlim = limits.x, ylim = limits.y, ...
  )
  box()

  for (i in 1:J) {
    plot(y.mra[, i],
      type = "l", axes = FALSE, xlab = NA, ylab = paste0(ylab[1], i), xaxs = "i",
      xlim = limits.x, ylim = limits.y
    )
    box()
  }

  plot(y.mra[, J + 1],
    type = "l", axes = FALSE, xlab = NA, ylab = paste0(ylab[2], J), xaxs = "i",
    xlim = limits.x, ylim = limits.y
  )
  box()
  axis(side = 1, at = seq(limits.x[1], limits.x[2], by = 50), labels = seq(limits.x[1], limits.x[2], by = 50))

  return(recordPlot())
}


#' Plot function: Variance structure before and after variance transformation
#'
#' @param dwt.data  Output data from variance transformation function
#'
#' @return A plot with variance structure before and after variance transformation.
#' @export
#' @import ggplot2
#'
#' @examples
#' data("data.HL")
#' data("data.SW1")
#'
#' # variance transfrom
#' dwt.SW1 <- dwt.vt(data.SW1[[1]],
#'   wf = "d4", J = 7, method = "dwt",
#'   pad = "zero", boundary = "periodic", cov.opt = "auto"
#' )
#'
#' # plot
#' fig1 <- fig.dwt.vt(dwt.SW1)
#' fig1
#'
#' # variance transfrom
#' dwt.HL <- dwt.vt(data.HL[[1]],
#'   wf = "d4", J = 7, method = "dwt",
#'   pad = "zero", boundary = "periodic", cov.opt = "auto"
#' )
#'
#' # plot
#' fig2 <- fig.dwt.vt(dwt.HL)
#' fig2
fig.dwt.vt <- function(dwt.data) {
  x <- dwt.data$x
  dp <- dwt.data$dp
  dp.n <- dwt.data$dp.n
  n <- nrow(dp)
  ndim <- ncol(dp)

  wf <- dwt.data$wavelet
  method <- dwt.data$method
  boundary <- dwt.data$boundary
  pad <- dwt.data$pad
  if (wf != "haar") v <- as.integer(as.numeric(substr(wf, 2, 3)) / 2) else v <- 1
  J <- ceiling(log(n / (2 * v - 1)) / log(2))

  # variance structure of response in spectrum
  mra.x <- waveslim::mra(padding(x, pad), wf = wf, J = J, method = method, boundary = boundary)
  idwt.x <- lapply(mra.x, function(z) z[1:n])

  x.Dj <- c(unlist(lapply(idwt.x, var)) / var(x))
  x.Dj <- data.frame(matrix(rep(x.Dj, ndim), ncol = ndim))
  colnames(x.Dj) <- paste0("X", 1:ndim)
  x.Dj$Group <- 1

  # variance structure of predictors in spectrum
  idwt.dp <- lapply(1:ndim, function(i) waveslim::mra(padding(dp[, i], pad), wf = wf, J = J, method = method, boundary = boundary))
  idwt.dp.n <- lapply(1:ndim, function(i) waveslim::mra(padding(dp.n[, i], pad), wf = wf, J = J, method = method, boundary = boundary))

  dpred.dp <- lapply(seq_len(length(idwt.dp)), function(i) lapply(idwt.dp[[i]], function(z) z[1:n]))
  dpred.dp.n <- lapply(seq_len(length(idwt.dp.n)), function(i) lapply(idwt.dp.n[[i]], function(z) z[1:n]))

  dp.Dj <- data.frame(vapply(seq_len(length(dpred.dp)), function(i) c(unlist(lapply(dpred.dp[[i]], var))) / var(dp[, i]), numeric(J + 1)))
  dp.Dj$Group <- 2
  dp.Dj.n <- data.frame(vapply(seq_len(length(dpred.dp.n)), function(i) c(unlist(lapply(dpred.dp.n[[i]], var))) / var(dp.n[, i]), numeric(J + 1)))
  dp.Dj.n$Group <- 3

  df <- cbind(Level = rep(seq_len(length(idwt.x)), 3), rbind(x.Dj, dp.Dj, dp.Dj.n))
  df.n <- tidyr::gather(df, Predictor, Value, 2:(ndim + 1))
  df.n$Predictor <- factor(df.n$Predictor, levels = paste0("X", 1:ndim))

  # barplot+lineplot
  fig <- ggplot2::ggplot(df.n[df.n$Group != 1, ], aes(factor(Level), Value)) +
    geom_bar(aes(fill = factor(Group)), position = "dodge", stat = "identity") +
    geom_line(data = df.n[df.n$Group == 1, ], aes(x = factor(Level), y = Value), group = 1) +
    facet_grid(Predictor ~ ., scales = "free", space = "free") +
    scale_y_continuous(breaks = seq(0, 1, 0.5), limits = c(0, 1)) +
    scale_fill_manual(values = c("red", "blue")) +
    xlab("Decomposition level") +
    ylab("Variance (percent)") +
    theme_bw() +
    theme(
      text = element_text(size = 8),
      plot.margin = unit(c(1, 1, 1, 1), "cm"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      strip.text.y = element_text(angle = 90),
      legend.position = "none",
      legend.title = element_blank()
    )

  return(fig)
}
