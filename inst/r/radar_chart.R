#' Radar Chart
#'
#' @param df The data frame to be used to draw radarchart. If maxmin is TRUE, this must include maximum values as row 1 and minimum values as row 2 for each variables, and actual data should be given as row 3 and lower rows. The number of columns (variables) must be more than 2.
#' @param axistype The type of axes, specified by any of 0:5. 0 means no axis label. 1 means center axis label only. 2 means around-the-chart label only. 3 means both center and around-the-chart (peripheral) labels. 4 is *.** format of 1, 5 is *.** format of 3. Default is 0.
#' @param seg The number of segments for each axis (default 4).
#' @param pty A vector to specify point symbol: Default 16 (closed circle), if you don't plot data points, it should be 32. This is repeatedly used for data series.
#' @param pcol A vector of color codes for plot data: Default 1:8, which are repeatedly used.
#' @param plty A vector of line types for plot data: Default 1:6, which are repeatedly used.
#' @param plwd A vector of line widths for plot data: Default 1, which is repeatedly used.
#' @param pdensity A vector of filling density of polygons: Default NULL, which is repeatedly used.
#' @param pangle A vector of the angles of lines used as filling polygons: Default 45, which is repeatedly used.
#' @param pfcol A vector of color codes for filling polygons: Default NA, which is repeatedly usd.
#' @param cglty Line type for radar grids: Default 3, which means dotted line.
#' @param cglwd Line width for radar grids: Default 1, which means thinnest line.
#' @param cglcol Line color for radar grids: Default "navy".
#' @param axislabcol Color of axis label and numbers: Default "blue"
#' @param title if any, title should be typed.
#' @param maxmin Logical. If true, data frame includes possible maximum values as row 1 and possible minimum values as row 2. If false, the maximum and minimum values for each axis will be calculated as actual maximum and minimum of the data. Default TRUE.
#' @param na.itp Logical. If true, items with NA values are interpolated from nearest neighbor items and connect them. If false, items with NA are treated as the origin (but not pointed, only connected with lines). Default FALSE.
#' @param centerzero Logical. If true, this function draws charts with scaling originated from (0,0). If false, charts originated from (1/segments). Default FALSE.
#' @param vlabels Character vector for the names for variables. If NULL, the names of the variables as colnames(df) are used. Default NULL.
#' @param vlcex Font size magnification for vlabels. If NULL, the font size is fixed at text()'s default. Default NULL.
#' @param caxislabels Character vector for center axis labels, overwriting values specified in axistype option. If NULL, the values specified by axistype option are used. Default is NULL.
#' @param calcex Font size magnification for caxislabels. If NULL, the font size is fixed at text()'s default. Default NULL.
#' @param paxislabels Character vector for around-the-chart (peripheral) labels, overwriting values specified in axistype option. If NULL, the values specified by axistype option are used. Default is NULL.
#' @param palcex Font size magnification for paxislabels. If NULL, the font size is fixed at text()'s default. Default NULL.
#' @param ... Miscellaneous arguments to be given for plot.default().
#'
#' @return No
#' @export
#'
#' @examples
#'
radar.chart <- function (df, axistype = 0, seg = 4, pty = 16, pcol = 1:8, plty = 1:6,
          plwd = 1, pdensity = NULL, pangle = 45, pfcol = NA, cglty = 3,
          cglwd = 1, cglcol = "navy", axislabcol = "blue", title = "",
          maxmin = TRUE, na.itp = TRUE, centerzero = FALSE, vlabels = NULL,
          vlcex = NULL, caxislabels = NULL, calcex = NULL, paxislabels = NULL,
          palcex = NULL, ...)
{
  if (!is.data.frame(df)) {
    cat("The data must be given as dataframe.\n")
    return()
  }
  if ((n <- length(df)) < 3) {
    cat("The number of variables must be 3 or more.\n")
    return()
  }
  if (maxmin == FALSE) {
    dfmax <- apply(df, 2, max)
    dfmin <- apply(df, 2, min)
    df <- rbind(dfmax, dfmin, df)
  }
  plot(c(-1.2, 1.2), c(-1.2, 1.2), type = "n", frame.plot = FALSE,
       axes = FALSE, xlab = "", ylab = "", main = title, asp = 1,
       ...)
  theta <- seq(90, 450, length = n + 1) * pi/180
  theta <- theta[1:n]
  xx <- cos(theta)
  yy <- sin(theta)
  CGap <- ifelse(centerzero, 0, 1)
  for (i in 0:seg) {
    polygon(xx * (i + CGap)/(seg + CGap), yy * (i + CGap)/(seg +
                                                             CGap), lty = cglty, lwd = cglwd, border = cglcol)
    if (axistype == 1 | axistype == 3)
      CAXISLABELS <- paste(i/seg * 100, "(%)")
    if (axistype == 4 | axistype == 5)
      CAXISLABELS <- sprintf("%3.2f", i/seg)
    if (!is.null(caxislabels) & (i < length(caxislabels)))
      CAXISLABELS <- caxislabels[i + 1]
    if (axistype == 1 | axistype == 3 | axistype == 4 |
        axistype == 5) {
      if (is.null(calcex))
        text(-0.05, (i + CGap)/(seg + CGap), CAXISLABELS,
             col = axislabcol)
      else text(-0.05, (i + CGap)/(seg + CGap), CAXISLABELS,
                col = axislabcol, cex = calcex)
    }
  }
  if (centerzero) {
    arrows(0, 0, xx * 1, yy * 1, lwd = cglwd, lty = cglty,
           length = 0, col = cglcol)
  }
  else {
    arrows(xx/(seg + CGap), yy/(seg + CGap), xx * 1, yy *
             1, lwd = cglwd, lty = cglty, length = 0, col = cglcol)
  }
  PAXISLABELS <- df[1, 1:n]
  if (!is.null(paxislabels))
    PAXISLABELS <- paxislabels
  if (axistype == 2 | axistype == 3 | axistype == 5) {
    if (is.null(palcex))
      text(xx[1:n]*1.05, yy[1:n]*1.05, PAXISLABELS, col = axislabcol)
    else text(xx[1:n], yy[1:n], PAXISLABELS, col = axislabcol,
              cex = palcex)
  }
  VLABELS <- colnames(df)
  if (!is.null(vlabels))
    VLABELS <- vlabels
  if (is.null(vlcex))
    text(xx * 1.2, yy * 1.2, VLABELS)
  else text(xx * 1.2, yy * 1.2, VLABELS, cex = vlcex)
  series <- length(df[[1]])
  SX <- series - 2
  if (length(pty) < SX) {
    ptys <- rep(pty, SX)
  }
  else {
    ptys <- pty
  }
  if (length(pcol) < SX) {
    pcols <- rep(pcol, SX)
  }
  else {
    pcols <- pcol
  }
  if (length(plty) < SX) {
    pltys <- rep(plty, SX)
  }
  else {
    pltys <- plty
  }
  if (length(plwd) < SX) {
    plwds <- rep(plwd, SX)
  }
  else {
    plwds <- plwd
  }
  if (length(pdensity) < SX) {
    pdensities <- rep(pdensity, SX)
  }
  else {
    pdensities <- pdensity
  }
  if (length(pangle) < SX) {
    pangles <- rep(pangle, SX)
  }
  else {
    pangles <- pangle
  }
  if (length(pfcol) < SX) {
    pfcols <- rep(pfcol, SX)
  }
  else {
    pfcols <- pfcol
  }
  for (i in 3:series) {
    xxs <- xx
    yys <- yy
    scale <- CGap/(seg + CGap) + (df[i, ] - df[2, ])/(df[1,
                                                         ] - df[2, ]) * seg/(seg + CGap)
    if (sum(!is.na(df[i, ])) < 3) {
      cat(sprintf("[DATA NOT ENOUGH] at %d\n%g\n", i,
                  df[i, ]))
    }
    else {
      for (j in 1:n) {
        if (is.na(df[i, j])) {
          if (na.itp) {
            left <- ifelse(j > 1, j - 1, n)
            while (is.na(df[i, left])) {
              left <- ifelse(left > 1, left - 1, n)
            }
            right <- ifelse(j < n, j + 1, 1)
            while (is.na(df[i, right])) {
              right <- ifelse(right < n, right + 1,
                              1)
            }
            xxleft <- xx[left] * CGap/(seg + CGap) +
              xx[left] * (df[i, left] - df[2, left])/(df[1,
                                                         left] - df[2, left]) * seg/(seg + CGap)
            yyleft <- yy[left] * CGap/(seg + CGap) +
              yy[left] * (df[i, left] - df[2, left])/(df[1,
                                                         left] - df[2, left]) * seg/(seg + CGap)
            xxright <- xx[right] * CGap/(seg + CGap) +
              xx[right] * (df[i, right] - df[2, right])/(df[1,
                                                            right] - df[2, right]) * seg/(seg +
                                                                                            CGap)
            yyright <- yy[right] * CGap/(seg + CGap) +
              yy[right] * (df[i, right] - df[2, right])/(df[1,
                                                            right] - df[2, right]) * seg/(seg +
                                                                                            CGap)
            if (xxleft > xxright) {
              xxtmp <- xxleft
              yytmp <- yyleft
              xxleft <- xxright
              yyleft <- yyright
              xxright <- xxtmp
              yyright <- yytmp
            }
            xxs[j] <- xx[j] * (yyleft * xxright - yyright *
                                 xxleft)/(yy[j] * (xxright - xxleft) -
                                            xx[j] * (yyright - yyleft))
            yys[j] <- (yy[j]/xx[j]) * xxs[j]
          }
          else {
            xxs[j] <- 0
            yys[j] <- 0
          }
        }
        else {
          xxs[j] <- xx[j] * CGap/(seg + CGap) + xx[j] *
            (df[i, j] - df[2, j])/(df[1, j] - df[2,
                                                 j]) * seg/(seg + CGap)
          yys[j] <- yy[j] * CGap/(seg + CGap) + yy[j] *
            (df[i, j] - df[2, j])/(df[1, j] - df[2,
                                                 j]) * seg/(seg + CGap)
        }
      }
      if (is.null(pdensities)) {
        polygon(xxs, yys, lty = pltys[i - 2], lwd = plwds[i -
                                                            2], border = pcols[i - 2], col = pfcols[i -
                                                                                                      2])
      }
      else {
        polygon(xxs, yys, lty = pltys[i - 2], lwd = plwds[i -
                                                            2], border = pcols[i - 2], density = pdensities[i -
                                                                                                              2], angle = pangles[i - 2], col = pfcols[i -
                                                                                                                                                         2])
      }
      points(xx * scale, yy * scale, pch = ptys[i - 2],
             col = pcols[i - 2])
    }
  }
}
