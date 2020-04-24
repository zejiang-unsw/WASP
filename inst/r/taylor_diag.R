#' Taylor diagram
#'
#' @param ref             numeric vector - the reference values.
#' @param model           numeric vector - the predicted model values.
#' @param add             whether to draw the diagram or just add a point.
#' @param col             the color for the points displayed.
#' @param pch             the type of point to display.
#' @param pos.cor         whether to display only positive (TRUE) or all values of correlation (FALSE).
#' @param xlab            plot x axis labels.
#' @param ylab            plot y axis labels.
#' @param main            title for the plot.
#' @param show.gamma      whether to display standard deviation arcs around the reference point (only for pos.cor=TRUE).
#' @param ngamma          the number of gammas to display (default=3).
#' @param gamma.col       color to use for the gamma arcs (only with pos.cor=TRUE).
#' @param sd.arcs         whether to display arcs along the standard deviation axes (see Details).
#' @param cex.max         scale of range of maximum sd.
#' @param sd.method       whether to use the sample or estimated population SD.
#' @param grad.corr.lines the values for the radial lines for correlation values (see Details).
#' @param pcex            character expansion/size for the plotted points.
#' @param cex.axis        character expansion/size for the axis text.
#' @param normalize       whether to normalize the models so that the reference has a standard deviation of 1.Default TRUE.
#' @param mar             margins - only applies to the pos.cor=TRUE plot.
#' @param ...             Additional arguments passed to plot.
#'
#' @details
#' The Taylor diagram is used to display the quality of model predictions against the reference values, typically direct observations.
#' A diagram is built by plotting one model against the reference, then adding alternative model points. If normalize=TRUE when plotting the first model, remember to set it to TRUE when plotting additional models.
#' Two displays are available. One displays the entire range of correlations from -1 to 1. Setting pos.cor to FALSE will produce this display. The -1 to 1 display includes a radial grid for the correlation values. When pos.cor is set to TRUE, only the range from 0 to 1 will be displayed. The gamma lines and the arc at the reference standard deviation are optional in this display.
#' Both the standard deviation arcs and the gamma lines are optional in the pos.cor=TRUE version. Setting sd.arcs or grad.corr.lines to zero or FALSE will cause them not to be displayed. If more than one value is passed for sd.arcs, the function will try to use the values passed, otherwise it will call pretty to calculate the values.
#'
#' @return  The values of par that preceded the function. This allows the user to add points to the diagram, then restore the original values. This is only necessary when using the 0 to 1 correlation range.
#' @export
#'
#' @author Olivier Eterradossi with modifications by Jim Lemon
#' @references Taylor, K.E. (2001) Summarizing multiple aspects of model performance in a single diagram. Journal of Geophysical Research, 106: 7183-7192.
#'
#' @examples
#' # fake some reference data
#' ref<-rnorm(30,sd=2)
#' # add a little noise
#' model = model1<-ref+rnorm(30)/2
#' # add more noise
#' model2<-ref+rnorm(30)
#' par(pty="s")
#' # display the diagram with the better model
#' oldpar<-taylor.diag(ref,model1,cex.max = 1.2)
#' # now add the worse model
#' taylor.diag(ref,model2,add=TRUE,col="blue")
#' # get approximate legend position
#' lpos<-1.5*sd(ref)
#' # add a legend
#' legend(lpos,lpos,legend=c("Better","Worse"),pch=19,col=c("red","blue"))
#' # now restore par values
#' par(oldpar)
#' # show the "all correlation" display
#' taylor.diag(ref,model1,pos.cor=FALSE,cex.max = 1.2)
#' taylor.diag(ref,model2,add=TRUE,col="blue")

taylor.diag <- function (ref, model, add = FALSE, col = "red", pch = 19, pos.cor = TRUE,
          xlab = "Standard deviation", ylab = "", main = "Taylor Diagram",
          show.gamma = TRUE, ngamma = 3, gamma.col = 8, sd.arcs = 0, cex.max=1.5,
          sd.method = "sample", grad.corr.lines = c(0.2,0.4, 0.6, 0.8, 0.9),
          pcex = 1, cex.axis = 1, normalize = TRUE, #should be true otherwise sd of obs won't at the same location
          mar = c(4, 3, 4, 3), ...)
{
  grad.corr.full <- c(0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99, 1)
  R <- cor(ref, model, use = "pairwise")
  if (is.list(ref))
    ref <- unlist(ref)
  if (is.list(model))
    ref <- unlist(model)
  SD <- function(x, subn) { #RMSE
    meanx <- mean(x, na.rm = TRUE)
    devx <- x - meanx
    ssd <- sqrt(sum(devx * devx, na.rm = TRUE)/(length(x[!is.na(x)]) - subn)) #Sn or Sn-1
    return(ssd)
  }
  subn <- sd.method != "sample"
  sd.r <- SD(ref, subn)
  sd.f <- SD(model, subn)
  if (normalize) {
    sd.f <- sd.f/sd.r
    sd.r <- 1
  }
  maxsd <- cex.max * max(sd.f, sd.r)
  oldpar <- par("mar", "xpd", "xaxs", "yaxs")
  if (!add) {
    par(mar = mar)
    if (pos.cor) {
      if (nchar(ylab) == 0)
        ylab = "Standard deviation"
      plot(0, xlim = c(0, maxsd * 1.1), ylim = c(0, maxsd *
                                                   1.1), xaxs = "i", yaxs = "i", axes = FALSE,
           main = main, xlab = "", ylab = ylab, type = "n",
           cex = cex.axis, ...)
      mtext(xlab, side = 1, line = 2.3)
      if (grad.corr.lines[1]) {
        for (gcl in grad.corr.lines) lines(c(0, maxsd *
                                               gcl), c(0, maxsd * sqrt(1 - gcl^2)), lty = 3)
      }
      segments(c(0, 0), c(0, 0), c(0, maxsd), c(maxsd,
                                                0))
      axis.ticks <- pretty(c(0, maxsd))
      axis.ticks <- axis.ticks[axis.ticks <= maxsd]
      axis(1, at = axis.ticks, cex.axis = cex.axis)
      axis(2, at = axis.ticks, cex.axis = cex.axis)
      if (sd.arcs[1]) {
        if (length(sd.arcs) == 1)
          sd.arcs <- axis.ticks
        for (sdarc in sd.arcs) {
          xcurve <- cos(seq(0, pi/2, by = 0.03)) * sdarc
          ycurve <- sin(seq(0, pi/2, by = 0.03)) * sdarc
          lines(xcurve, ycurve, col = "blue", lty = 3)
        }
      }
      if (show.gamma[1]) {
        if (length(show.gamma) > 1)
          gamma <- show.gamma
        else gamma <- pretty(c(0, maxsd), n = ngamma)[-1]
        if (gamma[length(gamma)] > maxsd)
          gamma <- gamma[-length(gamma)]
        labelpos <- seq(45, 70, length.out = length(gamma))
        for (gindex in 1:length(gamma)) {
          xcurve <- cos(seq(0, pi, by = 0.03)) * gamma[gindex] +
            sd.r
          endcurve <- which(xcurve < 0)
          endcurve <- ifelse(length(endcurve), min(endcurve) -
                               1, 105)
          ycurve <- sin(seq(0, pi, by = 0.03)) * gamma[gindex]
          maxcurve <- xcurve * xcurve + ycurve * ycurve
          startcurve <- which(maxcurve > maxsd * maxsd)
          startcurve <- ifelse(length(startcurve), max(startcurve) +
                                 1, 0)
          lines(xcurve[startcurve:endcurve], ycurve[startcurve:endcurve], lty=3, #"dotted"
                col = gamma.col)
          if (xcurve[labelpos[gindex]] > 0)
            boxed.labels(xcurve[labelpos[gindex]], ycurve[labelpos[gindex]],
                         gamma[gindex], border = FALSE)
        }
      }
      xcurve <- cos(seq(0, pi/2, by = 0.01)) * maxsd
      ycurve <- sin(seq(0, pi/2, by = 0.01)) * maxsd
      lines(xcurve, ycurve)
      bigtickangles <- acos(seq(0.1, 0.9, by = 0.1))
      medtickangles <- acos(seq(0.05, 0.95, by = 0.1))
      smltickangles <- acos(seq(0.91, 0.99, by = 0.01))
      segments(cos(bigtickangles) * maxsd, sin(bigtickangles) *
                 maxsd, cos(bigtickangles) * 0.97 * maxsd, sin(bigtickangles) *
                 0.97 * maxsd)
      par(xpd = TRUE)
      points(sd.r, 0, cex = pcex, pch = 22, bg = "darkgreen")
      text(cos(c(bigtickangles, acos(c(0.95, 0.99)))) *
             1.05 * maxsd, sin(c(bigtickangles, acos(c(0.95,
                                                       0.99)))) * 1.05 * maxsd, c(seq(0.1, 0.9, by = 0.1),
                                                                                  0.95, 0.99), cex = cex.axis)
      text(maxsd * 0.8, maxsd * 0.8, "Correlation", srt = 315,
           cex = cex.axis)
      segments(cos(medtickangles) * maxsd, sin(medtickangles) *
                 maxsd, cos(medtickangles) * 0.98 * maxsd, sin(medtickangles) *
                 0.98 * maxsd)
      segments(cos(smltickangles) * maxsd, sin(smltickangles) *
                 maxsd, cos(smltickangles) * 0.99 * maxsd, sin(smltickangles) *
                 0.99 * maxsd)
    }
    else {
      x <- ref
      y <- model
      R <- cor(x, y, use = "pairwise.complete.obs")
      E <- mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)
      xprime <- x - mean(x, na.rm = TRUE)
      yprime <- y - mean(y, na.rm = TRUE)
      sumofsquares <- (xprime - yprime)^2
      Eprime <- sqrt(sum(sumofsquares)/length(complete.cases(x)))
      E2 <- E^2 + Eprime^2
      if (add == FALSE) {
        maxray <- cex.max * max(sd.f, sd.r)
        plot(c(-maxray, maxray), c(0, maxray), type = "n",
             asp = 1, bty = "n", xaxt = "n", yaxt = "n",
             xlim = c(-1.1 * maxray, 1.1 * maxray), xlab = xlab,
             ylab = ylab, main = main, cex = cex.axis)
        discrete <- seq(180, 0, by = -1)
        listepoints <- NULL
        for (i in discrete) {
          listepoints <- cbind(listepoints, maxray *
                                 cos(i * pi/180), maxray * sin(i * pi/180))
        }
        listepoints <- matrix(listepoints, 2, length(listepoints)/2)
        listepoints <- t(listepoints)
        lines(listepoints[, 1], listepoints[, 2])
        lines(c(-maxray, maxray), c(0, 0))
        lines(c(0, 0), c(0, maxray))
        for (i in grad.corr.lines) {
          lines(c(0, maxray * i), c(0, maxray * sqrt(1 -
                                                       i^2)), lty = 3)
          lines(c(0, -maxray * i), c(0, maxray * sqrt(1 -
                                                        i^2)), lty = 3)
        }
        for (i in grad.corr.full) {
          text(1.05 * maxray * i, 1.05 * maxray * sqrt(1 -
                                                         i^2), i, cex = cex.axis, adj = cos(i)/2)
          text(-1.05 * maxray * i, 1.05 * maxray * sqrt(1 -
                                                          i^2), -i, cex = cex.axis, adj = 1 - cos(i)/2)
        }
        seq.sd <- seq.int(0, 2 * maxray, by = (maxray/10))[-1]
        for (i in seq.sd) {
          xcircle <- sd.r + (cos(discrete * pi/180) *
                               i)
          ycircle <- sin(discrete * pi/180) * i
          for (j in 1:length(xcircle)) {
            if ((xcircle[j]^2 + ycircle[j]^2) < (maxray^2)) {
              points(xcircle[j], ycircle[j], col = "darkgreen",
                     pch = ".")
              if (j == 10)
                text(xcircle[j], ycircle[j], signif(i,
                                                    2), cex = cex.axis, col = "darkgreen",
                     srt = 90)
            }
          }
        }
        seq.sd <- seq.int(0, maxray, length.out = 5)
        for (i in seq.sd) {
          xcircle <- cos(discrete * pi/180) * i
          ycircle <- sin(discrete * pi/180) * i
          if (i)
            lines(xcircle, ycircle, lty = 3, col = "blue")
          text(min(xcircle), -0.06 * maxray, signif(i,
                                                    2), cex = cex.axis, col = "blue")
          text(max(xcircle), -0.06 * maxray, signif(i,
                                                    2), cex = cex.axis, col = "blue")
        }
        text(0, -0.14 * maxray, "Standard Deviation",
             cex = cex.axis, col = "blue")
        text(0, -0.22 * maxray, "Centered RMS Difference",
             cex = cex.axis, col = "darkgreen")
        points(sd.r, 0, pch = 22, bg = "darkgreen",
               cex = pcex)
        text(0, 1.2 * maxray, "Correlation Coefficient",
             cex = cex.axis)
      }
      S <- (2 * (1 + R))/(sd.f + (1/sd.f))^2
    }
  }
  points(sd.f * R, sd.f * sin(acos(R)), pch = pch, col = col, cex = pcex)
  invisible(oldpar)
}


