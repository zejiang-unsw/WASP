#' Generate predictor and response data: Rossler system
#' @description
#' Generates a 3-dimensional time series using the Rossler equations.
#' @details
#' The Rossler system is a system of ordinary differential equations defined as:
#' \deqn{\dot{x} = -(y + z)}{dx/dt = -(y + z)}
#' \deqn{\dot{y} = x+a \cdot y}{dy/dt = x + a*y}
#' \deqn{\dot{z} = b + z*(x-w)}{dz/dt = b + z*(x-w)}
#' The default selection for the system parameters (\emph{a} = 0.2, \emph{b} = 0.2, \emph{w} = 5.7) is known to
#' produce a deterministic chaotic time series.
#'
#' @param start A 3-dimensional numeric vector indicating the starting point for the time series. Default: c(-2, -10, 0.2).
#' @param a The \emph{a} parameter. Default:0.2.
#' @param b The \emph{b} parameter. Default: 0.2.
#' @param w The \emph{w} parameter. Default: 5.7.
#' @param time The temporal interval at which the system will be generated.
#' Default: time=seq(0,50,length.out = 5000).
#'
#' @return A list with four vectors named \emph{time}, \emph{x}, \emph{y}
#' and \emph{z} containing the time, the x-components, the
#' y-components and the z-components of the Rossler system, respectively.
#' @export
#'
#' @note Some initial values may lead to an unstable system that will tend to infinity.
#' @references RÖSSLER, O. E. 1976. An equation for continuous chaos. Physics Letters A, 57, 397-398.
#' @examples
#' ###synthetic example - Rossler
#' ts.r <- data.gen.Rossler(a = 0.2, b = 0.2, w = 5.7, start = c(-2, -10, 0.2),
#'                 time = seq(0, 50, length.out = 1000))
#'
#' #add noise
#' ts.r$x <- ts(ts.r$x + rnorm(length(ts.r$time),mean=0, sd=1))
#' ts.r$y <- ts(ts.r$y + rnorm(length(ts.r$time),mean=0, sd=1))
#' ts.r$z <- ts(ts.r$z + rnorm(length(ts.r$time),mean=0, sd=1))
#'
#' ts.plot(ts.r$x,ts.r$y,ts.r$z, col=c("black","red","blue"))

data.gen.Rossler <- function(a = 0.2, b = 0.2, w = 5.7, start=c(-2, -10, 0.2),
                             time = seq(0, 50, length.out = 5000)) {
  params = c(a, b, w)
  rosslerEquations = function(coord, t, params) {
    x = coord[[1]]
    y = coord[[2]]
    z = coord[[3]]
    a = params[[1]]
    b = params[[2]]
    w = params[[3]]
    c(-y - z, x + a * y, b + z * (x - w))
  }
  r = rungeKutta(rosslerEquations, start, time, params)

  list(time = time, x = r[, 1], y = r[, 2], z = r[, 3])
}

# Runge-Kutta method for solving differential equations. It is used to generate
# both Lorenz and  Rossler systems.
rungeKutta = function(func, initial.condition, time, params) {
  n.samples = length(time)
  h = time[[2]] - time[[1]]
  y = matrix(ncol = length(initial.condition), nrow = n.samples)
  y[1,] = initial.condition
  for (i in 2:n.samples) {
    k1 = h * func(y[i - 1, ], time[[i - 1]], params)
    k2 = h * func(y[i - 1, ] + k1 / 2 , time[[i - 1]] + h / 2, params)
    k3 = h * func(y[i - 1, ] + k2 / 2 , time[[i - 1]] + h / 2, params)
    k4 = h * func(y[i - 1, ] + k3 , time[[i - 1]] + h, params)

    y[i, ] = y[i - 1, ] + (k1 + 2 * k2 + 2 * k3 + k4) / 6
  }
  y
}
