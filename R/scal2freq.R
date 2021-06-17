#' Scale to frequency by Matlab
#'
#' @param wf wavelet name
#' @param scale a scale
#' @param delta the sampling period.
#'
#' @return
#' @export
#'
#' @examples
scal2freqM <- function(wf, scale, delta){

  if(TRUE){
    # central frequency for other wavelet filters can be obtained by centfrq(wf,8,'plot') in matlab
    centfrq <- switch(wf,
           haar= 0.9961,
           d4=0.6667,
           d6=0.8000,
           d8=0.7143,
           d10=0.6667,
           d12=0.7273,
           d14=0.6923,
           d16=0.6667,
           d18=0.7059,
           d20=0.6842,
           warning("Only daubechies wavelet family is available for now!"))

  }

  freq <- centfrq/(scale*delta)
  period <- 1/freq

  return(list(frequency=freq,
              period=period))
}


#' Scale to frequency by R
#'
#' @param wf wavelet name
#' @param scale a scale
#' @param delta the sampling period.
#'
#' @return
#' @export
#'
#' @examples
scal2freqR <- function(wf, scale, delta){

  if(wf %in% c("haar", "d4", "d6", "d8", "d16")){

    wavelet <- wavelet.filter(wf, "HLLLLL")
    n <- length(wavelet)
    xmax <- ifelse(wf!='haar', as.integer(readr::parse_number(wf))-1, 1)
    xmin <- 0
    xval <- seq(xmin, xmax, length.out=n)

    wavelet.fft <- fft(scale(wavelet, scale=F))
    mod <- Mod(wavelet.fft)
    sp <- abs(mod)
    indmax <- which(mod==max(sp))[1]
    if(indmax>n/2) indmax = n-indmax+2

    per = (xmax-xmin)/(indmax-1)
    centfrq = 1/per
  } else {
    # only daubechies wavelet of haar, d4, d6, d8, d16 available in waveslim
    stop("Invalid selection for wave.filter in waveslim package")
  }

  freq <- centfrq/(scale*delta)
  period <- 1/freq

  return(list(frequency=freq,
              period=period))
}
