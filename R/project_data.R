#' @title Delta mapping
#' @description
#' Map the rank delta in model between calibration period `mc` and projection
#' period `mp` on corrected projection period `mph` to correct for climate
#' trends. For climate variables bounded at 0, a multiplicative delta scheme is
#' applied when the projection quantiles are lower than the calibration
#' quantiles.
#'
#' @param mp [double]\cr
#'  uncorrected climate variable in projection period.
#' @param mph [double]\cr
#'  corrected climate variable in projection period.
#' @param mp_kde [kde1d::kde1d]\cr
#'  a kernel density estimation of the climate variable in projection period.
#' @param mc_kde [kde1d::kde1d]\cr
#'  a kernel density estimation of the climate variable in calibration period.
#' @param xmin double(1)\cr
#'  A vector indicating if `xmin` is a ratio variable type or any other type (NA).
#' @return A climate variable that is corrected by the climate trend in the
#' model between correction and projection period.
map_delta <- function(mp, mph, mp_kde, mc_kde, xmin = NA) {
  mpu <- pkde1d(mp, mp_kde)
  mc_p <- qkde1d(mpu, mc_kde)
  if(is.na(xmin)) {
    delta <- mp - mc_p
    return(mph + delta)
  } else {
    delta_rat <- mp / mc_p
    delta_rat[is.na(delta_rat)] <- 1
    mph_c <- mph * delta_rat
    delta_add <- mp - mc_p
    mph_c[delta_rat > 1] <- mph[delta_rat > 1] + delta_add[delta_rat > 1]
    return(mph_c)
  }
}
