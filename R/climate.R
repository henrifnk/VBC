#' Sample CRCM5 and SDCLIREF data
#'
#' Sample CRCM5-LE member "kba" and SDCLIREF v1 data for the city of Munich
#' (48°8'14.75"N, 11°34'31.76"E).
#'
#' Data are available in 3 hourly resolution from:
#'
#' 01.01.1991 - 31.12.2000 for calibration
#'
#' 01.01.2001 - 31.12-2010 for projection.
#'
#' With 5 climate variables to correct:
#'
#' tas: average surface temperature (deg. C)
#'
#' pr: precipitation (mm day-1)
#'
#' dew: dewpoint temperature (deg. C)
#'
#' rsds: surface downwelling shortwave radiation (W m-2)
#'
#' sfcWind: surface wind speed (m s-1)
#'
#' time: date and time in three hourly resolution
#'
#' @references
#'
#' Leduc, M., Mailhot, A., Frigon, A., Martel, J. L., Ludwig, R.,
#' Brietzke, G. B., ... & Scinocca, J. (2019). The ClimEx project: A 50-member
#' ensemble of climate change projections at 12-km resolution over Europe and
#' northeastern North America with the Canadian Regional Climate Model (CRCM5).
#' Journal of Applied Meteorology and Climatology, 58(4), 663-693.
#'
#' Wood, R. R., Willkofer, F., Schmid, F. J., Trentini, F., Komischke, H., &
#' Ludwig, R. (2017, April). SDCLIREF-A sub-daily reference dataset.
#' In EGU General Assembly Conference Abstracts (p. 15739).
#'
#' @format A list with four elements
#' \describe{
#'   \item{mp}{Simulation data from the model during projection period.}
#'   \item{mc}{Simulation data from the model during calibration period.}
#'   \item{rp}{Measured observations during projection period.}
#'   \item{rc}{Measured observations during calibration period.}
#' }
"climate"
