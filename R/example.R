#' \dontrun{
#' data("climate")
#' climate_sub = lapply(climate, function(x) x[, c("pr", "tas", "dew")])
#' margins_controls = list(xmin = c(0, NaN, NaN), type = c("zi", "c", "c"))
#' mp_vbc = vine_correct(climate_sub$oc, climate_sub$mc, climate_sub$mp,
#'                       margins_controls = margins_controls,
#'                       family_set = "tll", trunc_lvl = Inf)
#' summary(attr(mp_vbc, "vine_oc"))
#' summary(attr(mp_vbc, "vine_mp"))
#' calc_mci(climate_sub$mp, mp_vbc, margins_controls = margins_controls,
#'          family_set = "tll", trunc_lvl = Inf)
#' calc_wasserstein(climate_sub$oc, mp_vbc)
#' }
