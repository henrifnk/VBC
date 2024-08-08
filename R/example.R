#' \dontrun{
#' climate_sub = lapply(climate, function(x) x[, c("pr", "tas", "dew")])
#' margins_controls = list(xmin = c(0, NaN, NaN), type = c("zi", "c", "c"))
#' mp_vbc = vbc(climate_sub$oc, climate_sub$mc, climate_sub$mp,
#'              margins_controls = margins_controls, family_set = "tll",
#'              trunc_lvl = Inf)
#' class(mp_vbc)
#' summary(attr(mp_vbc, "vine_oc"))
#' summary(attr(mp_vbc, "vine_mp"))
#' calc_mci(climate_sub$mp, mp_vbc, time = climate$mp$time,
#'          margins_controls = margins_controls,
#'          family_set = "tll", trunc_lvl = Inf)
#' rbind(
#'  uncorrected = calc_wasserstein(climate_sub$oc[1:5000, ], climate_sub$mp[1:5000,]),
#'  corrected = calc_wasserstein(climate_sub$oc[1:5000, ], mp_vbc[1:5000, ])
#' )
#' }
