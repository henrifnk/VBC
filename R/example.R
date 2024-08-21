#' \dontrun{
#' data("climate")
#' climate_sub = lapply(climate, function(x) x[, c("pr", "tas", "dew", "time")])
#' margins_controls = list(xmin = c(0, NaN, NaN), type = c("zi", "c", "c"))
#' t_subs = list(day = list(hours = c(6, 9, 12, 15), month = 1:12),
#'               night = list(hours = c(18, 21, 0, 3), month = 1:12))
#' mp_vbc = vbc_tsub(climate_sub$rc, climate_sub$mc, climate_sub$mp,
#'                   t_subs = t_subs, margins_controls = margins_controls,
#'                   family_set = "tll", trunc_lvl = Inf)
#' class(mp_vbc)
#' attr(mp_vbc, "mvd")
#' 
#' calc_mci(climate_sub$mp[, -"time"], mp_vbc[, -"time"],
#'          time = climate$mp$time, margins_controls = margins_controls,
#'          family_set = "tll", trunc_lvl = Inf)
#' rbind(
#'  uncorrected = calc_wasserstein(climate_sub$rp[1:5000, -"time"],
#'                                 climate_sub$mp[1:5000, -"time"]),
#'  corrected = calc_wasserstein(climate_sub$rp[1:5000, -"time"],
#'                               mp_vbc[1:5000, -"time"])
#' )
#' }
