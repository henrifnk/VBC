#' \dontrun{
#' # Start with a simple example:
#' # There is one timeline that should be corrected and projection and 
#' # calibaration period are both in 2010.
#' set.seed(1234L)
#' library(VBC)
#' library(data.table)
#' library(patchwork)
#' data("climate")
#' margins_controls <- list(xmin = c(NaN, 0, NaN, 0, 0),
#'                          type = c("c", "zi", "c", "zi", "c"))
#'                   
#' climate_2010 = lapply(climate, function(data) data[year(time) == 2010, ])
#' 
#' mp_vbc = vbc(climate_2010$mp[, -"time"], climate_2010$mp[, -"time"],
#'              climate_2010$rp[, -"time"], margins_controls = margins_controls,
#'              family_set = "tll", trunc_lvl = Inf)                 
#' 
#' plot_tails(climate_2010$mp, "pr", scale_d = 0.1, mult = 4, xmin = 0) 
#' plot_tails(climate_2010$rp, "pr", scale_d = 1, mult = 4, xmin = 0)
#' plot_tails(round(mp_vbc, 3), "pr", scale_d = 1, mult = 3, xmin = 0)
#' calc_wasserstein(climate_2010$mp[, "pr"], climate_2010$rp[, "pr"])
#' calc_wasserstein(climate_2010$rp[, "pr"], mp_vbc[, "pr"])
#' 
#' calc_wasserstein(climate_2010$rp[, -"time"], climate_2010$mp[, -"time"])
#' calc_wasserstein(climate_2010$rp[, -"time"], mp_vbc)
#' 
#' # An example with temporal subsetting using the method of fragments
#'  
#' temp_subs <- list(DJF = list(hours = seq(0, 21, 3), month = c(12, 1, 2)),
#'                   MAM = list(hours = seq(0, 21, 3), month = 3:5),
#'                   JJA = list(hours = seq(0, 21, 3), month = 6:8),
#'                   SON = list(hours = seq(0, 21, 3), month = 9:11))  
#'                                  
#' mp_vbc = vbc_tsub(climate$mp, climate$mc, climate$rc, t_subs = temp_subs,
#'                   margins_controls = margins_controls, family_set = "tll",
#'                   trunc_lvl = Inf)
#' class(mp_vbc)
#' attr(mp_vbc, "mvd")
#' 
#' measure = lapply(temp_subs, function(sub) {
#'    idx = which(hour(mp_vbc$time) %in% sub$hour &
#'                   month(mp_vbc$time) %in% sub$month)
#'    idxrp = which(hour(climate$rp$time) %in% sub$hour &
#'                      month(climate$rp$time) %in% sub$month)
#'    mci = calc_mci(climate$mp[idx, -"time"], mp_vbc[idx, -"time"],
#'                     time = climate$mp$time[idx])
#'    mci = attr(mci, "global_mci")
#'    
#'    uncorrected = calc_wasserstein(climate$rp[idxrp, -"time"],
#'                                     climate$mp[idx, -"time"])
#'    corrected = calc_wasserstein(climate$rp[idxrp, -"time"],
#'                                 mp_vbc[idx, -"time"])
#'    data.table(
#'        ModelCorrectionInconsistancy = mci,
#'        Wasserstein2_uncorrected = uncorrected[2],
#'        Wasserstein2_corrected = corrected[2],
#'        Improvement_Wasserstein2 = uncorrected[2] - corrected[2], 
#'        month = paste(sub$month, collapse = "-")
#'        )
#' })
#' rbindlist(measure)
#' }
