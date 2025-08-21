#' Combine VBC and the Method of fragments
#' 
#' @param rc [data.table::data.table]\cr
#'  Measured (and interpolated) observations during calibration period. 
#'  Expects one column time.
#' 
#' @param mc [data.table::data.table]\cr
#'  Simulation data from a climate model during calibration period.
#'  Expects one column time.
#'
#' @param mp [data.table::data.table]\cr
#'  Simulation data from a climate model during projection period.
#'  Expects one column time.
#'  
#' @inheritParams vbc
#' 
#' @param t_subs [list]\cr
#' A list of two lists: hours and month. Each list contains the temporal
#' indicators to subset the data. The first list contains the hours of the day
#' the length of list marks the number of fragments. The inner list contains the
#' hour and monthly information to subset the data.
#' 
#' @param overlap [integer]\cr
#' The number of times to overlap the temporal indicators.
#' 
#' @param cores_t [integer]\cr
#' The number of cores to use for parallel processing of the temporal 
#' subsetting. Default is `NA` which means no parallel processing.
#' 
#' @param verbose [logical]\cr
#' Print messages during the temporal subsetting.
#' 
#' @importFrom parallel mclapply
#' @importFrom lubridate as.duration
#' 
#' @references 
#' Sharma, A.; Srikanthan, S. (2006): Continuous Rainfall Simulation: 
#' A Nonparametric Alternative. In: 30th Hydrology and Water Resources
#' Symposium 4-7 December 2006, Launceston, Tasmania.
#' 
#' Westra, S.; Mehrotra, R.; Sharma, A.; Srikanthan, R. (2012): Continuous
#' rainfall simulation. 1. A regionalized subdaily disaggregation approach. 
#' In: Water Resour. Res. 48 (1). DOI: 10.1029/2011WR010489.
#' 
#' @return [data.table::data.table]\cr
#'  The corrected projection period data in `mp`. Additionally the data frame
#'  contains the attribute `mvd` with one list per temporal subset. In each
#'  subset `vine_rc`, `kde_rc`, `vine_mp`, and `kde_mp` which store the 
#'  vine copula and kernel density estimation objects of the observed and model
#'  data. The time column is attached to the data. 
#' 
#' @example R/example.R
#' 
#' @importFrom stats setNames
#' 
#' @export
vbc_tsub = function(mp, mc, rc, var_names = colnames(rc),
                    margins_controls = list(
                      mult = NULL, xmin = NaN, xmax = NaN, bw = NA, deg = 2,
                      type = "c"),
                    t_subs = list(
                      list(hours = 0:23, month = 1:12)
                    ), overlap = 1, cores_t = NA, verbose = TRUE, ...) {
  if(is.na(cores_t)) {
    tmp_mph <- lapply(t_subs, function(t_sub) {
      rc_sub <- subset_time(rc, t_sub$hours, t_sub$month, overlap)
      mc_sub <- subset_time(mc, t_sub$hours, t_sub$month, overlap)
      if(is.data.frame(mp)) {
        mp_sub <- subset_time(mp, t_sub$hours, t_sub$month, overlap)
      } else {
        mp_sub <- lapply(mp, subset_time,  hrs = t_sub$hours, mnt = t_sub$month,
                         overlap = overlap)
      }
      mph <- vbc(mp_sub, mc_sub, rc_sub, margins_controls =  margins_controls,
                 ...)
      if(verbose) {
        message("subset for ", t_sub$hours, " hours and ", t_sub$month,
                " month done.")
      }
      if(is.data.frame(mp)) {
        final_idx <- attr(mp_sub, "final_idx")
        idx <- attr(mp_sub, "idx")
        mph[, "idx" := idx]
        mph[idx %in% final_idx, ]
      } else {
        mph <- mapply(function(mem_mph, mem_mp) {
          final_idx <- attr(mem_mp, "final_idx")
          idx <- attr(mem_mp, "idx")
          mem_mph[, "idx" := idx]
          mem_mph[idx %in% final_idx, ]
        }, mem_mph = mph, mem_mp = mp_sub, SIMPLIFY = FALSE)
        mph
      }
    })
  } else {
    tmp_mph <- mclapply(t_subs, function(t_sub) {
      rc_sub <- subset_time(rc, t_sub$hours, t_sub$month, overlap)
      mc_sub <- subset_time(mc, t_sub$hours, t_sub$month, overlap)
      if(is.list(mp)) {
        mp_sub <- lapply(mp, subset_time,  hrs = t_sub$hours, mnt = t_sub$month,
                         overlap = overlap)
      } else {
        mp_sub <- subset_time(mp, t_sub$hours, t_sub$month, overlap)
      }
      final_idx <- attr(mp_sub, "final_idx")
      idx <- attr(mp_sub, "idx")
      mph <- vbc(rc_sub, mc_sub, mp_sub, margins_controls =  margins_controls,
                 ...)
      if(verbose) {
        message("subset for ", t_sub$hours, " hours and ", t_sub$month,
                " month done.")
      }
      if(is.data.frame(mp)) {
        mph[, "idx" := idx]
        mph[idx %in% final_idx, ]
      } else {
        mph <- lapply(mph, function(x) {
          x[, "idx" := idx]
          x[idx %in% final_idx, ]
        })
        mph
      }
    }, mc.cores = cores_t)
  }
  if(is.data.frame(mp)) {
    mph <- rbindlist(tmp_mph)[order(rank(get("idx")))]
    attribs <- lapply(tmp_mph, function(x) {
      names <- c("vine_rc", "kde_rc", "vine_mp", "kde_mp")
      names_list <- setNames(names, names)
      lapply(names_list, function(y) attr(x, y))
    })
    time_mem = mp[["time"]][mph[["idx"]]]
    mph[, "idx" := NULL][, "time" := time_mem]
    class(mph) <- c("vbc_tsub", class(mph))
    attr(mph, "mvd") <- attribs
  } else {
    mem_names <- names(tmp_mph[[1]])
    tmp_mph <- invert_list(tmp_mph)
    names(tmp_mph) <- mem_names
    mph <- mapply(function(tmp_mem, mp_mem) {
      mph <- rbindlist(tmp_mem)[order(rank(get("idx")))]
      attribs <- lapply(tmp_mem, function(x) {
        names <- c("vine_rc", "kde_rc", "vine_mp", "kde_mp")
        names_list <- setNames(names, names)
        lapply(names_list, function(y) attr(x, y))
      })
      time_mem = mp_mem[["time"]][mph[["idx"]]]
      mph[, "idx" := NULL][, "time" := time_mem]
      class(mph) <- c("vbc_tsub", class(mph))
      attr(mph, "mvd") <- attribs
      mph
    }, tmp_mem = tmp_mph, mp_mem = mp, SIMPLIFY = FALSE)
  }
  mph
}

# utilities --------------------------------------------------------------------
#' @title subset time and final reading index
#' 
#' @inheritParams model_vine
#' @inheritParams vbc_tsub
#' 
#' @param hrs [integer]\cr
#' The hours of the day to be used for subsetting.
#' 
#' @param mnt [integer]\cr
#' The months of the year to be used for subsetting.
subset_time <- function(data, hrs, mnt, overlap = 1) {
  time <- data$time
  resolution_t <- as.duration(difftime(time[2], time[1]))
  resolution_hrs <- as.numeric(resolution_t, "hours")
  hrs_over <- calc_range(hrs, min = min(hour(time)), max = max(hour(time)),
                         res = resolution_hrs, overlap)
  mnt_over <- calc_range(mnt, min = 1, max = 12, res = 1, overlap)
  idx <- which(hour(time) %in% hrs_over & month(time) %in% mnt_over)
  final_idx <- which(hour(time) %in% hrs & month(time) %in% mnt)
  data_sub <- data[idx, ]
  #browser()
  # stop if no data is selected
  assert_that(nrow(data_sub) > 0,
              msg = paste("no data selected for hours ", hrs,
                          " and months ", mnt))
  data_sub$time <- NULL
  attr(data_sub, "final_idx") <- final_idx
  attr(data_sub, "idx") <- idx
  data_sub
}

#' @title calculate temporal range recursively by overlap.
#' @param v sorted vector of temporal indicators
#' @param min cyclic minimum
#' @param max cyclic maximum
#' @param res resolution of time
#' @inheritParams vbc_tsub
#' @return a vector of temporal indicators
calc_range = function(v, min, max, res, overlap = 1) {
  if(overlap == 0) return(unique(v))
  range <- c(v[1] - res, v, v[length(v)] + res)
  range[1] <- ifelse(range[1] < min, max, range[1])
  range[length(range)] <- ifelse(range[length(range)] > max,
                                 min, range[length(range)])
  calc_range(range, min, max, res, overlap = overlap - 1L)
}

#' @title Invert a list
#' 
#' @param list [list]\cr
#' A list of depth 2. The first depth indicates the temporal subset and the
#' second depth the members of the model ensemble.
#' 
#' @return [list]\cr
#' A list of depth 2. The first depth indicates the members of the model
#' ensemble and the second depth the temporal subset.
#'
invert_list = function(list) {
  lapply(seq_along(list[[1]]), function(i) {
    lapply(list, function(x) x[[i]])
  })
}
