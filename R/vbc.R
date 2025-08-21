# VBC --------------------------------------------------------------------------
#' @title Correct zero inflated climate model simulations via vine copulas
#'
#' @description
#' The multivariate distribution functions of model and observed are estimated
#' via vine copula estimation (see [rvinecopulib::vine]).
#' This method provides the tool to model the multivariate distribution function
#' of the climate data by Vine Copulas. The climate variables may be zero
#' inflated.
#' The quantiles are then mapped via (inverse) Rosenblatt transformation from
#' the model to the observed distribution. The corrected data are mapped to the
#' projection domain via delta mapping. The steps are equivalent to those in
#' univariate bias correction via quantile delta mapping.
#'
#' @param mp [data.table::data.table] OR [list]\cr
#'  Simulation data from a climate model during projection period. If this is a
#'  list, the algorithm expects each element to be a member of a model ensemble.
#'  Each list element is then a data table.
#'
#' @param mc [data.table::data.table]\cr
#'  Simulation data from a climate model during calibration period.
#'
#' @param rc [data.table::data.table]\cr
#'  Historical reference in calibration period.
#'
#' @param var_names [character]\cr
#'  Names of corrected climate data. Defaults to the column names of the
#'  observed data `rc`.
#'
#' @param margins_controls [list]\cr
#' A list with arguments to be passed to [kde1d::kde1d()]. Currently, there can be
#'   * `mult` numeric vector of length one or d; all bandwidths for marginal
#'   kernel density estimation are multiplied with `mult`. Defaults to
#'   `log(1 + d)` where `d` is the number of climate variables.
#'   * `xmin` numeric vector of length d; see [kde1d::kde1d()].
#'   * `xmax` numeric vector of length d; see [kde1d::kde1d()].
#'   * `bw` numeric vector of length d; see [kde1d::kde1d()].
#'   * `deg` numeric vector of length one or d; [kde1d::kde1d()].
#'   * `type` character vector of length one or d; must be one of
#'   c, cont, continuous for continuous variables, one of d, disc, discrete for
#'   discrete integer variables, or one of zi, zinfl, zero-inflated for
#'   zero-inflated variables.
#'
#' @param time_mp [numeric]\cr
#' Time vector of the projection period attached to the returned object as a 
#' column. Defaults to `NA`, which means no time vector is attached.
#'
#' @param ... \cr
#'  Arguments are passed to [rvinecopulib::vinecop] to specify the structure of
#'  vines and margins. Note that the ellipsis of observed and model data are
#'  specified with the same arguments.
#'
#' @return [data.table::data.table] OR [list]\cr
#'  The corrected projection period data in `mp`. Additionally the data frame
#'  contains the attributes `vine_rc`, `kde_rc`, `vine_mp`, and `kde_mp` which
#'  store the vine copula and kernel density estimation objects of the observed
#'  and model data.
#'  If `mp` is a list, the function returns a list of corrected data frames for
#'  each member of the model ensemble.
#'
#' @example R/example.R
#'
#' @references
#' Czado, C. (2019). Analyzing dependent data with vine copulas. Lecture Notes
#' in Statistics, Springer, 222.
#'
#' Rosenblatt, M. (1952). Remarks on a multivariate transformation. The annals
#' of mathematical statistics, 23(3), 470-472.
#'
#' Cannon, A.J., S.R. Sobie, and T.Q. Murdock, (2015). Bias correction of
#' simulated precipitation by quantile mapping: How well do methods preserve
#' relative changes in quantiles and extremes? Journal of Climate, 28:6938-6959.
#'
#' @import checkmate
#' @import rvinecopulib
#' @import kde1d
#' @import data.table
#' @export
vbc <- function(mp, mc, rc, var_names = colnames(rc), margins_controls = list(
  mult = NULL, xmin = NaN, xmax = NaN, bw = NA, deg = 2, type = "c"
), time_mp = NA, ...) {
  check_vbc_args(mp, mc, rc, var_names)
  if(ncol(mc != length(margins_controls[["type"]]))) {
    margins_controls <- expand_margin_controls(margins_controls, mc)
  }
  
  if(!is.data.frame(mp) & is.list(mp)) {
    return(vbc_ensemble(mp, mc, rc, var_names, margins_controls, time_mp, ...))
  }
  mc_kde <- attr(estimate_margins(mc, margins_controls), "kde")
  mpu <- model_vine(mp, margins_controls, ...)
  rcu <- model_vine(rc, margins_controls, ...)
  attr(rcu, "vine")$var_types = rep("c", times = ncol(mp))
  x_mph <- correct_rosenblatt(mpu, rcu)
  xmin = if(length(margins_controls$xmin) != ncol(rc)) {
    rep(NA, times = ncol(rc))
  } else {
    margins_controls$xmin
  }
  xproj <- mapply(map_delta, mp = mp, mph = data.frame(x_mph),
                  mp_kde = attr(mpu, "kde"), mc_kde = mc_kde, xmin = xmin,
                  SIMPLIFY = TRUE)
  xproj <- data.table(xproj)
  colnames(xproj) <- var_names
  attr(xproj, "vine_rc") <- attr(rcu, "vine")
  attr(xproj, "kde_rc") <- attr(rcu, "kde")
  attr(xproj, "vine_mp") <- attr(mpu, "vine")
  attr(xproj, "kde_mp") <- attr(mpu, "kde")
  class(xproj) <- c("vbc", class(xproj))
  if(any(!is.na(time_mp))) {
    xproj[, "time" := time_mp]
  }
  xproj
}

# utils ------------------------------------------------------------------------

#' @title Check arguments for vine correction
#' @inheritParams vbc
#' @import checkmate
check_vbc_args <- function(mp, mc, rc, var_names) {
  if(is.data.frame(mp)) {
    assert_data_frame(mp, types = "numeric", any.missing = FALSE,
                      ncols = ncol(mc))
  } else {
    lapply(mp, function(mem){
      assert_data_frame(mem, types = "numeric", any.missing = FALSE,
                        ncols = ncol(mc))
    })
  }
  assert_data_frame(mc, types = "numeric", any.missing = FALSE,
                    ncols = ncol(rc))
  assert_data_frame(rc, types = "numeric", any.missing = FALSE,
                    ncols = ncol(mp))
  assert_set_equal(apply(rc, 2, class), apply(mc, 2, class), ordered = TRUE)
  assert_character(var_names, any.missing = FALSE, unique = TRUE,
                   len = ncol(mp))
}

#' @title Correction by VBC for model ensembles
#' 
#' @inheritParams vbc
#' 
#' @return [list]\cr
#' A list of corrected data frames for each member of the model ensemble.
#' 
vbc_ensemble <- function(mp, mc, rc, var_names, margins_controls, time_mp, ...) {
  mc_kde <- attr(estimate_margins(mc, margins_controls), "kde")
  rcu <- model_vine(rc, margins_controls, ...)
  attr(rcu, "vine")$var_types = rep("c", times = ncol(rc))
  lapply(mp, function(member){
    mpu <- model_vine(member, margins_controls, ...)
    x_mph <- correct_rosenblatt(mpu, rcu)
    xmin = if(length(margins_controls$xmin) != ncol(rc)) {
      rep(NA, times = ncol(rc))
    } else {
      margins_controls$xmin
    }
    xproj <- mapply(map_delta, mp = member, mph = data.frame(x_mph),
                    mp_kde = attr(mpu, "kde"), mc_kde = mc_kde, xmin = xmin,
                    SIMPLIFY = TRUE)
    xproj <- data.table(xproj)
    colnames(xproj) <- var_names
    attr(xproj, "vine_rc") <- attr(rcu, "vine")
    attr(xproj, "kde_rc") <- attr(rcu, "kde")
    attr(xproj, "vine_mp") <- attr(mpu, "vine")
    attr(xproj, "kde_mp") <- attr(mpu, "kde")
    class(xproj) <- c("vbc", class(xproj))
    if(!is.na(time_mp)) {
      xproj[, "time" := time_mp]
    }
    message("An ensemble member is done.")
    xproj
  })  
}
