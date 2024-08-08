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
#' @param oc [data.table]\cr
#'  Measured (and interpolated) observations during calibration period.
#'
#' @param mc [data.table]\cr
#'  Simulation data from a climate model during calibration period.
#'
#' @param mp [data.table]\cr
#'  Simulation data from a climate model during projection period.
#'
#' @param var_names [character]\cr
#'  Names of corrected climate data. Defaults to the column names of the
#'  observed data `oc`.
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
#' @param ... \cr
#'  Arguments are passed to [rvinecopulib::vinecop] to specify the structure of
#'  vines and margins. Note that the ellipsis of observed and model data are
#'  specified with the same arguments.
#'
#' @return [data.table]\cr
#'  The corrected projection period data in `mp`. Additionally the data frame
#'  contains the attributes
#'
#'  - `vine_oc` the vine copula for the measured data `oc`.
#'  - `vine_mp` the vine copula for the model data from projection period `mp`.
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
vbc <- function(oc, mc, mp, var_names = colnames(oc), margins_controls = list(
  mult = NULL, xmin = NaN, xmax = NaN, bw = NA, deg = 2, type = "c"
), ...) {
  check_vbc_args(oc, mc, mp, var_names)
  mc_kde <- attr(estimate_margins(mc, margins_controls), "kde")
  mpu <- model_vine(mp, margins_controls, ...)
  ocu <- model_vine(oc, margins_controls, ...)
  attr(ocu, "vine")$var_types = rep("c", times = ncol(mp))
  x_mph <- correct_rosenblatt(mpu, ocu, any(margins_controls$type == "zi"))
  xmin = if(length(margins_controls$xmin) != ncol(oc)) {
    rep(NA, times = ncol(oc))
  } else {
    margins_controls$xmin
  }
  xproj <- mapply(map_delta, mp = mp, mph = data.frame(x_mph),
                  mp_kde = attr(mpu, "kde"), mc_kde = mc_kde, xmin = xmin,
                  SIMPLIFY = TRUE)
  xproj <- data.table(xproj)
  colnames(xproj) <- var_names
  attr(xproj, "vine_oc") <- attr(ocu, "vine")
  attr(xproj, "kde_oc") <- attr(ocu, "kde")
  attr(xproj, "vine_mp") <- attr(mpu, "vine")
  attr(xproj, "kde_mp") <- attr(mpu, "kde")
  class(xproj) <- c("vbc", class(xproj))
  xproj
}

# utils ------------------------------------------------------------------------

#' @title Check arguments for vine correction
#' @inheritParams vbc
#' @import checkmate
check_vbc_args <- function(oc, mc, mp, var_names) {
  if(is.list(oc)) {
    assert_data_frame(oc, types = "numeric", any.missing = FALSE,
                      ncols = ncol(mc))
  }
  assert_data_frame(mc, types = "numeric", any.missing = FALSE,
                    ncols = ncol(mp))
  assert_data_frame(mp, types = "numeric", any.missing = FALSE,
                    ncols = ncol(oc))
  assert_set_equal(apply(oc, 2, class), apply(mc, 2, class), ordered = TRUE)
  assert_set_equal(apply(mc, 2, class), apply(mp, 2, class), ordered = TRUE)
  assert_character(var_names, any.missing = FALSE, unique = TRUE,
                   len = ncol(mp))
}
