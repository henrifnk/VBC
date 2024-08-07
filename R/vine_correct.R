# VBC --------------------------------------------------------------------------
#' @title Correct zero inflated climate model simulations via vine copulas
#'
#' @description
#' The multivariate distribution functions of model and observed are estimated
#' via vine copula estimation (see [rvinecopulib::vine]).
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
vine_correct <- function(oc, mc, mp, var_names = colnames(oc),
                         margins_controls = list(mult = NULL, xmin = NaN,
                                                 xmax = NaN, bw = NA, deg = 2,
                                                 type = "c"), ...) {
  check_vbc_args(oc, mc, mp, var_names)
  ocu <- calculate_margins(oc, margins_controls)
  oc_kde <- attr(ocu, "kde")

  mpu <- calculate_margins(mp, margins_controls)
  mc_kde <- attr(calculate_margins(mc, margins_controls), "kde")
  if(any(margins_controls$type == "zi")) {
    vec = which(rep(margins_controls$type == "zi", times = 2))
    ocu[, -vec] = pseudo_obs(ocu[, -vec], ties_method = 'random')
    mpu[, -vec] = pseudo_obs(mpu[, -vec], ties_method = 'random')
  } else {
    mpu = pseudo_obs(mpu, ties_method = 'random')
    ocu = pseudo_obs(ocu, ties_method = 'random')
  }
  var_types = ifelse(margins_controls$type == "zi", "d", "c")
  vine_oc <- vinecop(ocu, var_types = var_types, ...)
  vine_oc$var_types = rep("c", times = ncol(mp))
  vine_mp <- vinecop(mpu, var_types = var_types, ...)
  if(any(margins_controls$type == "zi")) {
    mpu_m <- as.matrix(mpu)
    mpu_m <- flatten_zi_margins(mpu_m, margins_controls$type == "zi")
    u = rosenblatt_discrete(mpu_m, vine_mp)
  } else {
    u = rosenblatt(mpu, vine_mp)
  }
  u <- pseudo_obs(u, ties_method = 'average')
  u_mph = inverse_rosenblatt(u, vine_oc)
  u_mph <- pseudo_obs(u_mph, ties_method = 'average')
  x_mph <- mapply(function(u, kde) {
    qkde1d(u, kde)
  }, u = data.table(u_mph), kde = oc_kde)
  xmin = if(length(margins_controls$xmin) != ncol(oc)) {
    rep(NA, times = ncol(oc))
  } else {
    margins_controls$xmin
   }
  final <- mapply(map_delta, mp = mp, mph = data.frame(x_mph),
                  mp_kde = attr(mpu, "kde"), mc_kde = mc_kde, xmin = xmin,
                  SIMPLIFY = TRUE)
  final <- data.table(final)
  colnames(final) <- var_names
  attr(final, "vine_oc") <- vine_oc
  attr(final, "kde_oc") <- attr(ocu, "kde")
  attr(final, "vine_mp") <- vine_mp
  attr(final, "kde_mp") <- attr(mpu, "kde")
  final
}

# utils ------------------------------------------------------------------------

#' @title Check arguments for vine correction
#' @inheritParams vine_correct
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

#' @title Delta mapping
#' @description
#' Map the rank delta in model between calibration period `mc` and projection
#' period `mp` on corrected projection period `mph` to correct for climate
#' trends. For ratio scaled variables, multiplicative delta scheme is applied
#' while for other
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
#'  A vector indicating if xmin is a ratio variable type or any other type (NA).
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

#' @title Flatten zero inflated margins
#' @param margins [data.frame]\cr
#'  The margins of the climate data.
#' @param eps [double]\cr
#' A small number to avoid zero and one values in zero inflated margins.
#' @inheritParams vine_correct
#' @return Flatten zero inflated margins.
flatten_zi_margins <- function(margins, zero_inf, eps = 1e-10) {
  vec = which(rep(zero_inf, times = 2))
  margins_z = margins[, vec]
  margins_z[margins_z == 1] <- 1 - eps
  margins[, vec] = margins_z
  margins_z = margins[ , which(zero_inf)]
  margins_z[margins_z == 0] <- 0 + eps
  margins[ , which(zero_inf)] = margins_z
  margins
}
