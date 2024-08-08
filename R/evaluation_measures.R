#' @title Model Correction Inconcistency (MCI)
#'
#' @description
#' MCI provides a score to quantify how invasive the correction method
#' transformed the multivariate temporal rank structure of the model after
#' correcting the model bias. Therefore, this model compares non-excedance
#' probabilities of an event t before and after bias correction. A noninvaseive
#' method preserves the model-specific rank structure.
#'
#' @param model [data.frame]\cr
#' Simulation climate data from the climate model.
#'
#' @param model_correction [data.frame]\cr
#'  Simulation climate data `model` from the climate model corrected by an
#'  arbitrary bias correction method.
#'
#' @param time_p [date]\cr
#'  An optional vector with a time point for each observation/row in the
#'  (corrected) model data. Default is NA which indicates that no temporal
#'  information are given.
#'
#' @param p `numeric(1)`\cr
#'  The power of the global MCI. Default is 2.
#'
#' @inheritParams vine_correct
#'
#' @return [data.table]\cr
#' Contains the points of time if supplied in `time_p` and the respective
#' local MCI. Note that the last time step is truncated while derivation. The
#' attribute global contains the global MCI, which is the mean over all MCI
#' scores. The estimated multivariate distribution of model and correction are
#' placed in attributes in of the returned object.
#'
#' @example R/example.R
#'
#' @export
calc_mci <- function(model, model_correction, time_p = NA, p = 2,
                     margins_controls = list(mult = NULL, xmin = NaN,
                                             xmax = NaN, bw = NA, deg = 2,
                                             type = "c"), ...) {
  # Checks arguments
  assert_set_equal(dim(model), dim(model_correction), ordered = TRUE)
  assert_set_equal(colnames(model), colnames(model_correction), ordered = TRUE)
  if(!anyNA(time_p)) {
    assert_vector(time_p, len = nrow(model), any.missing = FALSE)
  } else time_p <- 1:nrow(model)
  pcop_model = calculate_pcop(model, margins_controls, ...)
  pcop_correction = calculate_pcop(model_correction, margins_controls, ...)
  probs = data.table(pcop_model, pcop_correction, time_p)
  colnames(probs) = c("nep_model", "nep_correction", "time")
  probs[, "mci" := get("nep_model") - get("nep_correction")]
  attr(probs, "global_mci") = mean(abs(probs$mci)^p)^1/p
  attr(probs, "kde_model") = attr(pcop_model, "kde")
  attr(probs, "kde_correction") = attr(pcop_correction, "kde")
  attr(probs, "vine_model") = attr(pcop_model, "vine")
  attr(probs, "vine_correction") = attr(pcop_correction, "vine")
  probs
}

#' @title Calculate Non-Exceedance Probability
#'
#' @inheritParams model_vine
#' 
#' @return A data frame with the non-exceedance probability for the events.
calculate_pcop <- function(data, margins_controls, ...) {
  u_data = model_vine(data, margins_controls, ...)
  z = pvinecop(u_data, attr(u_data, "vine"))
  attr(z, "kde") = attr(u_data, "kde")
  attr(z, "vine") = attr(u_data, "vine")
  z
}

#' @title Calculate Wasserstein Distance
#'
#' @description Calculate the Wasserstein distance of order 1 and 2 between
#' multivariate distributions of model and observed data. Can be used to
#' evaluate the performance of a bias correction method by comparing the
#' distribution of the model data after the correction.
#'
#' @inheritParams calc_mci
#'
#' @param observed [data.frame]\cr
#'  Measured (and interpolated) observations.
#'
#' @param n `integer(1)`\cr
#'  Length of samples to be drawn from model and observed data.
#'
#' @param scale_dta `logical(1)`\cr
#'  If `TRUE` the data is scaled before the calculation of the Wasserstein.
#'  Default is `TRUE`.
#'
#' @param ... \cr
#'  Further arguments passed to [transport::transport].
#'
#' @return A named vector of length 2.\cr
#' - Wasserstein_1: The Wasserstein distance of order 1
#' - Wasserstein_2: The Wasserstein distance of order 2
#'
#' @example R/example.R
#'
#' @seealso [transport::wasserstein()]
#'
#' @importFrom transport pp transport.pp wasserstein
#'
#' @export
calc_wasserstein <- function(observed, model, n = nrow(observed),
                             scale_dta = TRUE, ...) {
  assert_set_equal(ncol(observed), ncol(model))
  assert_set_equal(colnames(model), colnames(observed), ordered = TRUE)
  assert_integerish(n, lower = 1, len = 1, any.missing = FALSE)
  if(scale_dta) observed <- scale(observed); model = scale(model)
  samle_n <- min(nrow(observed), nrow(model), n)
  pp_o <- pp(observed[sample(nrow(observed), samle_n),])
  pp_m <- pp(model[sample(nrow(model), samle_n),])
  tplan_mo <- transport.pp(pp_o, pp_m, ...)
  c(
    "Wasserstein_1" = wasserstein(pp_o, pp_m, p = 1, tplan = tplan_mo),
    "Wasserstein_2" = wasserstein(pp_o, pp_m, p = 2, tplan = tplan_mo)
  )
}
