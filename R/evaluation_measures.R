#' @title Model Correction Inconsistency (MCI)
#'
#' @description
#' MCI provides a score to quantify how invasive the correction method
#' transformed the multivariate temporal rank structure of the model after
#' correcting the model bias. Therefore, this model compares non-exceedance
#' probabilities of an event t before and after bias correction. A non-invasive
#' method preserves the model-specific rank structure.
#'
#' @param model [data.frame] 
#' Simulation climate data from the climate model.
#'
#' @param model_correction [data.frame] 
#'  Simulation climate data `model` from the climate model corrected by an
#'  arbitrary bias correction method.
#'
#' @param time_p [date] 
#'  An optional vector with a time point for each observation/row in the
#'  (corrected) model data. Default is NA which indicates that no temporal
#'  information are given.
#'
#' @param p `numeric(1)` 
#'  The power of the global MCI. Default is 1.
#'
#' @return [data.table::data.table] 
#' Contains the points of time if supplied in `time_p` and the respective
#' local MCI. Note that the last time step is truncated while derivation. The
#' attribute global contains the global MCI, which is the mean over all MCI
#' scores. The estimated multivariate distribution of model and correction are
#' placed in attributes in of the returned object.
#'
#' @example R/example.R
#'
#' @export
calc_mci <- function(model, model_correction, time_p = NA, p = 1) {
  assert_set_equal(dim(model), dim(model_correction), ordered = TRUE)
  assert_set_equal(colnames(model), colnames(model_correction), ordered = TRUE)
  if(!anyNA(time_p)) {
    assert_vector(time_p, len = nrow(model), any.missing = FALSE)
  } else time_p <- 1:nrow(model)
  pcop_model = calculate_pemp(model)
  pcop_correction = calculate_pemp(model_correction)
  probs = data.table(pcop_model, pcop_correction, time_p)
  colnames(probs) = c("nep_model", "nep_correction", "time")
  probs[, "mci" := abs(get("nep_model") - get("nep_correction"))]
  attr(probs, "global_mci") = mean(probs$mci^p)^1/p
  probs
}

#' @title Calculate Non-Exceedance Probability
#'
#' @param data [data.table::data.table] 
#' On each of the rows in the data.table, the non-exceedance probability is
#' calculated. 
#'  
#' @return A data frame with the non-exceedance probability for the events. The 
#' non-exceedance probability is the percentage of rows that satisfy the
#' condition that all values in the row are less or equal to the values in the
#' row.
calculate_pemp <- function(data) {
  apply(data, 1, function(row) {
    row_satisfies_condition = apply(data, 1, function(x) all(x <= row))
    # Calculate the percentage of rows that satisfy the condition
    mean(row_satisfies_condition)
  })
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
#' @param observed [data.frame] 
#'  Measured (and interpolated) observations.
#'
#' @param n `integer(1)` 
#'  Length of samples to be drawn from model and observed data.
#'
#' @param scale_dta `cahracter(1)` 
#' Scale the data before calculating the Wasserstein distance. Default is
#' "observed". Possible values are "both", "observed", and "sd". If "both" the
#' data are scaled to zero mean and unit variance. If "observed" the data are
#' scaled to zero mean and unit variance of the observed data. If "sd" the data
#' are scaled to zero mean and unit variance.
#'
#' @param ...  
#'  Further arguments passed to [transport::transport].
#'
#' @return A named vector of length 2. 
#' - Wasserstein_1: The Wasserstein distance of order 1
#' - Wasserstein_2: The Wasserstein distance of order 2
#'
#' @example R/example.R
#'
#' @seealso [transport::wasserstein()]
#'
#' @importFrom transport pp transport.pp wasserstein wasserstein1d
#'
#' @export
calc_wasserstein <- function(observed, model, n = nrow(observed),
                             scale_dta = "observed", ...) {
  assert_set_equal(ncol(observed), ncol(model))
  assert_set_equal(colnames(model), colnames(observed), ordered = TRUE)
  assert_integerish(n, lower = 1, len = 1, any.missing = FALSE)
  samle_n <- min(nrow(model), n)
  if(ncol(observed) == 1) {
    return(
      c("Wasserstein_1" = wasserstein1d(unlist(observed), unlist(model), p = 1),
        "Wasserstein_2" = wasserstein1d(unlist(observed), unlist(model), p = 2))
    )
  }
  if(scale_dta == "both") {
    observed <- scale(observed)
    model <- scale(model)
  } else if(scale_dta == "observed") {
    observed <- scale(observed)
    model <- scale(model, center = attr(observed, "scaled:center"),
                   scale = attr(observed, "scaled:scale"))
  } else if(scale_dta == "sd") {
    model <- scale(model, center = FALSE)
    observed <- scale(observed, center = FALSE)
  }

  pp_o <- pp(observed[sample(nrow(observed), samle_n),])
  pp_m <- pp(model[sample(nrow(model), samle_n),])
  tplan_mo <- transport.pp(pp_o, pp_m, ...)
  c(
    "Wasserstein_1" = wasserstein(pp_o, pp_m, p = 1, tplan = tplan_mo),
    "Wasserstein_2" = wasserstein(pp_o, pp_m, p = 2, tplan = tplan_mo)
  )
}
