#' @title Model Multivariate Distribution by Vine Copula
#' 
#' @description
#' Estimate the multivariate distribution of the model data via vine copula
#' estimation (see [rvinecopulib::vine]).
#' 
#' @param data [data.frame]
#' Data to estimate the multivariate distribution.
#' 
#' @inheritParams vbc
#' 
#' @return The PIT-transformed margins from [estimate_margins()]. Additionally
#' the data frame contains the attribute `vine` with the vine copula model and
#' the attribute `kde` with the kernel density estimation of the data.
model_vine <- function(data, margins_controls, ...) {
  u_data <- estimate_margins(data, margins_controls)
  if (any(margins_controls$type == "zi")) {
    vec <- which(rep(margins_controls$type, times = 2) == "zi")
    u_data[, -vec] <- pseudo_obs(u_data[, -vec], ties_method = 'random')
  } else {
    u_data <- pseudo_obs(u_data, ties_method = 'random')
  }
  var_types <- ifelse(margins_controls$type == "zi", "d", "c")
  vine_model <- vinecop(u_data, var_types = var_types, ...)
  attr(u_data, "vine") <- vine_model
  u_data
}