#' Rosenblatt transformation
#' @description
#' The Rosenblatt transformation is a method to transform uniform 
#' pseudo-observations into the original scale of the data, preserving the
#' copula structure. Adapted from [rvinecopulib::rosenblatt()] and
#' [rvinecopulib::inverse_rosenblatt()].
#' 
#' @param x A numeric matrix or data.frame of observations. If not on \eqn{[0,1]^d},
#'   they will be converted via \code{rvinecopulib::pseudo_obs()} unless already
#'   provided as pseudo-observations.
#' @param model A \code{vinecop} object (from \pkg{rvinecopulib}) describing the
#'   dependence structure.
#' @param randomize_discrete Logical; if \code{TRUE}, apply jittering/randomized
#'   PIT for discrete/mixture margins (Rosenblatt for non-continuous margins).
#'   Default: \code{FALSE}.
#' @param cores Integer; number of threads to use. Default: \code{1}.
#' @return A numeric matrix of the same dimension as \code{x}, containing
#'   PIT-transformed values in \eqn{[0,1]}.
#'
#' @importFrom assertthat assert_that is.number
rosenblatt <- function(x, model, cores = 1, randomize_discrete = TRUE) {
  assert_that(
    inherits(model, c("bicop_dist", "vinecop_dist", "vine_dist")),
    is.number(cores)
  )
  
  if (inherits(model, "bicop_dist")) {
    model <- vinecop_dist(
      list(list(model)),
      cvine_structure(1:2),
      var_types = model$var_types
    )
  }
  
  if (inherits(model, "vine_dist")) {
    x <- .rvine_internal$expand_factors(x)
    if (!is.null(model$names)) {
      x <- x[, model$names, drop = FALSE]
    }
    x <- .rvine_internal$compute_pseudo_obs(x, model)
    model <- model$copula
  }
  
  # model is now a vinecop_dist
  assert_that(all((x >= 0) & (x <= 1)))
  x <- as.matrix(x)
  x <- pmin(pmax(x, 1e-10), 1 - 1e-10)
  x <- .rvine_internal$if_vec_to_matrix(x, dim(model)[1] == 1)
  x <- .rvine_internal$vinecop_rosenblatt_cpp(x, model, cores, randomize_discrete, .rvine_internal$get_seeds())
  colnames(x) <- unique(model$names)
  
  x
}

inverse_rosenblatt <- function(u, model, cores = 1) {
  assert_that(
    all((u > 0) & (u < 1)),
    inherits(model, c("bicop_dist", "vinecop_dist", "vine_dist")),
    is.number(cores)
  )
  
  to_col <- if (inherits(model, "bicop_dist")) FALSE else (dim(model)[1] == 1)
  u <- .rvine_internal$if_vec_to_matrix(u, to_col)
  
  if (inherits(model, "bicop_dist")) {
    model <- vinecop_dist(
      list(list(model)),
      cvine_structure(1:2),
      var_types = model$var_types
    )
  }
  
  if (inherits(model, "vinecop_dist")) {
    u <- .rvine_internal$vinecop_inverse_rosenblatt_cpp(u, model, cores)
  } else {
    u <- .rvine_internal$vinecop_inverse_rosenblatt_cpp(u, model$copula, cores)
    u <- .rvine_internal$dpq_marg(u, model, "q")
  }
  colnames(u) <- unique(model$names)
  
  u
}
