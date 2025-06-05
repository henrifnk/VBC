
#' 
#' 
#' 
rosenblatt <- function(x, model, cores = 1, randomize_discrete = TRUE) {
  assertthat::assert_that(
    inherits(model, c("bicop_dist", "vinecop_dist", "vine_dist")),
    assertthat::is.number(cores)
  )
  
  if (inherits(model, "bicop_dist")) {
    model <- vinecop_dist(
      list(list(model)),
      cvine_structure(1:2),
      var_types = model$var_types
    )
  }
  
  if (inherits(model, "vine_dist")) {
    x <- expand_factors(x)
    if (!is.null(model$names)) {
      x <- x[, model$names, drop = FALSE]
    }
    x <- compute_pseudo_obs(x, model)
    model <- model$copula
  }
  
  # model is now a vinecop_dist
  assertthat::assert_that(all((x >= 0) & (x <= 1)))
  x <- as.matrix(x)
  x <- pmin(pmax(x, 1e-10), 1 - 1e-10)
  x <- rvinecopulib:::if_vec_to_matrix(x, dim(model)[1] == 1)
  x <- rvinecopulib:::vinecop_rosenblatt_cpp(x, model, cores, randomize_discrete, rvinecopulib:::get_seeds())
  colnames(x) <- unique(model$names)
  
  x
}

inverse_rosenblatt <- function(u, model, cores = 1) {
  assertthat::assert_that(
    all((u > 0) & (u < 1)),
    inherits(model, c("bicop_dist", "vinecop_dist", "vine_dist")),
    assertthat::is.number(cores)
  )
  
  to_col <- if (inherits(model, "bicop_dist")) FALSE else (dim(model)[1] == 1)
  u <- rvinecopulib:::if_vec_to_matrix(u, to_col)
  
  if (inherits(model, "bicop_dist")) {
    model <- vinecop_dist(
      list(list(model)),
      cvine_structure(1:2),
      var_types = model$var_types
    )
  }
  
  if (inherits(model, "vinecop_dist")) {
    u <- rvinecopulib:::vinecop_inverse_rosenblatt_cpp(u, model, cores)
  } else {
    u <- rvinecopulib:::vinecop_inverse_rosenblatt_cpp(u, model$copula, cores)
    u <- rvinecopulib:::dpq_marg(u, model, "q")
  }
  colnames(u) <- unique(model$names)
  
  u
}
