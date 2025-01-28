
#' 
rosenblatt <- function (x, model, cores = 1, randomize_discrete = TRUE) 
{
  assertthat::assert_that(inherits(model, c("bicop_dist", "vinecop_dist", 
                                "vine_dist")), assertthat::is.number(cores))
  to_col <- if (inherits(model, "bicop_dist")) 
    FALSE
  else (dim(model)[1] == 1)
  x <- rvinecopulib:::if_vec_to_matrix(x, to_col)
  col_names <- colnames(x)
  if (inherits(model, "bicop_dist")) {
    model <- vinecop_dist(list(list(model)), cvine_structure(1:2), 
                          var_types = model$var_types)
  }
  if (inherits(model, "vinecop_dist")) {
    assertthat::assert_that(all((x >= 0) & (x <= 1)))
    x <- pmin(pmax(x, 1e-10), 1 - 1e-10)
    x <- rvinecopulib:::vinecop_rosenblatt_cpp(x, model, cores, 
                                               randomize_discrete, 
                                               rvinecopulib:::get_seeds())
  }
  else {
    stop("Model must be of class bicop_dist or vinecop_dist")
  }
  colnames(x) <- col_names[1:ncol(x)]
  x
}
