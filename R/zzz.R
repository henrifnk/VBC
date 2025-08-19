#' @importFrom utils getFromNamespace
.rvine_internal <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  getNS <- utils::getFromNamespace
  .rvine_internal$vinecop_rosenblatt_cpp <- getNS("vinecop_rosenblatt_cpp",
                                                  "rvinecopulib")
  .rvine_internal$vinecop_inverse_rosenblatt_cpp <- 
    getNS("vinecop_inverse_rosenblatt_cpp",
          "rvinecopulib")
  .rvine_internal$dpq_marg <- getNS("dpq_marg", "rvinecopulib")
  .rvine_internal$get_seeds <- getNS("get_seeds", "rvinecopulib")
  .rvine_internal$if_vec_to_matrix <- getNS("if_vec_to_matrix", "rvinecopulib")
}
