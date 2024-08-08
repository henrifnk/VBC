#' Correct data by Rosenblatt transformation
#' 
#' @param mpu [data.table]\cr
#' Model data from [vine_correct()].
#' 
#' @param ocu [data.table]\cr
#' Observed data from [vine_correct()].
#' 
#' @param z_inf `logical(1)`\cr
#' If `TRUE` at least one margin is zero inflated.
#' 
#' @return Corrected data by Rosenblatt transformation.
correct_rosenblatt <- function(mpu, ocu, z_inf = FALSE) {
  mp_vine <- attr(mpu, "vine")
  oc_vine <- attr(ocu, "vine")
  oc_kde <- attr(ocu, "kde")
  if(z_inf) {
    mpu_m <- as.matrix(mpu)
    u = rosenblatt_discrete(mpu_m, mp_vine)
  } else {
    u = rosenblatt(mpu, mp_vine)
  }
  u <- pseudo_obs(u, ties_method = 'average')
  u_mph = inverse_rosenblatt(u, oc_vine)
  u_mph <- pseudo_obs(u_mph, ties_method = 'average')
  x_mph <- mapply(function(u, kde) {
    qkde1d(u, kde)
  }, u = data.table(u_mph), kde = oc_kde)
  x_mph
}