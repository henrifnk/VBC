#' Correct data by Rosenblatt transformation
#' 
#' @param mpu [data.table]\cr
#' Model data from [vbc()].
#' 
#' @param rcu [data.table]\cr
#' Observed data from [vbc()].
#' 
#' @param z_inf `logical(1)`\cr
#' If `TRUE` at least one margin is zero inflated.
#' 
#' @return Corrected data by Rosenblatt transformation.
correct_rosenblatt <- function(mpu, rcu, z_inf = FALSE) {
  mp_vine <- attr(mpu, "vine")
  rc_vine <- attr(rcu, "vine")
  rc_kde <- attr(rcu, "kde")
  if(z_inf) {
    mpu_m <- as.matrix(mpu)
    u = rosenblatt_discrete(mpu_m, mp_vine)
  } else {
    u = rosenblatt(mpu, mp_vine)
  }
  u <- pseudo_obs(u, ties_method = 'average')
  u_mph = inverse_rosenblatt(u, rc_vine)
  u_mph <- pseudo_obs(u_mph, ties_method = 'average')
  x_mph <- mapply(function(u, kde) {
    qkde1d(u, kde)
  }, u = data.table(u_mph), kde = rc_kde)
  x_mph
}