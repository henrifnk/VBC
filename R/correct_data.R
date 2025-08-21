#' Correct data by Rosenblatt transformation
#' 
#' @param mpu [data.table::data.table]\cr
#' Model data from [vbc()].
#' 
#' @param rcu [data.table::data.table]\cr
#' Observed data from [vbc()].
#' 
#' @return Corrected data by Rosenblatt transformation.
correct_rosenblatt <- function(mpu, rcu) {
  mp_vine <- attr(mpu, "vine")
  len <- length(mp_vine$var_types)
  rc_vine <- attr(rcu, "vine")
  rc_kde <- attr(rcu, "kde")
  u <- rosenblatt(mpu, mp_vine)
  u <- pseudo_obs(u, ties_method = 'average')
  u_mph <- inverse_rosenblatt(u, rc_vine)
  u_mph <- pseudo_obs(u_mph, ties_method = 'average')
  x_mph <- mapply(function(u, kde) {
    qkde1d(u, kde)
  }, u = data.table(u_mph), kde = rc_kde)
  x_mph
}