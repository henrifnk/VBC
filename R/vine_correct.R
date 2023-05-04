vine_correct = function(oc, mc, mp, var_names = colnames(oc), digits = NA,
                        zero_inf = rep(NA, times = nrow(oc)),
                        margins_controls = list(mult = NULL, xmin = NaN, 
                                                xmax = NaN, bw = NA, deg = 2), ...) {
  # jitter zero inf for numerical stability
  ocj = mapply(function(o, inf){
    if(is.na(inf)) return(o)
    o[o == 0] = rexp(n = sum(o == 0), rate = inf)
    o
  }, o = data.frame(oc), inf = zero_inf)
  
  mpj = mapply(function(m, inf){
    if(is.na(inf)) return(m)
    m[m == 0] = rexp(n = sum(m == 0), rate = inf)
    m
  }, m = data.frame(mp), inf = zero_inf)
  
  ## VBC step
  # oc data copula estimation
  vine_oc = vine(ocj, margins_controls = margins_controls, ...)
  vine_mp = vine(mpj, margins_controls = margins_controls, ...)
  
  # rosenblatt transform
  u_mph = rosenblatt(mp, vine_mp)
  x_mph = inverse_rosenblatt(u_mph, vine_oc) 
  
  ## Delta step and random truncation
  final = mapply(function(x, y, xh, inf){
    x_q = unname(quantile(x, probs = seq(0, 1, length.out = length(y))))
    y_sort = sort(y, index.return = TRUE)
    delta = (y_sort$x - x_q)[order(y_sort$ix)]
    xhd = xh + delta
    if(is.na(inf)) return(xhd)
    
    sz = min(round(sum(oc == 0) / nrow(oc) * nrow(mp), 0), sum(xhd < (1/inf)))
    sz = sz - sum(xhd < 0)
    id = which(xhd < 0)
    id = c(id, sample(which(xhd < (1 / inf)), size = sz))
    xhd[id] = 0
    xhd
  }, x = data.frame(mc), y = data.frame(mp), xh = data.frame(x_mph),
  inf = zero_inf, SIMPLIFY = TRUE)
  
  if(!is.na(digits)) final = round(final, digits)
  
  # Return Shaping
  colnames(final) = var_names
  attr(final, "calibration") = x_mph
  attr(final, "vine_oc") = vine_oc
  attr(final, "vine_mp") = vine_mp
  #attr(final, "hyman") = models
  final
}
