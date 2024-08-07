# utils ------------------------------------------------------------------------
#' Calculate Margins
#'
#' @inheritParams vine_correct
#' @param dta [data.table]\cr
#' Data to calculate the margins.
#' @param kde [list]\cr
#' kernel density estimations.
#'
#' @return A data frame with the marginal PIT.
calculate_margins = function(dta, margins_controls, kde = NULL) {
  if (is.null(kde)) {
    margins_controls <- expand_margin_controls(margins_controls, dta)
    kde <- mapply(kde1d, x = dta, mult = margins_controls$mult,
                  xmin = margins_controls$xmin, xmax = margins_controls$xmax,
                  bw = margins_controls$bw, deg = margins_controls$deg,
                  type = margins_controls$type, SIMPLIFY = FALSE)
  }
  u = mapply(pkde1d, q = data.frame(dta), obj = kde)
  if(any(margins_controls$type == "zi")) {
    u_sub <- u
    u_sub[, which(margins_controls$type == "zi")] <- do.call(
      "cbind", lapply(which(margins_controls$type == "zi"), function(zero_v) {
        u_sub[dta[, ..zero_v] == 0, zero_v] <- pkde1d(-1e300, kde[[zero_v]])
        u_sub[, zero_v]
      }))
    u_final = cbind.data.frame(u, u_sub)
  } else {
    u_final = u
  }
  attr(u_final, "kde") <- kde
  u_final
}

#' Expand margins
#'
#' @inheritParams calculate_margins
#' @param controls margins_controls
#' @importFrom utils modifyList
#' @return Expanded margins.
expand_margin_controls = function (controls, dta) {
  default_controls <- list(mult = NULL, xmin = NaN, xmax = NaN,
                           bw = NA, deg = 2)
  controls <- modifyList(default_controls, controls)
  if (is.null(controls[["mult"]])) controls[["mult"]] <- log(1 + ncol(dta))
  for (par in names(controls)) {
    if (length(controls[[par]]) != ncol(dta))
      controls[[par]] <- rep(controls[[par]], ncol(dta))
  }
  controls
}
