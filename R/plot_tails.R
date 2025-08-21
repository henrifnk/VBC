#' @title Plot Tails
#' @description This function plots the tails probabilities of a marginal
#' distribution.
#' @param data [data.table::data.table]\cr
#' A data.table containing the data to be plotted.
#' @param var [character]\cr
#' The variable to be plotted. Must be a column of `data`.
#' @param xmin [numeric]\cr
#' The minimum value of the variable to be plotted. Passed to [kde1d::kde1d].
#' Default is `NA` indicating no lower limit.
#' @param scale_d [numeric]\cr
#' The scaling factor for the density plot between continuous and relative axis.
#' Default is `3`. Regulate the factor, if continuous density estimates and 
#' relative frequency are too distinct on the y-axis. Only used for zero
#' inflated variables.
#' @param mult [numeric]\cr
#' Bin width multiplier. Passed to [kde1d::kde1d].
#' Default is 1 indicating no multiplication.
#' @param length_out [numeric]\cr
#' The length of the sequence to be used for the density plot. Default is `1e4`.
#' @return A ggplot that plots the continuous density estimate of the variable 
#' and the relative frequency of inflation.
#' @importFrom scales gradient_n_pal trans_new
#' @import ggplot2
#' @export
plot_tails = function(data, var, xmin = NA, scale_d = 3, mult = 1,
                      length_out = 1e4) {
  is_inf = sum(data[, get(var) == 0])/ nrow(data) > 0.01
  zero_data = data.table("zero_freq" = data[, sum(get(var) == 0) / .N])
  zero_data[, "x" := 0]
  if(is_inf) {
    xmin = 0
    seq = seq(from = log(min(data[, get(var)][data[, get(var)] != 0])),
              to = log(max(data[, get(var)])), length.out = length_out)
    seq = exp(seq)
  } else {
    seq = seq(from = min(data[, get(var)]), to = max(data[, get(var)]),
              length.out = length_out)
  }
  density_plot <- data[get(var) != 0]
  kde <- kde1d(data[, get(var)], xmin = xmin, type = ifelse(is_inf, "zi", "c"),
               mult = mult)
  kde_data <- data.table(x = seq)
  kde_data[, `:=`("density" = dkde1d(kde_data$x, kde),
                    "cumulative_density" = pkde1d(kde_data$x, kde)
  )]
  kde_data <- kde_data[, "tail_probs" := 0.5 - abs(
    0.5 - get("cumulative_density"))]
  kde_data[, "scaled_density" := get("density") / max(get("density")) *
             min(get("cumulative_density")) * scale_d]
  colour_bar = gradient_n_pal(
    c("blue", "yellow"), c(0, 0.5)
    )(c(0.5 - abs(0.5 - zero_data$zero_freq)))
  bar <- if(is_inf) {
    geom_bar(data = zero_data, aes(x = get("x"), y = get("zero_freq")),
             fill = colour_bar, stat = "identity", show.legend = FALSE,
             just = 1, width = 0.2)
  } else NULL
  log_trafo <- if(is_inf) {
    scale_x_continuous(trans = trans_new(name = "log_plus_one",
                       transform = function(x) log(x + 1),
                       inverse = function(x) exp(x)- 1))
    } else NULL
  sec_axis <- if(is_inf) {
    scale_y_continuous(
      name = "relative frequency of inflation",,
      sec.axis = sec_axis(~ . / (scale_d), name = "continous density estimate")
    )
  } else ylab("continous density estimate")
  ggplot(data = kde_data) + bar +
    geom_segment(aes(x = get("x"), y = 0, xend = get("x"),
                     yend = scale_d * get("density"), 
                     color = get("tail_probs")), linewidth = 0.5) + sec_axis +
    scale_color_gradient(low = "blue", high = "yellow", name = "Tail Prob.",
                         limits = c(0, 0.5)) +
    log_trafo +
    labs(x = var) + 
    theme(strip.text.y  = element_text(face = "bold"),
          panel.grid.minor = element_blank())
}
