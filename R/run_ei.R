ei_estimate <- function(row_margins, col_margins, model){

  R <- ncol(row_margins)
  C <- ncol(col_margins)

  out <- rstan::sampling(stanmodel$ssMND,
                         data = standata)
}
