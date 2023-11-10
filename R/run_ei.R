ei_estimate <- function(row_margins, col_margins, model, iter = 2000){

  R <- ncol(row_margins)
  C <- ncol(col_margins)
  n_areas <- nrow(row_margins)

  standata <- list(
    n_areas = n_areas,
    R = R,
    C = C,
    row_margins = row_margins,
    col_margins = col_margins
  )

  out <- rstan::sampling(stanmodels$ssMND,
                         data = standata,
                         cores = 4,
                         iter = iter)
}
