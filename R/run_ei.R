#' Run ssEI ecological inference model.
#'
#' @param row_margins A data frame with the row margins data, with one row for each area, and with columns giving the row margins in the area.
#' @param col_margins A data frame with the column margins data, with one row for each area, and with columns giving the column margins in the area.
#' @param model Which ssEI model to run (either contextual or average)
#' @param verbose Print some information about the data and model before running
#' @param ... other arguments to be passed to rstan::sampling
#'
#' @return
#' @export
#'
#' @examples
ei_estimate <- function(row_margins, col_margins, model = "contextual",
                        verbose = TRUE, ...){

  R <- ncol(row_margins)
  C <- ncol(col_margins)
  n_areas <- nrow(row_margins)
  if(verbose){
    print(R)
    print(C)
    print(n_areas)

  }

  standata <- list(
    n_areas = n_areas,
    R = R,
    C = C,
    row_margins = row_margins,
    col_margins = col_margins
  )

  if(model=="contextual"){
    mod <- stanmodels$ssContextual
  } else{
    mod <- stanmodels$ssMND
  }
  if(verbose){
    print(paste("now running model", mod@model_name))
  }

  out <- rstan::sampling(mod,
                         data = standata,
                         cores = 4,
                         ...)
}
