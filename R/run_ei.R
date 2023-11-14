ei_estimate <- function(row_margins, col_margins, model = "contextual", iter = 2000,
                        verbose = TRUE){

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
                         iter = iter)
}
