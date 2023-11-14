
#' Obtain cell value estimates from ssEI ecological inference models.
#'
#' @param ei_stanfit Usually the result of running the `ei_estimate` command. A `stanfit` object that results from running the one of the ssEI models.
#' @param pars the name of the parameters to extract (usually `cell_values`)
#' @param probs probabilities to pass to rstan::summary.stanfit
#' @param ...
#'
#' @return A tibble with estimates of cell values for each combination of area number, row number and column number.
#' @export
#'
#' @examples
ei_cv_summary <- function(ei_stanfit, pars = "cell_values", probs = c(0.025, 0.25, 0.50, 0.75, 0.975), ...){

  cv_summary <- rstan::summary(ei_stanfit,
                               pars = pars,
                               probs = probs,
                               ...
  )
  cv_summary$summary |>
    as.matrix() |>
    as_tibble(rownames="param") |>
    separate("param",
             into=c("param_name", "area_no", "row_no", "col_no", NA),
             sep="[\\[,\\]]",
             remove=FALSE)
}


#' Obtain row rate estimates from ssEI ecological inference models.
#'
#' @param ei_stanfit Usually the result of running the `ei_estimate` command. A `stanfit` object that results from running the one of the ssEI models.
#' @param pars the name of the parameters to extract (usually `row_rate`)
#' @param probs probabilities to pass to rstan::summary.stanfit
#' @param ...
#'
#' @return A tibble with estimates of row rates (the percentage of the row total in each column) for each combination of area number, row number and column number. For rows which sum to zero (and hence consist only of zero cells) a value of -1 is returned.
#' @export
#'
#' @examples
ei_row_rate_summary <- function(ei_stanfit, pars = "row_rate", probs = c(0.025, 0.25, 0.50, 0.75, 0.975), ...){

  cv_summary <- rstan::summary(ei_stanfit,
                               pars = pars,
                               probs = probs,
                               ...
  )
  cv_summary$summary |>
    as.matrix() |>
    as_tibble(rownames="param") |>
    separate("param",
             into=c("param_name", "area_no", "row_no", "col_no", NA),
             sep="[\\[,\\]]",
             remove=FALSE)
}

#' Obtain overall row rate estimates (total across all areas) from an ssEI model.
#'
#' @param ei_stanfit Usually the result of running the `ei_estimate` command. A `stanfit` object that results from running the one of the ssEI models.
#' @param pars the name of the parameters to extract (usually `overall_row_rates`)
#' @param probs probabilities to pass to rstan::summary.stanfit
#' @param ...
#'
#' @return A tibble with estimates of row rates (the percentage of the row total in each column) for each combination of row number and column number (summary total across all areas). For rows which sum to zero (and hence consist only of zero cells) a value of -1 is returned.
#' @export
#'
#' @examples
ei_overall_row_rate_summary <- function(ei_stanfit, pars = "overall_row_rates", probs = c(0.025, 0.25, 0.50, 0.75, 0.975), ...){

  cv_summary <- rstan::summary(ei_stanfit,
                               pars = pars,
                               probs = probs,
                               ...
  )
  cv_summary$summary |>
    as.matrix() |>
    as_tibble(rownames="param") |>
    separate("param",
             into=c("param_name", "row_no", "col_no", NA),
             sep="[\\[,\\]]",
             remove=FALSE)
}

#' Obtain overall total cell values estimates (total across all areas) from an ssEI model.
#'
#' @param ei_stanfit Usually the result of running the `ei_estimate` command. A `stanfit` object that results from running the one of the ssEI models.
#' @param pars the name of the parameters to extract (usually `overall_row_rates`)
#' @param probs probabilities to pass to rstan::summary.stanfit
#' @param ...
#'
#' @return A tibble with estimates of total in each cell summing across all areas for each combination of row number and column number (summary total across all areas).
#' @export
#'
#' @examples
ei_overall_cell_values_summary <- function(ei_stanfit, pars = "overall_cell_values", probs = c(0.025, 0.25, 0.50, 0.75, 0.975), ...){

  cv_summary <- rstan::summary(ei_stanfit,
                               pars = pars,
                               probs = probs,
                               ...
  )
  cv_summary$summary |>
    as.matrix() |>
    as_tibble(rownames="param") |>
    separate("param",
             into=c("param_name", "row_no", "col_no", NA),
             sep="[\\[,\\]]",
             remove=FALSE)
}

ei_to_cvdraws <- function(ei_posterior){
  cv <- ei_posterior$cell_values
  cvdf <- as.data.frame(cv)
  cvdf <- cvdf |>
    dplyr::mutate(iter = dplyr::row_number())|>
    tidyr::pivot_longer(cols=-iter, values_to = "est_cell_value") |>
    tidyr::separate("name", into = c("area_no", "row_no", "col_no")) |>
    dplyr::mutate(area_no=as.numeric(area_no),
                  row_no = as.numeric(row_no),
                  col_no = as.numeric(col_no))
  return(cvdf)
}

cvdraws_to_rowprops <- function(cv_draws){
  #' @description
    #' Adds row proportion calculation to a cell values draws tibble.
    #'
    #' @param cv_draws a cv_draws tibble created by the `ei_to_cvdraws` function.
    #' @return cv_draws tibble augmented with the row proportion of each cell in each iteration.
    #'
  cv_draws |> dplyr::group_by(iter, area_no, row_no) |>
    dplyr::mutate(est_row_prop = est_cell_value/sum(est_cell_value))
}

summarize_cvdraws <- function(cv_draws){
  if("est_row_prop" %in% names(cv_draws)) {
    ret <- cv_draws |>
      dplyr::group_by(area_no, row_no, col_no) |>
      dplyr::summarise(
        est_cell_value_max90 = quantile(est_cell_value,probs =.95),
        est_cell_value_min90=quantile(est_cell_value, probs =.05),
        est_cell_value = median(est_cell_value),
        est_row_prop_max90 = quantile(est_row_prop, probs = .95, na.rm=TRUE),
        est_row_prop_min90 = quantile(est_row_prop, probs = .05, na.rm = TRUE),
        est_row_prop = median(est_row_prop, na.rm = TRUE))
  } else {
    ret <- cv_draws |>
      dplyr::group_by(area_no, row_no, col_no) |>
      dplyr::summarise(
        est_cell_value_max90 = quantile(est_cell_value,probs =.95),
        est_cell_value_min90=quantile(est_cell_value, probs =.05),
        est_cell_value = median(est_cell_value))
  }
}
