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
        est_row_prop_max90 = quantile(est_row_prop, probs = .95),
        est_row_prop_min90 = quantile(est_row_prop, probs = .05),
        est_row_prop = median(est_row_prop))
  } else {
    ret <- cv_draws |>
      dplyr::group_by(area_no, row_no, col_no) |>
      dplyr::summarise(
        est_cell_value_max90 = quantile(est_cell_value,probs =.95),
        est_cell_value_min90=quantile(est_cell_value, probs =.05),
        est_cell_value = median(est_cell_value))
  }
}
