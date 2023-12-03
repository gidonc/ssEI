#' Ecological Inference Error Index
#'
#' Calculated the Ecological Inference Error Index which is the difference between actual and predicted values in each cell divided by the total of all cells (times 50 to avoid double counting errors).
#'
#' @param actual the actual value of the cell
#' @param predicted the predicted value of the cell
#'
#' @return real number error index
#' @export
#'
#' @examples
ei_error_index<- function(actual, predicted){
  N <- sum(actual)
  50 * sum(abs(predicted - actual))/N
}

mae <- function(actual, predicted){
  n <- length(actual)
  if(n == length(predicted)){
    ae <- abs(actual - predicted)
    res <- sum(ae)/n
  } else {
    print("actual and predicted vectors must be of same length")
    res <- NA
  }
  res
}

rmse <- function(actual, predicted){
  n <- length(actual)
  if(n == length(predicted)){
    se <- (actual - predicted)^2
    res <- sqrt(sum(se)/n)
  } else {
    print("actual and predicted vectors must be of same length")
    res <- NA
  }
  res
}
