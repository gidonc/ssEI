#' Generate Row Margins for Crosstable Simulation
#'
#' @param n_areas The number of areas to generate row margins for.
#' @param n_row The number of rows to generate.
#' @param max_row_margin The maximum size of any row.
#'
#' @return
#' @export
#'
#' @examples
mk_row_margins <- function(n_areas, n_row, max_row_margin=1000){
  row_margins <- data.frame(matrix(nrow = n_areas, ncol = n_row))
  names(row_margins) <- paste0("row_no.", 1:n_row)
  for(j in 1:n_areas){
    row_margins[j,] <- sample(1:max_row_margin, n_row, replace = TRUE)
  }
  row_margins
}

#' Generate a Variance-Covariance Matrix for Crosstable Simulation
#'
#' Generate a variance-covariance matrix, for crosstable simulation. Parameters control the variance and covariance of the matrix. The variance-covariance matrix $\Sigma$ is constructed from a vector of length $K$ of standard deviations ($\sigma$) and a correlation matrix $\Omega$.
#' The vector of standard deviations $\sigma$ is constructed by from an `rnorm` call where $\log(\sigma) \sim N(mu_{\sigma}, \tau_{sigma}$.
#' The correlation matrix $\Omega$  is constructed using LKJ (2009) methods.
#' @param K The size of the square matrix.
#' @param scale_sd The mean from which the log standard deviations are drawn.
#' @param sd_sd The standard deviation from which
#' @param eta See `rlkjcorr` function.
#'
#' @return
#' @export
#'
#' @examples
mk_Sigma <- function(K, scale_sd = -2, sd_sd = 1, lkj_eta = 1){

  sd_diag <- exp(rnorm(K, scale_sd, sd_sd))
  Omega <- rlkjcorr(1, K, lkj_eta)
  Sigma <- diag(sd_diag) %*% Omega %*% diag(sd_diag)

  Sigma
}

mk_mu_raw <- function(n_row, n_col, mu_sigma = 1){
  mu_raw <- matrix(nrow = n_row, ncol = n_col -1)
  for (r in 1:n_row){
    mu_raw[r,] <- rnorm(n_col - 1, 0, mu_sigma)
  }
  mu_raw
}

mk_eta <- function(n_areas, mu_raw, n_row, n_col, Sigma){
  mysoftmax <- function(x){
    x <- c(0, x)
    exp(x)/sum(exp(x))
  }
  K <- n_row * (n_col - 1)
  eta_raw <- vector(length = n_areas, mode = "list")
  eta <- vector(length = n_areas, mode = "list")
  for (j in 1:n_areas){
    eta_raw[[j]] <- matrix(MASS::mvrnorm(1, rep(0, K), Sigma), nrow = n_row, ncol = n_col - 1)
    eta[[j]] <- matrix(nrow = n_row, ncol = n_col)
    for(r in 1:n_row){
      eta[[j]][r,] <- mysoftmax(mu_raw[r,] + eta_raw[[j]][r,])
    }
  }
  eta
}

sim_tables_from_probs <- function(n_areas, row_margins, n_col, eta){
  n_row <- ncol(row_margins)
  the_sim <- vector(length=n_areas, mode = "list")
  for(j in 1:n_areas){
    the_area_df <- data.frame(matrix(nrow = n_row, ncol = n_col))
    for(r in 1:n_row){
      the_area_df[r,] <- rmultinom(n =1, size = row_margins[j, r], prob = eta[[j]][r,])
    }
    names(the_area_df) <- paste0("col_no.", 1:ncol(the_area_df))
    the_area_df$area_no <- j
    the_area_df$row_no <- 1:n_row
    the_sim[[j]] <- the_area_df
  }
  the_sim_df <- dplyr::bind_rows(the_sim)

  col_margins <- the_sim_df |>
    tidyr::pivot_longer(cols = contains("col_no"),
                        names_to = "col_name") |>
    dplyr::group_by(area_no, col_name) |>
    dplyr::summarise(value = sum(value)) |>
    tidyr::pivot_wider(names_from = col_name, values_from = value) |>
    dplyr::ungroup() |>
    dplyr::select(-area_no)
  sim_df_long <- the_sim_df |>
    tidyr::pivot_longer(-c(area_no, row_no),
                        values_to = "actual_cell_value") |>
    tidyr::separate(name, into=c(NA, "col_no"), sep="[.]") |>
    dplyr::mutate(col_no = as.numeric(col_no))

  row_rates <- sim_df_long |>
    dplyr::left_join(row_margins |>
                       dplyr::mutate(area_no = dplyr::row_number()) |>
                       tidyr::pivot_longer(cols = -area_no, values_to = "row_margin") |>
                       tidyr::separate(name, into = c(NA, NA, "row_no"))  |>
                       dplyr::mutate(row_no = as.numeric(row_no))) |>
    dplyr::mutate(actual_row_rate = actual_cell_value/row_margin)

  list(
    sim_df_wide = the_sim_df,
    sim_df_long = sim_df_long,
    row_margins= row_margins,
    col_margins = col_margins,
    row_rates = row_rates,
    eta = eta
  )
}

#' Generate Crosstables from a number of different areas with related structures.
#'
#' @param n_areas Number of areas to generate.
#' @param n_row Number of rows in each area.
#' @param n_col Number of columns in each area.
#' @param max_row_margin Maximum number of cases in each row.
#' @param scale_sd The mean from which the log standard deviations are drawn.
#' @param sd_sd The standard deviation
#' @param contextual_effects Are there within and between row correlations in
#' @param lkj_eta Concentration parameter of matrix. Must be >0. `eta = 1` then the density is uniform over correlation matricies. `eta >1` the identify matrix is the modal correlation matrix (favours lower correlations). `0 < eta , 1` the density has a trough at the identity matrix (favours higher correlations).
#' @param mu_sigma
#'
#' @return
#' @export
#'
#' @examples
mk_sim_tables <- function(n_areas, n_row,  n_col,
                          max_row_margin=1000,
                          scale_sd = -2,
                          sd_sd = 1,
                          contextual_effects = TRUE,
                          lkj_eta = 1,
                          mu_sigma = 1){
  K <- n_row * (n_col - 1)

  row_margins <- mk_row_margins(n_areas, n_row, max_row_margin)

  if(contextual_effects){
    Sigma <- mk_Sigma(K, scale_sd = scale_sd, sd_sd = sd_sd, lkj_eta = lkj_eta)
  } else {
    sd_diag <- exp(rnorm(K, scale_sd, sd_sd))
    Sigma <- diag(sd_diag^2)
  }


  mu_raw <- mk_mu_raw(n_row, n_col, mu_sigma = mu_sigma)

  eta <- mk_eta(n_areas = n_areas, mu_raw = mu_raw, n_row = n_row, n_col = n_col, Sigma = Sigma)

  sim_tables <- sim_tables_from_probs(n_areas, row_margins, n_col, eta)

  sim_tables$Sigma <- Sigma

  sim_tables
}




#' Generate Correlation Matrix Using LKJ Method
#'
#' Ben Bolker's rLKJ function copied from Richard McElreath rethinking package.
#'
#' @param n Number of matricies to simulate
#' @param K Dimension of matrix
#' @param eta Concentration parameter of matrix. Must be >0. `eta = 1` then the density is uniform over correlation matricies. `eta >1` the identify matrix is the modal correlation matrix (favours lower correlations). `0 < eta , 1` the density has a trough at the identity matrix (favours higher correlations).
#'
#' @return
#' @export
#'
#' @examples
rlkjcorr <- function ( n , K , eta = 1 ) {


  stopifnot(is.numeric(K), K >= 2, K == as.integer(K))
  stopifnot(eta > 0)
  #if (K == 1) return(matrix(1, 1, 1))

  f <- function() {
    alpha <- eta + (K - 2)/2
    r12 <- 2 * rbeta(1, alpha, alpha) - 1
    R <- matrix(0, K, K) # upper triangular Cholesky factor until return()
    R[1,1] <- 1
    R[1,2] <- r12
    R[2,2] <- sqrt(1 - r12^2)
    if(K > 2) for (m in 2:(K - 1)) {
      alpha <- alpha - 0.5
      y <- rbeta(1, m / 2, alpha)

      # Draw uniformally on a hypersphere
      z <- rnorm(m, 0, 1)
      z <- z / sqrt(crossprod(z)[1])

      R[1:m,m+1] <- sqrt(y) * z
      R[m+1,m+1] <- sqrt(1 - y)
    }
    return(crossprod(R))
  }
  R <- replicate( n , f() )
  if ( dim(R)[3]==1 ) {
    R <- R[,,1]
  } else {
    # need to move 3rd dimension to front, so conforms to array structure that Stan uses
    R <- aperm(R,c(3,1,2))
  }
  return(R)
}
