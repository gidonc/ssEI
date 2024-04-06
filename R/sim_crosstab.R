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

mk_E_j_all <- function(n_areas, E_mu, Sigma){
  E_j_all <- matrix(nrow = n_areas, ncol=nrow(Sigma))
  for(j in 1:n_areas){
    E_j_all[j,]<- MASS::mvrnorm(1, E_mu, Sigma)
  }
  E_j_all
}

sim_tables_from_llmod <- function(log_e_cell_values){
  J <- dim(log_e_cell_values)[1]
  R <- dim(log_e_cell_values)[2]
  C <- dim(log_e_cell_values)[3]
  cell_values_flat <- rpois(J*R*C, exp(c(log_e_cell_values)))
  cell_values <- array(dim=c(J, R, C))
  counter = 0
  sim_df_long <- data.frame(area_no = numeric(), row_no = numeric(), col_no=numeric(), actual_cell_value=numeric())
  for(c in 1:C){
    for(r in 1:R){
      for(j in 1:J){
        counter = counter + 1
        cell_values[j, r, c] = cell_values_flat[counter]

        sim_df_long <- bind_rows(sim_df_long, data.frame(area_no = j, row_no = r, col_no = c, actual_cell_value =cell_values_flat[counter]))
      }
    }
  }
  the_sim_df <- sim_df_long |> dplyr::mutate(col_no = paste0("col_no_", col_no))
  sim_df_long <- sim_df_long |>
    dplyr::group_by(area_no) |>
    dplyr::mutate(actual_area_total = sum(actual_cell_value)) |>
    dplyr::group_by(area_no, col_no) |>
    dplyr::mutate(actual_col_margin = sum(actual_cell_value),
                  actual_col_rate = actual_cell_value/actual_col_margin) |>
    dplyr::group_by(area_no, row_no) |>
    dplyr::mutate(actual_row_margin = sum(actual_cell_value),
                  actual_row_rate = actual_cell_value/actual_row_margin) |>
    dplyr::ungroup()

  sim_df_wide <- sim_df_long |> dplyr::select(-actual_row_rate, -actual_col_rate, -actual_col_margin) |>
    dplyr::mutate(col_no = paste0("col_no_", col_no)) |>
    tidyr::pivot_wider(names_from = col_no, values_from = actual_cell_value)

  row_margins <- sim_df_long |>
    dplyr::group_by(area_no, row_name = paste0("row_no_", row_no)) |>
    dplyr::summarise(actual_row_margin = mean(actual_row_margin)) |>
    pivot_wider(names_from="row_name", values_from = actual_row_margin) |>
    dplyr::ungroup() |>
    dplyr::select(-area_no)
  col_margins <- sim_df_long |>
    dplyr::group_by(area_no, col_name = paste0("col_no_", col_no)) |>
    dplyr::summarise(actual_col_margin = mean(actual_col_margin)) |>
    pivot_wider(names_from="col_name", values_from = actual_col_margin) |>
    dplyr::ungroup() |>
    dplyr::select(-area_no)

  list(
    sim_df_wide = sim_df_wide,
    sim_df_long = sim_df_long,
    row_margins= row_margins,
    col_margins = col_margins,
    log_e_cell_values = log_e_cell_values
  )
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
mk_sim_tables_old <- function(n_areas, n_row,  n_col,
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


  mu_raw <- rnorm(K, mu_sigma, sd_sigma)

  E_j_all <- mk_E_j_all()

  eta <- mk_eta(n_areas = n_areas, mu_raw = mu_raw, n_row = n_row, n_col = n_col, Sigma = Sigma)

  sim_tables <- sim_tables_from_probs(n_areas, row_margins, n_col, eta)

  sim_tables$Sigma <- Sigma

  sim_tables
}

#' Generate Crosstables from a number of different areas with related structures.
#'
#' @param n_areas Number of areas to generate.
#' @param R Number of rows in each area.
#' @param C Number of columns in each area.
#' @param loc_jr location (mean) of row effects (in EI this is often large relative to column and cell reffects)
#' @param loc_jc location (mean) of col effects (in EI this is often large relative to cell effects, but small relative to row effects)
#' @param loc_jrc location (mean) of cell effects (in EI these are often small relative to row and col effects)
#' @param contextual_effects Are there within and between row correlations in
#' @param lkj_eta Concentration parameter of matrix. Must be >0. `eta = 1` then the density is uniform over correlation matricies. `eta >1` the identify matrix is the modal correlation matrix (favours lower correlations). `0 < eta , 1` the density has a trough at the identity matrix (favours higher correlations).
#' @param mu_sigma
#'
#' @return
#' @export
#'
#' @examples
mk_sim_tables <- function(n_areas, R,  C,
                          loc_jr = 0,
                          loc_jc = -1,
                          loc_jrc = -1,
                          loc_jr_sd = 0,
                          loc_jc_sd = -1,
                          loc_jrc_sd = -1,
                          jr_sd_sd = 1,
                          jc_sd_sd = 1,
                          jrc_sd_sd = 1,
                          av_area_pop = 1000,
                          contextual_effects = TRUE,
                          lkj_eta = 1,
                          mu_sigma = 1){
  K <- R * C -1

  E_jr_sd_raw <- rnorm(R -1, loc_jr_sd, jr_sd_sd)
  E_jc_sd_raw <- rnorm(C - 1, loc_jc_sd, jc_sd_sd)
  E_jrc_sd_raw <- rnorm((R - 1)* (C - 1), loc_jrc_sd, jrc_sd_sd)

  sd_diag <- exp(c(E_jr_sd_raw, E_jc_sd_raw, E_jrc_sd_raw))


  if(contextual_effects){
    Omega <- rlkjcorr(1, K, lkj_eta)
    Sigma <- diag(sd_diag) %*% Omega %*% diag(sd_diag)
  } else {
    Sigma <- diag(sd_diag^2)
  }
  E_mu_jr <- rnorm(R - 1, loc_jr, mu_sigma)
  E_mu_jc <- rnorm(C - 1, loc_jc, mu_sigma)
  E_mu_jrc <- rnorm((R - 1)*(C - 1), loc_jrc, mu_sigma)
  E_mu <- c(E_mu_jr, E_mu_jc, E_mu_jrc)
  E_j_all <- mk_E_j_all(n_areas, E_mu, Sigma)

  cell_values <- array(rep(0, n_areas*R*C), dim=c(n_areas, R, C))
  log_r_cell_value <- array(rep(0, n_areas*R*C), dim=c(n_areas, R, C))
  log_e_cell_value <- array(rep(0, n_areas*R*C), dim=c(n_areas, R, C))
  E_j_adj <- log(rpois(n_areas, av_area_pop) + .01)
  E_jr <- matrix(rep(rep(0, R), n_areas), nrow = n_areas, ncol=R)
  E_jc <- matrix(rep(rep(0, C), n_areas), nrow = n_areas, ncol=C)
  E_jrc <- array(rep(0, n_areas*R*C), dim=c(n_areas, R, C))
  E_j <- rnorm(n_areas, 3, 1)

  jrc_rstart <- vector(length = R)
  for(r in 1:R - 1){
    jrc_rstart[r] <- R + C- 2 + ((r - 1) *(C - 1))
  }
  for(j in 1:n_areas){
    E_jr[j, 1:R - 1] <- E_j_all[j, 1:R - 1]
    E_jc[j, 1:C - 1] <- E_j_all[j, R:(R + C - 2)]
    for(r in 1:(R - 1)){
      E_jrc[j, r, 1:C - 1] <- E_j_all[j, jrc_rstart[r]:(jrc_rstart[r] + C - 2)]
    }
    for(r in 1:R){
      for(c in 1:C){
        log_r_cell_value[j, r, c] <- E_jr[j, r] + E_jc[j, c] + E_jrc[j, r, c]
      }
    }
    E_j[j] <- E_j_adj[j] - matrixStats::logSumExp(c(log_r_cell_value[j, 1:R, 1:C]))

    for(r in 1:R){
      for(c in 1:C){
        log_e_cell_value[j, r, c] <- log_r_cell_value[j, r, c] + E_j[j]
        cell_values[j, r, c] <- rpois(1, exp(log_e_cell_value[j, r, c]))
      }
    }
  }


  sim_tables <- sim_tables_from_llmod(log_e_cell_value)
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
