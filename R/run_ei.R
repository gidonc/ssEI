#' Run ssEI ecological inference model.
#'
#' @param row_margins A data frame with the row margins data, with one row for each area, and with columns giving the row margins in the area.
#' @param col_margins A data frame with the column margins data, with one row for each area, and with columns giving the column margins in the area.
#' @param use_dist Which distribution to use to model the cell values relative to the mean estimate, options are: pois (Poisson), negbinom (Negative Binomial), multinom (Multinomial), multinomdirich (Multinomial Dirichlet - not yet implemented). Note Poisson/Multinomial and Negative Binomial/Multinomial Dirichlet are alternative paramaterizations of the same models.
#' @param area_re Are there a random effects in the model of the area means? Options are: none, normal, multinormal
#' @param ll_rep Which log-linear parameters are used to represent tables? Options are: "ALR 1" (C - 1 Additive log-ratios relative to the final cell in each of the R -1 free rows), "LOR" (Log-odds ratios)
#' @param inc_rm Include row margin log-ratios in the random effects model of the area means. Options are TRUE, FALSE
#' @param mod_cols Model columns by examining the whole table (v model conditional on rows looking at row to column rates)
#' @vary_sd Is the standard deviation of cell parameters error terms shared across the whole table (FALSE), does it vary by cell (TRUE) or is there a shared model ("partial")
#' @llmod_const What constraints (conditions) are placed on the log-linear model. Fewer constraints on the log-linear model, gives more opportunity for the model to learn from the data. Options are "none", "area_tot" (conditions the model to ensure that the area totals match the area totals in the data).
#' @llmod_structure_omit Which terms to omit from the table model log-linear structure. Options are none (saturated model), area*r*c (omit three-way interaction), area*c (omit area column interaction - so column effects are constant across areas)
#' @param predictors_rm Include row margin log-ratios in the model of area means. Options are TRUE, FALSE
#' @param prior_lkj Prior of the LKJ correlation matrix. Parameter must be real number greater than zero, where 1 is uniform across correlations.
#' @param prior_mu_ce_sigma Prior of the scale for the mean column effects parameter
#' @param prior_mu_re_sigma Prior of the scale for the mean row effects parameter
#' @param verbose Print some information about the data and model before running
#' @param method Either 'sampling' for NUTS or 'optimizing'
#' @param cores To be passed to rstan::sampling
#' @param chains To be passed to rstan::sampling
#' @param ... other arguments to be passed to rstan::sampling
#'
#' @return
#' @export
#'
#' @examples
ei_estimate <- function(row_margins, col_margins, E_rc_prior, known_cell_values,
                        E_rc_fixed = 0, sigma_jrc_fixed =0,
                        use_known_cells = 0,
                        fix_E_rc = 0, fix_sigma_jrc = 0,
                        zeros_structural = FALSE,
                        use_dist = "pois",
                        area_re = "normal",
                        ll_rep = "ALR 1",
                        V_ilr = NULL,
                        ROT = NULL,
                        ROT_E_rc = NULL,
                        V_dcorr_mat = NULL,
                        n_ilr_rows = NULL,
                        inc_rm = FALSE,
                        vary_sd = FALSE,
                        mod_cols = TRUE,
                        rotate_llrep = TRUE,
                        rotate_E_rc = FALSE,
                        neutral_logit = "table",
                        llmod_omit_jr = FALSE,
                        llmod_omit_jc = FALSE,
                        llmod_omit_jrc = FALSE,
                        pin_row_effect = FALSE,
                        predictors_cm = FALSE,
                        noncentred = "noncentred",
                        margin_reduction = NULL,
                        E_rc_hier = FALSE,
                        lambda_raw_offset = TRUE,
                        noncentred_mat = matrix(1, nrow=3, ncol=2),
                        family = "lognormal",
                        raw_seq_cell_weights = FALSE,
                        sigma_floor = 0,
                        hinge_delta_floor = 1e-8,
                        hinge_delta_min = 1e-8,
                        slack_tol = 1e-8,
                        prior_lkj = 2,
                        prior_mu_ce_scale = 2,
                        prior_mu_re_scale = 2,
                        prior_sigma_c_scale = 1,
                        prior_sigma_c_mu_scale = 1,
                        prior_sigma_mu = 0,
                        prior_sigma_ce_scale = 1,
                        prior_sigma_re_scale = 1,
                        prior_cell_effect_scale = 1,
                        prior_lambda_raw_scale = 3,
                        prior_gamma_rate = 2,
                        prior_gamma_shape = .5,
                        sample_optim = "sampling",
                        cores = 4,
                        chains = 4,
                        verbose = TRUE, ...){

  if(ll_rep == "ILR 3" & is.null(V_ilr)){
    stop("A V_ilr basis matrix is required for the ILR 3 representation of the log-linear model.")
  }
  if(ll_rep=="ILR 3" & is.null(n_ilr_rows)){
    if(is.null(n_ilr_rows)|nrow(V_ilr)<n_ilr_rows){
      stop("n_ilr_rows is required if for the ILR 3 representation of the log-linear model and must be less or equal to the number of rows in the V_ilr basis matrix.")
    }
  }


  standata <- prep_data_stan(row_margins, col_margins, V_ilr, n_ilr_rows)
  standata <- modifyList(standata,
                         prep_options_stan(
                           use_dist = use_dist,
                           area_re = area_re,
                           inc_rm = inc_rm,
                           vary_sd = vary_sd,
                           ll_rep = ll_rep,
                           llmod_omit_jr = llmod_omit_jr,
                           llmod_omit_jc = llmod_omit_jc,
                           llmod_omit_jrc = llmod_omit_jrc,
                           predictors_cm = predictors_cm,
                           noncentred = noncentred,
                           noncentred_mat = noncentred_mat,
                           raw_seq_cell_weights = raw_seq_cell_weights,
                           family = family,
                           E_rc_hier = E_rc_hier,
                           neutral_logit = neutral_logit,
                           lambda_raw_offset = lambda_raw_offset,
                           rotate_llrep = rotate_llrep,
                           rotate_E_rc = rotate_E_rc
                         ))
  standata <- modifyList(standata,
                         prep_priors_stan(
                           prior_lkj = prior_lkj,
                           prior_mu_ce_scale = prior_mu_ce_scale,
                           prior_mu_re_scale = prior_mu_re_scale,
                           prior_sigma_c_scale = prior_sigma_c_scale,
                           prior_sigma_c_mu_scale = prior_sigma_c_mu_scale,
                           prior_sigma_mu = prior_sigma_mu,
                           prior_sigma_ce_scale = prior_sigma_ce_scale,
                           prior_sigma_re_scale = prior_sigma_re_scale,
                           prior_cell_effect_scale = prior_cell_effect_scale,
                           prior_lambda_raw_scale = prior_lambda_raw_scale
                         ))
  if(zeros_structural == TRUE){
    zero_rm <- row_margins
    zero_cm <- col_margins
  } else{
    zero_rm <- row_margins + 1
    zero_cm <- col_margins + 1
  }
  standata <- modifyList(standata,
                         prep_zeros(zero_rm, zero_cm))
  standata <- modifyList(standata,
                         list(E_rc_prior = E_rc_prior,
                              known_cell_values = known_cell_values,
                              use_known_cells = use_known_cells,
                              E_rc_fixed = E_rc_fixed,
                              sigma_jrc_fixed = sigma_jrc_fixed,
                              lflag_fix_E_rc = fix_E_rc,
                              lflag_fix_sigma_jrc = fix_sigma_jrc,
                              sigma_floor = sigma_floor))
  if(standata$lflag_ll_rep %in% c(0, 3)){
    R_ll = standata$R - 1
    C_ll = standata$C - 1
  } else if(standata$lflag_ll_rep %in% c(1, 2)){
    R_ll = standata$R
    C_ll = standata$C - 1
  } else if(standata$lflag_ll_rep %in% c(4)){
    R_ll = n_ilr_rows
    C_ll = standata$C - 1
  }
  if(is.null(V_dcorr_mat)){
    V_dcorr_mat <- vector("list", standata$n_areas)
    for(j in 1:standata$n_areas){
      V_dcorr_mat[[j]] <- diag((standata$R - 1)*(standata$C - 1))
    }
  }
  if(!is.null(margin_reduction)){
    standata <- modifyList(standata, margin_reduction)
  }

  standata <- modifyList(standata,
                         list(
                           R_ll = R_ll,
                           C_ll = C_ll,
                           prior_gamma_shape = prior_gamma_shape,
                           prior_gamma_rate = prior_gamma_rate,
                           lflag_pin_row_effect = pin_row_effect,
                           ROT = ROT,
                           ROT_red = margin_reduction$ROT,
                           ROT_E_rc = ROT_E_rc,
                           V_dcorr_mat = V_dcorr_mat,
                           hinge_delta_floor = hinge_delta_floor,
                           hinge_delta_min = hinge_delta_min,
                           slack_tol = slack_tol
                         ))

  if(verbose){
    print(standata$R)
    print(standata$C)
    print(standata$n_areas)
    print(standata$V_ilr)
    print(standata$n_ilr_rows)
    print(standata$lflag_ll_rep)

  }
  if(mod_cols){
    mod <- stanmodels$ssEItable
  } else{
    mod <- stanmodels$ssEIrow
  }

  if(verbose){
    print(paste("now running model", mod@model_name))
  }

  if(sample_optim == "optim"){
    out <- rstan::optimizing(mod, data = standata, ...)
    class(out) <- c("ei_optim", "list")
  } else {
    out <- rstan::sampling(mod, data = standata, cores = cores, chains = chains, ...)
  }
  out
}

#' @export
build_margin_reduction <- function(V_ilr, row_margins, R, C, eps = 0.5,
                                   noncentred_mat = NULL) {
  margin_structure <- detect_margin_columns(V_ilr, R, C)
  free_dims        <- margin_structure$free
  Dm1_model        <- length(free_dims)
  n_areas          <- nrow(row_margins)
  D                <- R * C

  V_ilr_full  <- V_ilr
  V_ilr_model <- V_ilr[, free_dims]

  senc_smooth  <- as.matrix(row_margins) + eps
  agg_row_prop <- senc_smooth / rowSums(senc_smooth)

  # rebuild ROT for reduced dimension
  make_helmert_basis <- function(N) {
    V <- matrix(0, N, N - 1)
    for (k in 1:(N - 1)) {
      norm_const <- sqrt(k * (k + 1))
      V[1:k, k]  <-  1 / norm_const
      V[k+1,  k] <- -k / norm_const
    }
    V
  }

  V_area <- make_helmert_basis(n_areas)
  j_area <- matrix(1 / sqrt(n_areas), n_areas, 1)
  j_cell <- matrix(1 / sqrt(D), D, 1)

  B_agg <- kronecker(j_area, V_ilr_model)
  B_vol <- kronecker(V_area, j_cell)
  B_dev <- kronecker(V_area, V_ilr_model)
  V_nested <- cbind(B_agg, B_vol, B_dev)

  V_block_diag <- matrix(0, n_areas * D, n_areas * Dm1_model)
  for (j in 1:n_areas) {
    rows <- ((j-1)*D + 1):(j*D)
    cols <- ((j-1)*Dm1_model + 1):(j*Dm1_model)
    V_block_diag[rows, cols] <- V_ilr_model
  }
  V_flat   <- cbind(V_block_diag, B_vol)
  ROT_full <- t(V_flat) %*% V_nested

  Dtot_m1_reduced <- n_areas * Dm1_model + (n_areas - 1)
  stopifnot(all.equal(t(ROT_full) %*% ROT_full,
                      diag(Dtot_m1_reduced), tolerance = 1e-6))

  noncentred_mat_model <- if (!is.null(noncentred_mat)) {
    stopifnot(ncol(noncentred_mat) >= max(free_dims))
    noncentred_mat[, free_dims, drop = FALSE]
  } else {
    matrix(1L, nrow = n_areas, ncol = Dm1_model)   # default: all noncentred
  }

  list(
    n_agg_derived    = length(margin_structure$derived),
    agg_derived_dim  = as.array(margin_structure$derived),
    agg_free_dim     = as.array(free_dims),
    Dm1_model        = Dm1_model,
    n_agg_free       = Dm1_model,
    agg_row_prop     = agg_row_prop,
    V_ilr_full       = V_ilr_full,
    V_ilr_model      = V_ilr_model,
    ROT              = ROT_full,
    lflag_noncentred_mat = noncentred_mat_model   # n_areas x Dm1_model
  )
}
