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
#' @param return_data If TRUE, return the list of data for the Stan model instead of fitting it (for example to run the model in cmdstanr)
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
                        E_rc_node_logit = FALSE,
                        V_ilr = NULL,
                        ROT = NULL,
                        ROT_E_rc = NULL,
                        V_dcorr_mat = NULL,
                        n_ilr_rows = NULL,
                        inc_rm = FALSE,
                        mod_cols = TRUE,
                        rotate_llrep = TRUE,
                        rotate_E_rc = FALSE,
                        link_E_rc = NULL,
                        rm_ref = NULL,
                        rotate_lambda = "none",
                        ROT_lambda = NULL,
                        row_decompose = FALSE,
                        fit_type = "hybrid",
                        neutral_logit = "E_rc",
                        llmod_omit_jr = FALSE,
                        llmod_omit_jc = FALSE,
                        llmod_omit_jrc = FALSE,
                        pin_row_effect = FALSE,
                        predictors_cm = FALSE,
                        noncentred = "noncentred",
                        margin_reduction = NULL,
                        E_rc_hier = FALSE,
                        E_rc_ncp = FALSE,
                        E_rc_group_id = NULL,      # if NULL, built from E_rc_hier/E_rc_ncp defaults
                        E_rc_group_mode = NULL,
                        E_rc_group_fixed_value = 0,
                        E_rc_group_prior_a = NULL, E_rc_group_prior_b = NULL,
                        E_rc_group_tau_a = NULL,   E_rc_group_tau_b = NULL,
                        sigma_group_prior_a = NULL, sigma_group_prior_b = NULL,
                        sigma_group_tau_a = NULL,   sigma_group_tau_b = NULL,
                        vary_sd = FALSE,
                        sigma_ncp = FALSE,
                        sigma_group_id = NULL,
                        sigma_group_mode = NULL,
                        sigma_group_fixed_value = 0,
                        lambda_centred = TRUE,
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
                        verbose = TRUE, return_data = FALSE, ...){

  if(ll_rep == "ILR 3" & is.null(V_ilr)){
    stop("A V_ilr basis matrix is required for the ILR 3 representation of the log-linear model.")
  }
  if(ll_rep=="ILR 3" & is.null(n_ilr_rows)){
    if(is.null(n_ilr_rows)|nrow(V_ilr)<n_ilr_rows){
      stop("n_ilr_rows is required if for the ILR 3 representation of the log-linear model and must be less or equal to the number of rows in the V_ilr basis matrix.")
    }
  }
  if (fit_type == "cell-poisson" && neutral_logit != "llrep") {
    stop("fit_type = 'cell-poisson' requires neutral_logit = 'llrep'")
  }

  standata <- prep_data_stan(row_margins, col_margins, V_ilr, n_ilr_rows)
  standata <- modifyList(standata,
                         prep_options_stan(
                           use_dist = use_dist,
                           area_re = area_re,
                           inc_rm = inc_rm,
                           vary_sd = vary_sd,
                           ll_rep = ll_rep,
                           fit_type = fit_type,
                           llmod_omit_jr = llmod_omit_jr,
                           llmod_omit_jc = llmod_omit_jc,
                           llmod_omit_jrc = llmod_omit_jrc,
                           predictors_cm = predictors_cm,
                           noncentred = noncentred,
                           noncentred_mat = noncentred_mat,
                           raw_seq_cell_weights = raw_seq_cell_weights,
                           family = family,
                           E_rc_hier = E_rc_hier,
                           E_rc_node_logit = E_rc_node_logit,
                           neutral_logit = neutral_logit,
                           lambda_centred = lambda_centred,
                           rotate_llrep = rotate_llrep,
                           rotate_E_rc = rotate_E_rc,
                           rotate_lambda = rotate_lambda,
                           row_decompose = row_decompose
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

  if (rotate_lambda != "none") {
    if (is.null(ROT_lambda)) stop("ROT_lambda required when rotate_lambda is not 'none'")
    st <- attr(ROT_lambda, "stage")
    if (!is.null(st) && st != rotate_lambda)
      stop("ROT_lambda was built with stage '", st, "' but rotate_lambda is '",
           rotate_lambda, "'")
    nl <- attr(ROT_lambda, "neutral_logit_mode")
    if (!is.null(nl) && nl != standata$lflag_neutral_logit)
      stop("ROT_lambda was built for neutral_logit mode ", nl,
           " but this run uses ", standata$lflag_neutral_logit)
    n_rot_lam <- ncol(ROT_lambda)
  } else {
    n_rot_lam  <- 0
    ROT_lambda <- matrix(0, 0, 0)
  }
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
  if (!row_decompose && standata$Dm1_model != standata$R * standata$C - 1)
    stop("row_decompose = FALSE requires a full basis (Dm1_model == R*C - 1); ",
         "a margin-reduced basis has no coordinates for row proportions")
  if (row_decompose && standata$Dm1_model != standata$R * (standata$C - 1))
    stop("row_decompose = TRUE requires within-row dimensions only (Dm1_model == R*(C-1))")

  ## --- E_rc link ------------------------------------------------------------

  lflag_link_E_rc <- rep(0L, standata$Dm1_model)
  if (!is.null(link_E_rc)) {
    if (!all(link_E_rc %in% seq_len(standata$Dm1_model)))
      stop("link_E_rc must be indices between 1 and Dm1_model (", standata$Dm1_model, ")")
    lflag_link_E_rc[link_E_rc] <- 1L
  }
  standata$lflag_link_E_rc <- lflag_link_E_rc


  if (rotate_E_rc) {
    if (is.null(ROT_E_rc)) stop("ROT_E_rc required when rotate_E_rc = TRUE")
    if (!all(dim(ROT_E_rc) == standata$n_agg_free))
      stop("ROT_E_rc must be n_agg_free x n_agg_free (", standata$n_agg_free, "), got ",
           paste(dim(ROT_E_rc), collapse = "x"))
    if (!isTRUE(all.equal(t(ROT_E_rc) %*% ROT_E_rc, diag(ncol(ROT_E_rc)), tolerance = 1e-6)))
      stop("ROT_E_rc is not orthonormal")
  } else {
    ROT_E_rc <- diag(standata$n_agg_free)   # identity default, so it's always well-sized even when unused
  }

  standata$lflag_rot_E_rc <- as.integer(rotate_E_rc)
  standata$ROT_E_rc <- ROT_E_rc

  if (is.null(E_rc_group_id)) {
    d <- default_E_rc_groups(standata$n_agg_free, E_rc_hier, E_rc_ncp)
    E_rc_group_id <- d$group_id; E_rc_group_mode <- d$group_mode; E_rc_group_fixed_value <- d$group_fixed_value
  }
  ## after margin_reduction is merged and groups have defaults:
  erg <- resolve_groups(E_rc_group_id, E_rc_group_mode, E_rc_group_fixed_value,
                        E_rc_group_prior_a %||% 0,
                        E_rc_group_prior_b %||% prior_mu_re_scale,
                        E_rc_group_tau_a   %||% prior_gamma_shape,
                        E_rc_group_tau_b   %||% prior_gamma_rate,
                        standata$n_agg_free, "E_rc")

  lognormal <- family != "Gamma"
  sg <- resolve_groups(sigma_group_id, sigma_group_mode, sigma_group_fixed_value,
                       sigma_group_prior_a %||% (if (lognormal) prior_sigma_mu      else prior_gamma_shape),
                       sigma_group_prior_b %||% (if (lognormal) prior_sigma_c_scale else prior_gamma_rate),
                       sigma_group_tau_a   %||% prior_gamma_shape,
                       sigma_group_tau_b   %||% prior_gamma_rate,
                       standata$Dm1_model, "sigma_jrc")

  standata[c("E_rc_group_id", "E_rc_n_groups", "E_rc_group_mode", "E_rc_group_fixed_value",
             "E_rc_group_prior_a", "E_rc_group_prior_b", "E_rc_group_tau_a", "E_rc_group_tau_b")] <-
    erg[c("group_id", "n_groups", "group_mode", "group_fixed_value",
          "group_prior_a", "group_prior_b", "group_tau_a", "group_tau_b")]
  standata[c("sigma_group_id", "sigma_n_groups", "sigma_group_mode", "sigma_group_fixed_value",
             "sigma_group_prior_a", "sigma_group_prior_b", "sigma_group_tau_a", "sigma_group_tau_b")] <-
    sg[c("group_id", "n_groups", "group_mode", "group_fixed_value",
         "group_prior_a", "group_prior_b", "group_tau_a", "group_tau_b")]


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
                           ROT_lambda = ROT_lambda,
                           n_rot_lam = n_rot_lam,
                           V_dcorr_mat = V_dcorr_mat,
                           hinge_delta_floor = hinge_delta_floor,
                           hinge_delta_min = hinge_delta_min,
                           slack_tol = slack_tol
                         ))

  if(verbose){
    print(standata$R)
    print(standata$C)
    print(standata$n_areas)
    print(standata$V_ilr_model)

  }
  if(return_data){
    return(standata)
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

  out <- list(
    n_agg_derived    = length(margin_structure$derived),
    agg_derived_dim  = as.array(margin_structure$derived),
    agg_free_dim     = as.array(free_dims),
    Dm1_model        = Dm1_model,
    n_agg_free       = Dm1_model,
    agg_row_prop     = agg_row_prop,
    V_ilr_full       = V_ilr_full,
    V_ilr_model      = V_ilr_model,
    ROT              = ROT_full,
    lflag_rot_agg    = 0L,
    ROT_agg          = diag(Dm1_model),
    lflag_noncentred_mat = noncentred_mat_model   # n_areas x Dm1_model
  )
  modifyList(out, margin_param_defaults(Dm1_model, R, C))
}

#' @export
build_rot_E_rc <- function(pilot_fit, n_agg_free) {
  Er <- as.matrix(pilot_fit, pars = "E_rc_raw")
  stopifnot(ncol(Er) == n_agg_free)
  cv <- cov(Er)
  ev <- eigen(cv)
  stopifnot(all(ev$values > 1e-10))          # genuinely positive-definite, not near-singular
  ROT <- ev$vectors
  stopifnot(isTRUE(all.equal(t(ROT) %*% ROT, diag(n_agg_free), tolerance = 1e-6)))
  attr(ROT, "n_agg_free") <- n_agg_free
  ROT
}

#'@export
## 1. Defaults so that every existing call still supplies the new data fields.
##    Add these to the list returned by build_margin_reduction() (same place
##    ROT_agg / lflag_rot_agg are added), or merge them onto mr before ei_estimate().
margin_param_defaults <- function(Dm1_model, R = NULL, C = NULL) {
  out <- list(
    lflag_margin_param = 0L,
    n_pin    = 0L,
    mp_b_ref = array(0, c(0, Dm1_model)),
    mp_G     = array(0, c(0, 0, Dm1_model)),
    mp_P     = array(0, c(0, Dm1_model, 0)),
    mp_N     = array(0, c(0, Dm1_model, Dm1_model))
  )
  if (!is.null(R) && !is.null(C)) out <- c(out, rowwise_defaults(R, C, Dm1_model))
  out
}

#' Defaults for the data fields of the sequential (row-by-row) margin model
#'
#' Every switch is off and every optional input is empty, so a model built with
#' these behaves as it did before the fields existed. `build_margin_reduction()`
#' merges them in; change individual entries afterwards to switch a feature on.
#'
#' @param R,C number of rows and columns of the table
#' @param Dm1_model number of ILR coordinates in the hierarchy
#' @export
rowwise_defaults <- function(R, C, Dm1_model) {
  K <- (R - 1) * (C - 1)
  list(
    ## sequential expected table
    lflag_mp_seq = 0L, lflag_mp_seq_anchor = 0L, lflag_mp_exact = 0L,
    mp_newton_iters = 5L, lflag_mp_gamma_centred = 0L,
    lflag_seq_or_expected = 0L,                       # 3 = row by row, 4 = adjusted table (raking style)
    ## what areas may do beyond moving whole rows and columns: 0 = free interior, 1 = interior fixed at E_rc's (raking
    ## E_rc's table to the area's margins; needs lflag_seq_or_expected = 4), 2 = the rest scaled by kappa
    lflag_mp_interior = 0L, mp_kappa_fixed = 0, prior_mp_kappa_a = log(0.5), prior_mp_kappa_b = 0.5,
    mp_row_order = as.array(seq_len(R)),              # allocation order of the rows (last = remainder row)
    mp_rem_col   = as.array(pmin(seq_len(R), C)),     # reference column of each row
    kink_delta_expected = 0, kink_delta_realised = 0,
    ## scaling of the margin and interior parameters
    lflag_mp_vol_scale = 0L,
    lflag_mp_beta_whiten = 0L, mp_T  = array(0, c(0, C - 1, C - 1)),
    lflag_mp_row_whiten  = 0L, mp_Tr = array(0, c(0, R - 1, R - 1)),
    lflag_mp_seq_scale = 0L, lflag_mp_scale_fast = 0L,
    mp_A = array(0, c(0, K, Dm1_model)),
    mp_B = array(0, c(0, Dm1_model, K)),
    mp_nc_w = numeric(0), mp_sigma0 = numeric(0),
    ## HYBRID: realised table row by row
    lflag_real_rowwise = 0L, lflag_real_scale = 0L, lflag_cpois_norm = 0L,
    rl_S = array(0, c(0, K, K)),
    ## small-cell term (soft multinomial) and its lookup table
    lflag_soft_smallcell = 0L, smallcell_scale = 0,
    cz_n = 0L, cz_t0 = 0, cz_h = 1, cz_v = numeric(0), cz_d = numeric(0),
    ## experimental: covariate in the hierarchy mean, area-level column effects
    lflag_E_cov = 0L, E_cov_x = array(0, c(0, Dm1_model)), prior_E_cov_scale = 0.5,
    lflag_col_eff = 0L, lflag_col_tau_shared = 0L, prior_col_tau_scale = 0.5,
    n_col_groups = 0L, col_group = integer(0),
    col_eff_known = array(0, c(0, C)),
    lflag_col_cov = 0L, lflag_col_kappa_shared = 1L, col_x = array(0, c(0, C)),
    prior_col_kappa_scale = 2
  )
}

#' @export
## 2. Switch margin mode on, given per-area maps
add_margin_param <- function(mr, maps) {
  stopifnot(length(maps) == nrow(mr$lflag_noncentred_mat))   # one map per area
  arr <- function(f) aperm(simplify2array(lapply(maps, `[[`, f)), c(3, 1, 2))
  n_pin <- nrow(maps[[1]]$G)
  mr$lflag_margin_param <- 1L
  mr$n_pin    <- n_pin
  mr$mp_b_ref <- t(sapply(maps, `[[`, "b_ref"))          # n_areas x Dm1
  mr$mp_G     <- arr("G")                                 # n_areas x n_pin x Dm1
  mr$mp_P     <- arr("P")                                 # n_areas x Dm1 x n_pin
  mr$mp_N     <- arr("N")                                 # n_areas x Dm1 x (Dm1 - n_pin)
  ## guard against n_pin == 1 collapsing a dim in simplify2array
  stopifnot(identical(dim(mr$mp_G), c(length(maps), n_pin, mr$Dm1_model)),
            identical(dim(mr$mp_P), c(length(maps), mr$Dm1_model, n_pin)),
            identical(dim(mr$mp_N), c(length(maps), mr$Dm1_model, mr$Dm1_model - n_pin)))
  mr
}

#'@export
## 3. Inits: start each area at its raked reference (beta = 0 => b = b_ref, z = 0)
make_mp_init <- function(row_margins, n_pin, Dm1_model) {
  tot <- log(rowSums(as.matrix(row_margins)))
  n_areas <- length(tot)
  function() list(
    log_volume_raw  = tot[1],
    log_volume_rest = as.array(tot[-1]),
    mp_beta = matrix(0, n_areas, n_pin),
    mp_z    = matrix(0, n_areas, Dm1_model - n_pin)
  )
}
