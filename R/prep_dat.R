#' Title
#'
#' @param row_margins
#' @param col_margins
#'
#' @return
#' @export
#'
#' @examples
prep_data_stan <- function(row_margins, col_margins, V_ilr, n_ilr_rows){
  R <- ncol(row_margins)
  C <- ncol(col_margins)
  n_areas <- nrow(row_margins)

  if(is.null(V_ilr)){
    V_ilr = matrix(0, C, C -1)
    n_ilr_rows = 0;
  } else {
    # check if V'V=I
    if(isFALSE(all.equal(crossprod(V_ilr), diag(ncol(V_ilr)),tolerance = 1e-5))) {
      stop("V_ilr must be an orthonormal matrix but crossprod(V_ilr) does not equal an identity matrix.")
    }
    V_ilr = V_ilr
  }

  standata <- list(
    n_areas = n_areas,
    R = R,
    C = C,
    row_margins = row_margins,
    col_margins = col_margins,
    V_ilr_data = V_ilr,
    n_ilr_rows = n_ilr_rows
  )
  return(standata)
}

#' Title
#'
#' @param struc_zero_rm
#' @param struc_zero_cm
#'
#' @return
#' @export
#'
#' @examples
prep_zeros <- function(struc_zero_rm="", struc_zero_cm=""){
  J = nrow(struc_zero_rm)
  R = ncol(struc_zero_rm)
  C = ncol(struc_zero_cm)

  struc_zero_array = array(dim=c(J, R, C))
  for(j in 1:J){
    for(r in 1:R){
      for(c in 1:C){
        struc_zero_array[j,r,c] = ifelse(struc_zero_rm[j,r]==0|struc_zero_cm[j,c]==0, 1, 0)
      }
    }
  }
  list(
    structural_zeros = struc_zero_array
  )
}

#' Title
#'
#' @param use_dist
#' @param area_re
#' @param inc_rm
#' @param predictors_rm
#' @param vary_sd
#' @param llmod_const
#' @param llmod_structure_omit
#'
#' @return
#' @export
#'
#' @examples
prep_options_stan <- function(use_dist,
                              area_re,
                              inc_rm,
                              vary_sd,
                              ll_rep,
                              mod_corr,
                              llmod_omit_jr,
                              llmod_omit_jc,
                              llmod_omit_jrc,
                              predictors_cm,
                              noncentred,
                              raw_seq_cell_weights
                              ){
  if(!use_dist %in% c("pois", "multinom", "negbinom", "multinomdirich")){
    stop("use_dist must be one of: pois, multinom, negbinom, multinomdirich")
  }
  if(use_dist == "multinomdirch"){
    stop("multinomdirich not yet implemented. Use negbinom instead.")
  }
  if(!area_re %in% c("none", "normal", "multinormal", "multinormal2")){
    stop("area_re must be one of: none, normal, multinormal, multinormal2")
  }
  if(!inc_rm %in% c(TRUE, FALSE)){
    stop("inc_rm must be one of: TRUE, FALSE")
  }
  if(!vary_sd %in% c(FALSE, TRUE, "partial")){
    stop("vary_sd must be one of: TRUE, FALSE, partial")
  }
  if(!mod_corr %in% c(FALSE, TRUE)){
    stop("mod_corr must be one of: TRUE, FALSE")
  }
  if(!ll_rep %in% c("ALR 1", "ILR 1", "ILR 2", "ILR 3", "LOR")){
    stop("ll_rep must be one of: ALR 1, ILR 1, ILR 2, ILR 3, LOR")
  }
  if(!llmod_omit_jr %in% c(TRUE, FALSE)){
    stop("llmod_omit_jr must be one of TRUE, FALSE")
  }
  if(!llmod_omit_jc %in% c(TRUE, FALSE)){
    stop("llmod_omit_jc must be one of TRUE, FALSE")
  }
  if(!llmod_omit_jrc %in% c(TRUE, FALSE)){
    stop("llmod_omit_jrc must be one of TRUE, FALSE")
  }
  if(!predictors_cm %in% c(TRUE, FALSE)){
    stop("predictors_cm must be one of: TRUE, FALSE")
  }
  if(!noncentred %in% c(TRUE, FALSE)){
    stop("noncentred must be one of: TRUE, FALSE")
  }
  if(!raw_seq_cell_weights %in% c(TRUE, FALSE)){
    stop("raw_seq_cell_weight must be one of: TRUE, FALSE")
  }

  list(
    lflag_dist = dplyr::case_when(
      use_dist == "pois" ~ 0,
      use_dist == "multinom" ~ 1,
      use_dist == "negbinom" ~ 2,
      use_dist == "multinomdirich" ~ 3
    ),
    lflag_area_re = dplyr::case_when(
      area_re == "none" ~ 0,
      area_re == "normal" ~ 1,
      area_re == "multinormal" ~ 2,
      area_re == "multinormal2" ~ 3
    ),
    lflag_inc_rm = dplyr::case_when(
      inc_rm == FALSE ~ 0,
      inc_rm == TRUE ~ 1
    ),
    lflag_predictors_cm = dplyr::case_when(
      predictors_cm == FALSE ~ 0,
      predictors_cm == TRUE ~ 1
    ),
    lflag_mod_corr = dplyr::case_when(
      mod_corr == FALSE ~ 0,
      mod_corr == TRUE ~ 1
    ),
    lflag_noncentred = dplyr::case_when(
      noncentred == FALSE ~ 0,
      noncentred == TRUE ~ 1
    ),
    lflag_rawscw = dplyr::case_when(
      raw_seq_cell_weights == FALSE ~ 0,
      raw_seq_cell_weights == TRUE ~ 1
    ),
    lflag_vary_sd = dplyr::case_when(
      vary_sd == FALSE ~ 0,
      vary_sd == TRUE ~ 1,
      vary_sd == "partial" ~ 2
    ),
    lflag_ll_rep = dplyr::case_when(
      ll_rep == "ALR 1" ~ 0,
      ll_rep == "ILR 1" ~ 1,
      ll_rep == "ILR 2" ~ 2,
      ll_rep == "ILR 3" ~ 4,
      ll_rep == "LOR" ~ 3
     ),
    lflag_llmod_omit_jr = dplyr::case_when(
      llmod_omit_jr == TRUE ~ 1,
      llmod_omit_jr == FALSE ~ 0
    ),
    lflag_llmod_omit_jc = dplyr::case_when(
      llmod_omit_jc == TRUE ~ 1,
      llmod_omit_jc == FALSE ~ 0
    ),
    lflag_llmod_omit_jrc = dplyr::case_when(
      llmod_omit_jrc == TRUE ~ 1,
      llmod_omit_jrc == FALSE ~ 0
    )
  )
}

#' Title
#'
#' @param prior_lkj
#' @param prior_mu_ce_scale
#' @param prior_mu_re_scale
#' @param prior_sigma_c_scale
#' @param prior_sigma_c_mu_scale
#' @param prior_sigma_ce_scale
#' @param prior_sigma_re_scale
#' @param prior_cell_effect_scale
#'
#' @return
#' @export
#'
#' @examples
prep_priors_stan <- function(prior_lkj,
                             prior_mu_ce_scale,
                             prior_mu_re_scale,
                             prior_sigma_c_scale,
                             prior_sigma_c_mu_scale,
                             prior_sigma_mu,
                             prior_sigma_ce_scale,
                             prior_sigma_re_scale,
                             prior_cell_effect_scale,
                             prior_lambda_raw_scale){
  list(prior_lkj = prior_lkj,
       prior_mu_ce_scale = prior_mu_ce_scale,
       prior_mu_re_scale = prior_mu_re_scale,
       prior_sigma_c_scale = prior_sigma_c_scale,
       prior_sigma_c_mu_scale = prior_sigma_c_mu_scale,
       prior_sigma_mu = prior_sigma_mu,
       prior_sigma_ce_scale = prior_sigma_ce_scale,
       prior_sigma_re_scale = prior_sigma_re_scale,
       prior_cell_effect_scale = prior_cell_effect_scale,
       prior_lambda_raw_scale = prior_lambda_raw_scale)
}

prep_king <- function(rm, cm){
  formula <- as.formula(paste0("cbind(", paste(names(cm), collapse = ","), ") ~ cbind(", paste(names(rm), collapse=","), ")"))
  list(
    formula = formula,
    data = cbind(cm, rm)
  )
}

prep_gq <- function(rm, cm){
  names(cm) <- paste0("col_no.", 1:ncol(cm))
  names(rm) <- paste0("row_no.", 1:ncol(rm))
  formula <- paste0(paste0(names(cm), collapse = ","), "~", paste0(names(rm), collapse = ","))
  list(
    formula = formula,
    data = cbind(cm, rm)
  )
}

#' @export
build_gm_ilr_basis <- function(C, ref_col = C) {
  # Reorder so reference category is last
  cat_order <- c(setdiff(1:C, ref_col), ref_col)

  # Column 1: geometric mean direction  unaffected by ordering
  v_gm <- rep(1/sqrt(C), C)

  if (C == 2) {
    return(matrix(v_gm, nrow=C, ncol=1))
  }

  # Build Helmert ILR basis C×(C-1)
  V_helm <- matrix(0, nrow=C, ncol=C-1)
  for (i in 1:(C-1)) {
    V_helm[1:i, i]  <- (1/i) * sqrt(i/(i+1))
    V_helm[i+1, i]  <- -sqrt(i/(i+1))
  }

  # Take first C-2 columns of ILR basis
  V_ratios <- matrix(V_helm[, 1:(C-2)], nrow=C)
  V <- cbind(v_gm, V_ratios)

  # Reorder rows so reference category is last
  V[cat_order, ]
}

#' @export
#'
mk_E_rc <- function(known_cell_values, V_ilr){
  J <- dim(known_cell_values)[1]
  R <- dim(known_cell_values)[2]
  C <- dim(known_cell_values)[3]
  D <- C - 1
  Y <- array(NA, dim = c(J, R, D))
  E_rc <- matrix(NA, nrow = R, ncol = D)
  sigma_jrc <- matrix(NA, nrow=R, ncol = D)
  for (j in 1:J){
    for(r in 1:R){
      comp <- known_cell_values[j, r, ]
      if(any(comp<= 0)){
        comp[comp<=0] <- 1e-3
      }
      comp <- comp/sum(comp)
      Y[j,r, ] <- log(comp) %*% V_ilr
    }
  }
  for(r in 1:R){
    for(d in 1:D){
      E_rc[r, d] <- mean(Y[,r, d])
      sigma_jrc[r, d] <- sd(Y[,r, d])
    }
  }
  return(
    list(
      E_rc = E_rc,
      simga_jrc = sigma_jrc
    )
  )
}
