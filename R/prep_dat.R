#' Title
#'
#' @param row_margins
#' @param col_margins
#'
#' @return
#' @export
#'
#' @examples
prep_data_stan <- function(row_margins, col_margins){
  R <- ncol(row_margins)
  C <- ncol(col_margins)
  n_areas <- nrow(row_margins)

  standata <- list(
    n_areas = n_areas,
    R = R,
    C = C,
    row_margins = row_margins,
    col_margins = col_margins
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
                              predictors_rm,
                              vary_sd,
                              llmod_const,
                              llmod_structure_omit
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
  if(!predictors_rm %in% c(TRUE, FALSE)){
    stop("predictors_rm must be one of: TRUE, FALSE")
  }
  if(!vary_sd %in% c(FALSE, TRUE, "partial")){
    stop("vary_sd must be one of: TRUE, FALSE, partial")
  }
  if(!llmod_const %in% c("none", "area_tot")){
    stop("llmod_const must be one of none, area_tot")
  }
  if(!llmod_structure_omit %in% c("none", "area*r*c", "area*c")){
    stop("llmod_const must be one of none, area*r*c, area*c")
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
    lflag_predictors_rm = dplyr::case_when(
      predictors_rm == FALSE ~ 0,
      predictors_rm == TRUE ~ 1
    ),
    lflag_vary_sd = dplyr::case_when(
      vary_sd == FALSE ~ 0,
      vary_sd == TRUE ~ 1,
      vary_sd == "partial" ~ 2
    ),
    lflag_llmod_const = dplyr::case_when(
      llmod_const == "none" ~ 0,
      llmod_const == "area_tot" ~ 1
    ),
    lflag_llmod_structure = dplyr::case_when(
      llmod_structure_omit == "none" ~ 0,
      llmod_structure_omit == "area*r*c" ~ 1,
      llmod_structure_omit == "area*c" ~ 2
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
                             prior_sigma_ce_scale,
                             prior_sigma_re_scale,
                             prior_cell_effect_scale){
  list(prior_lkj = prior_lkj,
       prior_mu_ce_scale = prior_mu_ce_scale,
       prior_mu_re_scale = prior_mu_re_scale,
       prior_sigma_c_scale = prior_sigma_c_scale,
       prior_sigma_c_mu_scale = prior_sigma_c_mu_scale,
       prior_sigma_ce_scale = prior_sigma_ce_scale,
       prior_sigma_re_scale = prior_sigma_re_scale,
       prior_cell_effect_scale = prior_cell_effect_scale)
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
