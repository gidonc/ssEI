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

prep_options_stan <- function(use_dist,
                              area_re,
                              inc_rm,
                              predictors_rm){
  if(!use_dist %in% c("pois", "multinom", "negbinom", "multinomdirich")){
    stop("use_dist must be one of: pois, multinom, negbinom, multinomdirich")
  }
  if(use_dist == "multinomdirch"){
    stop("multinomdirich not yet implemented. Use negbinom instead.")
  }
  if(!area_re %in% c("none", "normal", "multinormal", "multinormal2")){
    stop("area_re must be one of: none, normal, multinormal", "multinormal2")
  }
  if(!inc_rm %in% c(TRUE, FALSE)){
    stop("inc_rm must be one of: TRUE, FALSE")
  }
  if(!predictors_rm %in% c(TRUE, FALSE)){
    stop("predictors_rm must be one of: TRUE, FALSE")
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
      area_re == "multinormal2" ~3
    ),
    lflag_inc_rm = dplyr::case_when(
      inc_rm == FALSE ~ 0,
      inc_rm == TRUE ~ 1
    ),
    lflag_predictors_rm = dplyr::case_when(
      predictors_rm == FALSE ~ 0,
      predictors_rm == TRUE ~ 1
    )
  )
}

prep_priors_stan <- function(prior_lkj){
  list(prior_lkj = prior_lkj)
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
