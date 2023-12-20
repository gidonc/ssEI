prep_king <- function(cm, rm){
  formula <- as.formula(paste0("cbind(", paste(names(cm), collapse = ","), ") ~ cbind(", paste(names(rm), collapse=","), ")"))
  list(
    formula = formula,
    data = cbind(cm, rm)
  )
}

prep_gq <- function(cm, rm){
  names(cm) <- paste0("col_no.", 1:ncol(cm))
  names(rm) <- paste0("row_no.", 1:ncol(rm))
  formula <- paste0(paste0(names(cm), collapse = ","), "~", paste0(names(rm), collapse = ","))
  list(
    formula = formula,
    data = cbind(cm, rm)
  )
}
