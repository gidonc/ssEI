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
                        use_known_cells = 0,
                        zeros_structural = FALSE,
                        use_dist = "pois",
                        area_re = "normal",
                        ll_rep = "ALR 1",
                        V_ilr = NULL,
                        n_ilr_rows = NULL,
                        inc_rm = FALSE,
                        vary_sd = FALSE,
                        mod_cols = TRUE,
                        llmod_omit_jr = FALSE,
                        llmod_omit_jc = FALSE,
                        llmod_omit_jrc = FALSE,
                        predictors_cm = FALSE,
                        noncentred = TRUE,
                        raw_seq_cell_weights = FALSE,
                        prior_lkj = 2,
                        prior_mu_ce_scale = 2,
                        prior_mu_re_scale = 2,
                        prior_sigma_c_scale = 1,
                        prior_sigma_c_mu_scale = 1,
                        prior_sigma_ce_scale = 1,
                        prior_sigma_re_scale = 1,
                        prior_cell_effect_scale = 1,
                        prior_lambda_raw_scale = 3,
                        sample_optim = "sampling",
                        cores = 4,
                        chains = 4,
                        verbose = TRUE, ...){

  if(ll_rep == "ALR 3" & is.null(V_ilr)){
    stop("A V_ilr basis matrix is required for the ALR 3 representation of the log-linear model.")
  }
  if(ll_rep=="ALR 3" & is.null(n_ilr_rows)){
    if(is.null(n_ilr_rows)|nrow(V_ilr)<n_ilr_rows){
      stop("n_ilr_rows is required if for the ALR 3 representation of the log-linear model and must be less or equal to the number of rows in the V_ilr basis matrix.")
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
                           raw_seq_cell_weights = raw_seq_cell_weights
                         ))
  standata <- modifyList(standata,
                        prep_priors_stan(
                          prior_lkj = prior_lkj,
                          prior_mu_ce_scale = prior_mu_ce_scale,
                          prior_mu_re_scale = prior_mu_re_scale,
                          prior_sigma_c_scale = prior_sigma_c_scale,
                          prior_sigma_c_mu_scale = prior_sigma_c_mu_scale,
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
                              use_known_cells = use_known_cells))

  if(verbose){
    print(standata$R)
    print(standata$C)
    print(standata$n_areas)

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
