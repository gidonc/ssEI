#' Run ssEI ecological inference model.
#'
#' @param row_margins A data frame with the row margins data, with one row for each area, and with columns giving the row margins in the area.
#' @param col_margins A data frame with the column margins data, with one row for each area, and with columns giving the column margins in the area.
#' @param use_dist Which distribution to use to model the cell values relative to the mean estimate, options are: pois (Poisson), negbinom (Negative Binomial), multinom (Multinomial), multinomdirich (Multinomial Dirichlet - not yet implemented). Note Poisson/Multinomial and Negative Binomial/Multinomial Dirichlet are alternative paramaterizations of the same models.
#' @param area_re Are there a random effects in the model of the area means? Options are: none, normal, multinormal
#' @param inc_rm Include row margin log-ratios in the random effects model of the area means. Options are TRUE, FALSE
#' @param mod_cols Model columns by examining the whole table (v model conditional on rows looking at row to column rates)
#' @vary_sd Is the standard deviation of cell parameters error terms shared across the whole table (FALSE), does it vary by cell (TRUE) or is there a shared model ("partial")
#' @param predictors_rm Include row margin log-ratios in the model of area means. Options are TRUE, FALSE
#' @param prior_lkj Prior of the LKJ correlation matrix. Parameter must be real number greater than zero, where 1 is uniform across correlations.
#' @param prior_mu_ce_sigma Prior of the scale for the mean column effects parameter
#' @param prior_mu_re_sigma Prior of the scale for the mean row effects parameter
#' @param verbose Print some information about the data and model before running
#' @param cores To be passed to rstan::sampling
#' @param chains To be passed to rstan::sampling
#' @param ... other arguments to be passed to rstan::sampling
#'
#' @return
#' @export
#'
#' @examples
ei_estimate <- function(row_margins, col_margins,
                        use_dist = "pois",
                        area_re = "normal",
                        inc_rm = FALSE,
                        predictors_rm = FALSE,
                        vary_sd = FALSE,
                        mod_cols = FALSE,
                        prior_lkj = 2,
                        prior_mu_ce_scale = 2,
                        prior_mu_re_scale = 2,
                        prior_sigma_c_scale = 1,
                        prior_sigma_c_mu_scale = 1,
                        prior_sigma_ce_scale = 1,
                        prior_sigma_re_scale = 1,
                        prior_cell_effect_scale = 1,
                        cores = 4,
                        chains = 4,
                        verbose = TRUE, ...){

  standata <- prep_data_stan(row_margins, col_margins)
  standata <- modifyList(standata,
                         prep_options_stan(
                           use_dist = use_dist,
                           area_re = area_re,
                           inc_rm = inc_rm,
                           vary_sd = vary_sd,
                           predictors_rm = predictors_rm
                         ))
  standat <- modifyList(standata,
                        prep_priors_stan(
                          prior_lkj = prior_lkj,
                          prior_mu_ce_scale = prior_mu_ce_scale,
                          prior_mu_re_scale = prior_mu_re_scale,
                          prior_sigma_c_scale = prior_sigma_c_scale,
                          prior_sigma_c_mu_scale = prior_sigma_c_mu_scale,
                          prior_sigma_ce_scale = prior_sigma_ce_scale,
                          prior_sigma_re_scale = prior_sigma_re_scale,
                          prior_cell_effect_scale = prior_cell_effect_scale
                        ))

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

  out <- rstan::sampling(mod,
                         data = standata,
                         cores = cores,
                         chains = chains,
                         ...)
}
