#' Run ssEI ecological inference model.
#'
#' @param row_margins A data frame with the row margins data, with one row for each area, and with columns giving the row margins in the area.
#' @param col_margins A data frame with the column margins data, with one row for each area, and with columns giving the column margins in the area.
#' @param use_dist Which distribution to use to model the cell values relative to the mean estimate, options are: pois (Poisson), negbinom (Negative Binomial), multinom (Multinomial), multinomdirich (Multinomial Dirichlet - not yet implemented). Note Poisson/Multinomial and Negative Binomial/Multinomial Dirichlet are alternative paramaterizations of the same models.
#' @param area_re Are there a random effects in the model of the area means? Options are: none, normal, multinormal, multinormal2 (multinormal uses a non-centered paramaterization, multinormal2 uses an onion construction of the correlation matrix and a non-centred paramaterization)
#' @param inc_rm Include row margin log-ratios in the random effects model of the area means. Options are TRUE, FALSE
#' @param mod_cols Model columns to row as well as row to column rates;
#' @param predictors_rm Include row margin log-ratios in the model of area means. Options are TRUE, FALSE
#' @param prior_lkj Prior of the LKJ correlation matrix. Parameter must be real number greater than zero, where 1 is uniform across correlations.
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
                        use_dist = "negbinom",
                        area_re = "none",
                        inc_rm = FALSE,
                        mod_cols = FALSE,
                        predictors_rm = FALSE,
                        prior_lkj = 2,
                        cores = 4,
                        chains = 4,
                        verbose = TRUE, ...){

  standata <- prep_data_stan(row_margins, col_margins)
  standata <- modifyList(standata,
                         prep_options_stan(
                           use_dist = use_dist,
                           area_re = area_re,
                           inc_rm = inc_rm,
                           mod_cols = mod_cols,
                           predictors_rm = predictors_rm
                         ))
  standat <- modifyList(standata,
                        prep_priors_stan(
                          prior_lkj = prior_lkj
                        ))

  if(verbose){
    print(standata$R)
    print(standata$C)
    print(standata$n_areas)

  }

  mod <- stanmodels$ssEIdev
  if(verbose){
    print(paste("now running model", mod@model_name))
  }

  out <- rstan::sampling(mod,
                         data = standata,
                         cores = cores,
                         chains = chains,
                         ...)
}
