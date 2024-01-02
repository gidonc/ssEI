#' Ecological Inference Error Index
#'
#' Calculated the Ecological Inference Error Index which is the difference between actual and predicted values in each cell divided by the total of all cells (times 50 to avoid double counting errors).
#'
#' @param actual the actual value of the cell
#' @param predicted the predicted value of the cell
#'
#' @return real number error index
#' @export
#'
#' @examples
ei_error_index<- function(actual, predicted){
  N <- sum(actual)
  50 * sum(abs(predicted - actual))/N
}

mae <- function(actual, predicted){
  n <- length(actual)
  if(n == length(predicted)){
    ae <- abs(actual - predicted)
    res <- sum(ae)/n
  } else {
    print("actual and predicted vectors must be of same length")
    res <- NA
  }
  res
}

rmse <- function(actual, predicted){
  n <- length(actual)
  if(n == length(predicted)){
    se <- (actual - predicted)^2
    res <- sqrt(sum(se)/n)
  } else {
    print("actual and predicted vectors must be of same length")
    res <- NA
  }
  res
}


mk_eval_plot <- function(cv_eval, rr_eval, cv_eval_stats, rr_eval_stats){

  cv_eval_stats_long <- cv_eval_stats |>
    pivot_longer(cols = - model) |>
    group_by(model) |>
    mutate(value = paste(name, round(value, 4)),
           y = row_number() *.05 *max(cv_eval$mean))

  p1 <- cv_eval |>
    ggplot(aes(actual_cell_value, mean)) +
    facet_wrap(~model)+
    geom_point()+
    geom_text(data = cv_eval_stats_long, x =.8*max(cv_eval$actual_cell_value), aes(label = value, y = y)) +
    theme_bw() +
    labs(x ="actual cell values",
         y = "model estimated cell values") +
    geom_abline(intercept = 0, slope = 1)

  rr_eval_stats_long <- rr_eval_stats |>
    pivot_longer(cols = - model) |>
    group_by(model) |>
    mutate(value = paste(name, round(value, 4)),
           y = row_number() *.05 *max(rr_eval$mean))

  p2 <- rr_eval|>
    filter(mean > -1) |>
    ggplot(aes(actual_row_rate, mean)) +
    facet_wrap(~model) +
    geom_point()+
    geom_text(data = rr_eval_stats_long, x =.9, aes(y=y, label = value)) +
    theme_bw() +
    labs(x ="actual row proportions",
         y = "model estimated row proportions") +
    geom_abline(intercept = 0, slope = 1)


  cv_eval_stats_f <- cv_eval |>
    group_by(model, row_no, col_no) |>
    summarise(cor = cor(actual_cell_value, mean)) |>
    mutate(value = paste("cor", round(cor, 4)),
           y = .05 *max(cv_eval$mean))

  p3 <- cv_eval |>
    ggplot(aes(actual_cell_value, mean)) +
    facet_grid(row_no ~ col_no + model) +
    geom_point()+
    theme_bw() +
    labs(x ="actual cell values",
         y = "model estimated cell values") +
    geom_text(data = cv_eval_stats_f, x =.8*max(cv_eval$actual_cell_value), aes(label = value, y = y)) +
    geom_abline(intercept = 0, slope = 1)

  rr_eval_stats_f <- rr_eval |>
    filter(!is.na(actual_row_rate)) |>
    group_by(model, row_no, col_no) |>
    summarise(cor = cor(actual_row_rate, mean)) |>
    mutate(value = paste("cor", round(cor, 4)),
           y = .05 *max(rr_eval$mean, na.rm=TRUE))

  p4 <- rr_eval|>
    filter(mean > -1) |>
    ggplot(aes(actual_row_rate, mean)) +
    facet_grid(row_no ~ col_no + model) +
    geom_point()+
    geom_text(data = rr_eval_stats_f, x =.8*max(rr_eval$actual_row_rate, na.rm=TRUE), aes(label = value, y = y)) +
    theme_bw() +
    labs(x ="actual row proportions",
         y = "model estimated row proportions") +
    geom_abline(intercept = 0, slope = 1)

  list(p_cv = p1,
       p_rr = p2,
       pf_cv = p3,
       pf_rr = p4)
}


#' Summarise EI Model Performance
#'
#' Summarise the performance of results from a set of Ecological Inference models against known actual cell and row rates (in long format).
#'
#' @param mod_list The models to evaluate
#' @param actual_long
#'
#' @return
#' @export
#'
#' @examples
mods_summary <- function(mod_list, actual_long){
  model_names <- names(mod_list)
  n_mod <- length(mod_list)

  cv_res <- list(length = n_mod)
  rr_res <- list(length = n_mod)

  for (n in 1:n_mod){
    this_class <- class(mod_list[[n]])
    if(this_class =="eiMD"){
      rr_res[[n]] <- md.bayes.join.rr(mod_list[[n]], actual_long)
      cv_res[[n]] <- md.bayes.cv(rr_res[[n]])

    } else if(this_class[[1]]=="stanfit"){
      cv_res[[n]] <- ei_cv_summary(mod_list[[n]])
      rr_res[[n]] <- ei_row_rate_summary(mod_list[[n]])
    } else if(this_class=="mcmc.list"){
      tmp_res <- link.gq.res(mod_list[[n]], actual_long)
      rr_res[[n]] <- tmp_res
      cv_res[[n]] <- tmp_res |> dplyr::mutate(mean = value)
    }
  }
  cv_eval <- mk_cv_eval(cv_list = cv_res,
                        actual_long = actual_long,
                        model_names = model_names)
  rr_eval <- mk_rr_eval(rr_list = rr_res,
                        actual_long = actual_long,
                        model_names = model_names)

  eval_plots <- mk_eval_plot(cv_eval$cv_eval, rr_eval$rr_eval, cv_eval$cv_eval_stats, rr_eval$rr_eval_stats)


  list(
    cv_res = cv_res,
    rr_res = rr_res,
    cv_eval = cv_eval,
    rr_eval = rr_eval,
    eval_plots = eval_plots,
    model_names = model_names
  )
}

mk_cv_eval <- function(cv_list, actual_long, model_names = c("average", "contextual")){

  for(m in 1:length(cv_list)){
    this_res <- cv_list[[m]] |>
      dplyr::left_join(actual_long) |>
      dplyr::mutate(model = model_names[[m]])
    if(m==1){
      cv_eval <- this_res
    } else{
      cv_eval <- dplyr::bind_rows(cv_eval, this_res)
    }
  }
  cv_eval_stats <- cv_eval |>
    dplyr::group_by(model) |>
    dplyr::summarise(cor = cor(actual_cell_value, mean),
                     ei_error_index= ei_error_index(actual_cell_value, mean),
                     mae = mae(mean, actual_cell_value),
                     rmse = rmse(actual_cell_value, mean))

  list(
    cv_eval = cv_eval,
    cv_eval_stats = cv_eval_stats
  )
}

mk_rr_eval <- function(rr_list, actual_long, model_names = c("average", "contextual")){

  for(m in 1:length(rr_list)){
    this_res <- rr_list[[m]] |>
      dplyr::left_join(actual_long) |>
      dplyr::mutate(model = model_names[[m]])
    if(m==1){
      rr_eval <- this_res
    } else{
      rr_eval <- dplyr::bind_rows(rr_eval, this_res)
    }
  }
  rr_eval_stats <- rr_eval |>
    dplyr::filter(!is.na(actual_row_rate)) |>
    dplyr::group_by(model) |>
    dplyr::summarise(cor = cor(actual_row_rate, mean, use="pairwise.complete.obs"),
                     ei_error_index = ei_error_index(actual_row_rate, mean),
                     mae = mae(actual_row_rate, mean),
                     rmse = rmse(actual_row_rate, mean))

  list(
    rr_eval = rr_eval,
    rr_eval_stats = rr_eval_stats
  )
}

md.bayes.join.rr <- function(md.bayes.res, actual_row_rate){
  k.out<-summary(md.bayes.res$draws$Beta)$statistics
  k.res <- k.out |>
    as.matrix() |>
    data.frame() |>
    tibble::rownames_to_column("rowname") |>
    tibble::tibble()
  k.res <- k.res |>
    tidyr::separate(rowname, into=c(NA, NA, NA, "row_no", NA, NA, "col_no", "area_no"), remove=FALSE) |>
    dplyr::mutate(
      col_no = as.numeric(col_no),
      row_no = as.numeric(row_no),
      area_no= as.numeric(area_no),
      mean = Mean
    ) |>
    dplyr::left_join(actual_row_rate, by=c("area_no", "row_no", "col_no"))

  k.res
}

md.bayes.cv <- function(md.bayes.rr){
  md.bayes.rr |>
    dplyr::mutate(mean = Mean*actual_row_margin)
}

link.gq.res <- function(gq.MCMC.list, actual_cv){
  NNs <- lapply(gq.MCMC.list, function(x){
    attr(x, "NN.internals")
  })
  for(x in 1:length(NNs)){
    NNs[[x]] <- NNs[[x]] |> as.data.frame ()|> mutate(chain = x)
  }
  NNs <- bind_rows(NNs)
  NNs.df <- as.data.frame(NNs) |>
    mutate(iter = row_number()) |>
    pivot_longer(cols = c(-iter, -chain)) |>
    separate(name, into = c(NA, NA, "area_no", NA, NA, "row_no", NA, NA, "col_no")) |>
    mutate(area_no = as.numeric(area_no),
           row_no = as.numeric(row_no),
           col_no = as.numeric(col_no)
    )
  NN.summary <- NNs.df |>
    group_by(area_no, row_no, col_no) |>
    summarise(value = mean(value)) |> left_join(actual_cv, by = c("area_no","col_no", "row_no")) |>
    dplyr::group_by(area_no, row_no) |>
    dplyr::mutate(actual_row_rate = actual_cell_value/sum(actual_cell_value),
                  est_row_rate = value/sum(value),
                  mean = est_row_rate) |>
    ungroup()

}

