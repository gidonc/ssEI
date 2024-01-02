library(ei)
library(RxCEcolInf)
library(tidyverse)


data(senc)

# senc2 <- dplyr::filter (senc, natam>0, white>0, black>0)
senc2 <- senc


vote_order <- tibble::tibble(vote = c("dem", "rep", "non"), col = c(1, 2, 3))
race_order <- tibble::tibble(race = c("wh", "bl", "natam"), row = c(1, 2, 3))

senc2_flat <- senc2 |> dplyr::select(bldem, blrep, blnon, whdem, whrep, whnon, natamdem, natamrep, natamnon) |>
  dplyr::mutate(area_no = dplyr::row_number()) |>
  tidyr::pivot_longer(-area_no, values_to = "actual_value") |>
  tidyr::separate(name, into=c("race", "vote"), sep=-3) |>
  dplyr::left_join(vote_order) |>
  dplyr::left_join(race_order)  |>
  dplyr::arrange(area_no, row, col) |>
  dplyr::mutate(CVN = dplyr::row_number(),
                row_names = race,
                col_names = vote) |>
  dplyr::group_by(area_no, row) |>
  dplyr::mutate(actual_row_total = sum(actual_value)) |>
  dplyr::group_by(area_no, col) |>
  dplyr::mutate(actual_col_total = sum(actual_value)) |>
  dplyr::ungroup() |>
  dplyr:: mutate(actual_row_prop = actual_value/actual_row_total) |>
  dplyr:: mutate(actual_col_prop = actual_value/actual_col_total)


senc_rm <- senc2 |> dplyr::select(white, black, natam)
senc_cm <- senc2 |> dplyr::select(dem, rep, non)
senc_n_areas <- nrow(senc2)
senc_R <- ncol(senc_rm)
senc_C <- ncol(senc_cm)

senc_numbers <- tibble::tribble(
  ~sencname, ~row_no, ~col_no, ~row_name, ~col_name,
  "whdem", 1, 1, "wh", "dem",
  "whrep", 1, 2, "wh", "rep",
  "whnon", 1, 3, "wh", "non",
  "bldem", 2, 1, "bl", "dem",
  "blrep", 2, 2, "bl", "rep",
  "blnon", 2, 3, "bl", "non",
  "natamdem", 3, 1, "natam", "dem",
  "natamrep", 3, 2, "natam", "rep",
  "natamnon", 3, 3, "natam", "non"
)



senc_long <- senc2 |>
  dplyr::select(whdem, whrep, whnon, bldem, blrep, blnon, natamdem, natamrep, natamnon) |>
  dplyr::mutate(area_no = dplyr::row_number()) |>
  tidyr::pivot_longer(cols = - area_no, names_to = "sencname", values_to = "actual_cell_value") |>
  dplyr::left_join(senc_numbers) |>
  dplyr::group_by(area_no, row_no) |>
  dplyr::mutate(actual_row_prop = actual_cell_value/sum(actual_cell_value))


ssei.path <- "C:/Users/gidon/OneDrive/Documents/ssEI"
smod.path <- paste0(ssei.path, "/inst/stan/")
mnb.path <- paste0(smod.path, "ssMNB.stan")
mnb2.path <- paste0(smod.path, "ssMNB2.stan")
co.new.path <- paste0(smod.path, "ssContextual2.stan")
flexmod.path <- paste0(smod.path, "ssEIdev.stan")

senc_for_stan <- list(
  n_areas = nrow(senc_rm),
  R = ncol(senc_rm),
  C = ncol(senc_cm),
  row_margins = senc_rm,
  col_margins = senc_cm,
  lkj_param = 2,
  lflag_mn = 1
)
senc_for_stanp <- list(
  n_areas = nrow(senc_rm),
  R = ncol(senc_rm),
  C = ncol(senc_cm),
  row_margins = senc_rm,
  col_margins = senc_cm,
  lkj_param = 2,
  lflag_mn = 0
)
senc_for_stannb <- list(
  n_areas = nrow(senc_rm),
  R = ncol(senc_rm),
  C = ncol(senc_cm),
  row_margins = senc_rm,
  col_margins = senc_cm,
  lkj_param = 2,
  lflag_mn = 2
)

senc_for_stannb_noarea_re <- list(
  n_areas = nrow(senc_rm),
  R = ncol(senc_rm),
  C = ncol(senc_cm),
  row_margins = senc_rm,
  col_margins = senc_cm,
  lkj_param = 2,
  lflag_mn = 2,
  lflag_area_re =0,
  lflag_inc_rm = 0,
  lflag_rm_predictors = 0
)
sencd <- list()
sencd[[1]]<-senc_for_stan_mn <- list(
  n_areas = nrow(senc_rm),
  R = ncol(senc_rm),
  C = ncol(senc_cm),
  row_margins = senc_rm,
  col_margins = senc_cm,
  lkj_param = 2,
  lflag_mn = 1,
  lflag_area_re =1,
  lflag_inc_rm = 0,
  lflag_predictors_rm = 0
)
sencd[[2]] <- senc_for_stan_p_inc_rm <- list(
  n_areas = nrow(senc_rm),
  R = ncol(senc_rm),
  C = ncol(senc_cm),
  row_margins = senc_rm,
  col_margins = senc_cm,
  lkj_param = 2,
  lflag_mn = 0,
  lflag_area_re =1,
  lflag_inc_rm = 1,
  lflag_predictors_rm = 0
)

sencd[[3]] <-senc_for_stannb_predictors_rm_noarea_re <- list(
  n_areas = nrow(senc_rm),
  R = ncol(senc_rm),
  C = ncol(senc_cm),
  row_margins = senc_rm,
  col_margins = senc_cm,
  lkj_param = 2,
  lflag_mn = 2,
  lflag_area_re =0,
  lflag_inc_rm = 0,
  lflag_predictors_rm = 1
)

sencd[[4]] <- senc_for_stannb_contextual_nb <- list(
  n_areas = nrow(senc_rm),
  R = ncol(senc_rm),
  C = ncol(senc_cm),
  row_margins = senc_rm,
  col_margins = senc_cm,
  lkj_param = 2,
  lflag_mn = 2,
  lflag_area_re =2,
  lflag_inc_rm = 0,
  lflag_predictors_rm = 0
)
sencd[[5]] <- senc_for_stannb_contextualonion_nb <- list(
  n_areas = nrow(senc_rm),
  R = ncol(senc_rm),
  C = ncol(senc_cm),
  row_margins = senc_rm,
  col_margins = senc_cm,
  lkj_param = 2,
  lflag_mn = 2,
  lflag_area_re =3,
  lflag_inc_rm = 0,
  lflag_predictors_rm = 0
)
sencd[[6]] <- senc_for_stannb_contextual_inc_rm <- list(
  n_areas = nrow(senc_rm),
  R = ncol(senc_rm),
  C = ncol(senc_cm),
  row_margins = senc_rm,
  col_margins = senc_cm,
  lkj_param = 2,
  lflag_mn = 2,
  lflag_area_re =2,
  lflag_inc_rm = 1,
  lflag_predictors_rm = 0
)
sencd[[7]] <- senc_for_stannb_contextual_rm_predictors <- list(
  n_areas = nrow(senc_rm),
  R = ncol(senc_rm),
  C = ncol(senc_cm),
  row_margins = senc_rm,
  col_margins = senc_cm,
  lkj_param = 2,
  lflag_mn = 2,
  lflag_area_re =2,
  lflag_inc_rm = 0,
  lflag_predictors_rm = 1
)
sencd[[8]] <-senc_for_stannb_predictors_rm_noarea_re <- list(
  n_areas = nrow(senc_rm),
  R = ncol(senc_rm),
  C = ncol(senc_cm),
  row_margins = senc_rm,
  col_margins = senc_cm,
  lkj_param = 2,
  lflag_mn = 2,
  lflag_area_re =0,
  lflag_inc_rm = 0,
  lflag_predictors_rm = 0
)

sencr<- list()

sencr[[6]] <- rstan::stan(file = flexmod.path,
                          data = sencd[[6]],
                          cores = 4
)
sencr[[7]] <- rstan::stan(file = flexmod.path,
                          data = sencd[[7]],
                          cores = 4
)

sencr[[1]] <- rstan::stan(file = flexmod.path,
                          data = sencd[[1]],
                          cores = 4,
                          )
sencr[[4]] <- rstan::stan(file = flexmod.path,
                          data = sencd[[4]],
                          cores = 4
)
sencr[[3]] <- rstan::stan(file = flexmod.path,
                          data = sencd[[3]],
                          cores = 4
)


sencr[[2]] <- rstan::stan(file = flexmod.path,
                          data = sencd[[2]],
                          cores = 4
)
sencr[[5]] <- rstan::stan(file = flexmod.path,
                          data = sencd[[5]],
                          cores = 4
)


st1 <- mods_summary(list("cont_nb" = sencr[[1]],
                         "contnb_rm" =sencr[[6]],
                         "contnb_rm2" =sencr[[7]]
                         ), senc2_flat |>
                      rename(col_no=col, row_no=row,
                             actual_cell_value = actual_value,
                             actual_row_rate = actual_row_prop))

st1$cv_eval
st1$rr_eval

st1$eval_plots$pf_rr
st1$eval_plots$p_rr + geom_point(aes(colour=paste(row_names, col_names), size=actual_cell_value))
st1$eval_plots$p_cv

sencr[[3]] <- rstan::stan(file = flexmod.path,
                          data = sencd[[3]],
                          chains = 1,
                          iter=10
)



senc_mnb <- rstan::stan(file = mnb.path,
                                    data = senc_for_stan,
                                    cores = 4,
                                    control = list(
                                      "adapt_delta" = .9)
)

senc_mnb2_mn <- rstan::stan(file = mnb2.path,
                        data = senc_for_stan,
                        cores = 4
)
senc_mnb2_p <- rstan::stan(file = mnb2.path,
                            data = senc_for_stanp,
                            cores = 4
)
senc_mnb2_nb <- rstan::stan(file = mnb2.path,
                            data = senc_for_stannb,
                            cores = 4
)
senc_mnb2_nb2 <- rstan::stan(file = mnb2.path,
                            data = senc_for_stannb_noarea_re,
                            cores = 4
)
senc_flex_contextual_nb <- rstan::stan(file = flexmod.path,
                             data = senc_for_stannb_contextual_nb,
                             cores = 4
)
senc_flex_contextualonion_nb <- rstan::stan(file = flexmod.path,
                                       data = senc_for_stannb_contextualonion_nb,
                                       cores = 4
)
senc_flex_contextual_nb_inc_rm <- rstan::stan(file = flexmod.path,
                                            data = senc_for_stannb_contextual_inc_rm,
                                            cores = 4
)
senc_co_new <- test_mnb <- rstan::stan(file = co.new.path,
                                    data = senc_for_stan,
                                    cores = 4,
                                    iter =50)

st1 <- mods_summary(list("mnb2_mn" = senc_mnb2_mn,
                         "mnb2_p" = senc_mnb2_p,
                         "mnb2_nb" = senc_mnb2_nb,
                         "mnb2_nb_nore" = senc_mnb2_nb2,
                         "cont_nb" = senc_flex_contextual_nb), senc2_flat |>
                      rename(col_no=col, row_no=row,
                             actual_cell_value = actual_value,
                             actual_row_rate = actual_row_prop))

st1 <- mods_summary(list("cont_nb" = senc_flex_contextual_nb,
                         "contonion_nb" =senc_flex_contextualonion_nb), senc2_flat |>
                      rename(col_no=col, row_no=row,
                             actual_cell_value = actual_value,
                             actual_row_rate = actual_row_prop))

st1$cv_eval
st1$rr_eval

st1$eval_plots$pf_rr
st1$eval_plots$p_rr + geom_point(aes(colour=paste(row_names, col_names), size=actual_cell_value))
st1$eval_plots$p_cv

senc_average <- ei_estimate(senc_rm, senc_cm, model="average")
senc_contextual <- ei_estimate(senc_rm, senc_cm)
senc_contextualonion <- ei_estimate(senc_rm, senc_cm, model="contextualOnion")
e_senc_average <- rstan::extract(senc_average)
e_senc_bb <- rstan::extract(bb)
e_senc_context <- rstan::extract(senc_contextual)
cv_draws <- ei_to_cvdraws(e_senc_average)
cv_draws <- ei_to_cvdraws(e_senc_contextual)
cv_draws <- ei_to_cvdraws(e_senc_bb)
cv_draws <-ei_cv_summary(senc_contextualonion)
cv_draws <-ei_cv_summary(senc_contextual)
cv_draws_a <- ei_cv_summary(senc_average)|>
  left_join(senc2_flat, by = c("area_no", "col_no"="col", "row_no"="row"))
cv_summary <- ssEI::ei_cv_summary(senc_contextual)|>
  left_join(senc2_flat, by = c("area_no", "col_no"="col", "row_no"="row"))
cv_draws_oo <-ei_cv_summary(senc_contextualonion) |>
  left_join(senc2_flat, by = c("area_no", "col_no"="col", "row_no"="row"))
rr_summary<- ei_row_rate_summary(senc_contextual)|>
  left_join(senc2_flat, by = c("area_no", "col_no"="col", "row_no"="row"))

rr_summary |>
  filter(mean> -1) |>
  ggplot(aes(actual_row_prop, mean, colour=paste(row_names, col_names))) +
  geom_point(aes(size=actual_value))


rr_summary |>
  filter(mean> -1) |>
  filter(actual_row_total>10) |>
  ggplot(aes(actual_row_prop, mean, colour=paste(row_names, col_names))) +
  facet_grid(rows = vars(row_names), cols = vars(col_names))+
  geom_point(aes(size=actual_value))

rr_summary |>
  filter(actual_row_total>10) |>
  group_by(row_names, col_names) |>
  summarise(cor = cor(mean, actual_row_prop, use="pairwise.complete.obs"))

rr_summary |>
  summarise(cor = cor(mean, actual_row_prop, use="pairwise.complete.obs"))

# rr_summary |>
#   filter(is.na(actual_row_prop)) |>
#   select(mean, actual_row_total) |> View()

rr_summary |>
  ggplot(aes(actual_row_prop)) +
  facet_wrap(vars(row_names)) +
  geom_histogram()

rr_summary |>
  summarise(cor = cor(mean, actual_row_prop, use="pairwise.complete.obs"))

ei_error_index(cv_summary$mean, cv_summary$actual_value)
cor(cv_summary$mean, cv_summary$actual_value)
cor(rr_summary$mean, rr_summary$actual_row_prop, use="complete.obs")
ei_error_index(cv_draws_oo$mean, cv_draws_oo$actual_value)
ei_error_index(cv_draws_a$mean, cv_draws_a$actual_value)
rmse(cv_draws_oo$mean, cv_draws_oo$actual_value)
rmse(cv_draws_a$mean, cv_draws_a$actual_value)
mae(cv_draws_oo$mean, cv_draws_oo$actual_value)
mae(cv_draws_a$mean, cv_draws_a$actual_value)


rp_draws <-ei_row_rate_summary(senc_contextualonion)
rp_draws <- cvdraws_to_rowprops(cv_draws)
rp_summary <- summarize_cvdraws(rp_draws)

rp_draws <-ei_row_rate_summary(senc_contextualonion)
rp_draws <- cvdraws_to_rowprops(cv_draws)
rp_summary <- summarize_cvdraws(rp_draws)

cv_draws |>
  left_join(senc2_flat, by = c("area_no", "col_no"="col", "row_no"="row")) |>
  ggplot(aes(mean, actual_value))+
  geom_point()

rp_draws |>
  filter(mean>-1) |>
  left_join(senc2_flat, by = c("area_no", "col_no"="col", "row_no"="row")) |>
  ggplot(aes(actual_row_prop, mean))+
  geom_point(aes(size=actual_value))

rp_draws |>
  filter(mean>-1) |>
  left_join(senc2_flat, by = c("area_no", "col_no"="col", "row_no"="row")) |>
  filter(mean>-1, actual_value>1) |>
  ggplot(aes(actual_row_prop, mean, colour=paste(row_names, col_names)))+
  geom_point()+
  geom_point(aes(size=actual_value)) +
  theme_bw() +
  geom_abline()

rp_draws |>
  left_join(senc2_flat, by = c("area_no", "col_no"="col", "row_no"="row")) |>
  filter(mean>-1) |>
  ggplot(aes(actual_row_prop, mean, colour=paste(row_names, col_names)))+
  facet_grid(row_names~col_names) +
  geom_point(aes(size=actual_value)) +
  theme_bw()

rp_draws |>
  filter(mean>-1) |>
  left_join(senc2_flat, by = c("area_no", "col_no"="col", "row_no"="row")) |>
  group_by(row_names, col_names) |>
  summarise(cor=cor(mean, actual_row_prop))

rp_sum_act <- rp_summary |> dplyr::left_join(senc_long)


library(ggplot2)

plt <- ggplot(rp_sum_act,
              aes(est_cell_value, actual_cell_value)) +
  geom_point()

rp_sum_act |> ggplot(aes(actual_row_prop, est_row_prop, colour = sencname,
                       size = actual_cell_value)) +
  geom_point()

plot(rp_sum_act$actual_row_prop, rp_sum_act$est_row_prop)
cor(rp_sum_act$actual_row_prop, rp_sum_act$est_row_prop, use = "pairwise.complete.obs")
cor(rp_sum_act$actual_cell_value, rp_sum_act$est_cell_value)


standata <- list(
  n_areas = nrow(senc_rm),
  R = ncol(senc_rm),
  C = ncol(senc_cm),
  row_margins = senc_rm,
  col_margins = senc_cm
)
aa <- sampling(stanmodels$ssMND,
               data=standata,
               chains = 1,
               iter = 100)
bb <- sampling(stanmodels$ssContextual,
               data=standata,
               cores = 4)

sscontextpath <- paste0(getwd(), "/inst/stan/ssContextual.stan")
bb <- rstan::stan(file = sscontextpath,
     data = standata,
     cores = 4)


##---summary model tests----

##----otherempirical-----
library(ei.Datasets)
data("ei_SCO_2007")

# ei_SCO_2007[[3]][[1]] |> View()

ei_SCO_2007 |> select(Number_of_district, District, Votes_to_candidates) |>
  unnest(Votes_to_candidates)

nosmall_ei_SCO_2007 <- merge_small_options(ei_SCO_2007, min.party = 10, min.candidate = 10)
nosmall_ei_SCO_2007 <- ei_SCO_2007
which_slice <- 25

test_cand <- nosmall_ei_SCO_2007 |>
  dplyr::slice(which_slice) |>
  dplyr::select(Votes_to_candidates) |>
  tidyr::unnest(Votes_to_candidates) |>
  dplyr::select(-Polling, -Address)

test_party <- nosmall_ei_SCO_2007 |>
  dplyr::slice(which_slice) |>
  dplyr::select(Number_of_district, District, Votes_to_parties) |>
  tidyr::unnest(Votes_to_parties) |>
  dplyr::select(-Number_of_district, -District, -Polling, -Address)

test_res <- nosmall_ei_SCO_2007 |>
  dplyr::slice(which_slice) |>
  dplyr::select(District, District_cross_votes) |>
  tidyr::unnest(District_cross_votes)

test_res2 <- nosmall_ei_SCO_2007 |>
  dplyr::slice(which_slice) |>
  dplyr::select(Number_of_district, District, District_cross_percentages) |>
  tidyr::unnest(District_cross_percentages)


scotest_R <- ncol(test_cand)
scotest_C <- ncol(test_party)
scotest_n_areas <- nrow(test_cand)
scotest_row_margins <- test_cand
scotest_col_margins <- test_party

test_sco_to_stan <- list(n_areas = scotest_n_areas,
                         R = scotest_R,
                         C = scotest_C,
                         row_margins = scotest_row_margins,
                         col_margins = scotest_col_margins
)

scotest_average <- ei_estimate(row_margins = scotest_row_margins,
            col_margins = scotest_col_margins,
            model = "average")

scotest_areaspec_1e <- stan(file = dirmultinom_areaspecific_1e,
                            data = test_sco_to_stan,
                            cores = 4)

dirmultinom_1a <- paste0(ei_stan_path, "/test_development/ei_dirichlet_multinom_1a.stan")
scotest_dirmultinom_1a <- stan(file = dirmultinom_1a,
                               data = test_sco_to_stan,
                               cores = 4)

e_scotest <- rstan::extract(scotest_areaspec_1e)
e_scotest <- rstan::extract(scotest_dirmultinom_1a)
e_scotest_av <- rstan::extract(scotest_average)
scotest_cv_draws_av <- ei_to_cvdraws(e_scotest_av)
scotest_rp_draws_av <- cvdraws_to_rowprops(scotest_cv_draws_av)
scotest_rp_summary_av <- summarize_cvdraws(scotest_rp_draws_av)

scocells1<- purrr::array_tree(e_scotest$cell_values, 1)
# the_names <- vector(length = scotest_R*soctest_C)
# for(r in 1:scotest_R){
#   for (c in 1:scotest_C){
#     the_names[(r-1)*scotest_C + c] <- paste(names(scotest_row_margins)[r], names(scotest_col_margins)[c], sep="_")
#   }
# }

row_index <- tibble(row_index = 1:((test_sco_to_stan$n_areas)*test_sco_to_stan$R),
                    row_names = rep(names(scotest_row_margins),
                                    test_sco_to_stan$n_areas),
                    area_no = rep(1:test_sco_to_stan$n_areas,
                                  each = test_sco_to_stan$R))
scodistsummary1<- map(scocells1, function(x) {
  x<-as.data.frame(x)
  names(x) <- names(scotest_col_margins)
  x <- x |> mutate(row_index = row_number())
  x
})
scodistdraws <- tibble("draws" = scodistsummary1) |>
  mutate(iter = row_number()) |>
  unnest(draws) |>
  left_join(row_index) |>
  pivot_longer(cols = names(test_sco_to_stan$col_margins), names_to = "col_names") |>
  group_by(row_names, col_names, iter) |>
  summarise(dist_est_value = sum(value))

scocompres <- scodistdraws |>
  group_by(row_names, col_names) |>
  summarise(dist_est_max90 = quantile(dist_est_value, probs =.95),
            dist_est_min90 = quantile(dist_est_value, probs = .05),
            dist_est_value = quantile(dist_est_value, probs = .5)) |>
  left_join(
    test_res |>
      pivot_longer(cols = -c(1,2), names_to = "row_names") |>
      rename("col_names" = 2)
  )
scocompres|>
  ggplot(aes(value, dist_est_value)) +
  geom_point()+
  geom_errorbar(aes(ymax = dist_est_max90, ymin=dist_est_min90))+
  geom_abline(intercept = 0, slope =1)


scocompres |>
  ungroup()|>
  mutate(error_index_component = abs(value - dist_est_value)) |>
  summarise(area_tot = sum(value),
            error_index = 50*sum(error_index_component)/area_tot)


sco.king.formula <- cbind(dem, rep, non) ~ cbind(white, black, natam)
dbuf = ei(king.formula, data=senc2)


tmp_res<-data.frame(
  av = c(5.68, 6.91, 9.85, 4.16, 9.90, 2.66, 8.19, 13.3, 4.80, 6.77, 3.15, 6.85, 14.4, 2.23, 3.33, 5.26, 4.25, 8.2, 5.67, 4.47, 15.8, 4.2, 3.38),
  co = c(4.99, 8.66, 9.86, 3.99, 11.8, 4.35, 8.57, 13.5, 7.78, 7.42, 3.79, 9.35, 14.3, 2.5, 3.45, 6.17, 1.44, 10.4, 5.39, 6.65, 11.0, 4.05, 3.06))
library(ggplot2)
ggplot(tmp_res, aes(av, co)) +
  geom_point()+
  geom_abline(intercept = 0, slope =1) +
  stat_smooth(method="lm")
mean(tmp_res$av)
mean(tmp_res$co)
t.test(tmp_res$av, tmp_res$co, paired=TRUE)

ggplot(tmp_res, aes(av)) +
  geom_histogram(bins =6)
ggplot(tmp_res, aes(co)) +
  geom_histogram(bins =6)

##----
R <- ncol(row_margins_eth_wide)
C <- ncol(col_margins_ed_wide)
n_areas <- nrow(row_margins_eth_wide)


standata <- list(
  n_areas = n_areas,
  R = R,
  C = C,
  row_margins = row_margins_eth_wide,
  col_margins = col_margins_ed_wide,
  lkj_param = 2
)
R <- ncol(senc_rm)
C <- ncol(senc_cm)
n_areas <- nrow(senc_rm)

standata <- list(
  n_areas = n_areas,
  R = R,
  C = C,
  row_margins = senc_rm,
  col_margins = senc_cm,
  lkj_param = 2
)

tmpres <- rstan::sampling(stanmodels$ssContextualOnion,
                data=standata,
                chains = 4,
                cores = 4)

newstancode2 <- rstan::get_stancode(senc_contextual2)
R <- ncol(senc_rm)
C <- ncol(senc_cm)
n_areas <- nrow(senc_rm)
my_lkj_param = 2
standata <- list(
  n_areas = n_areas,
  R = R,
  C = C,
  row_margins = senc_rm,
  col_margins = senc_cm,
  lkj_param = my_lkj_param
)
senc_contextual2 <- rstan::stan(model_code = newstancode2,
            data = standata,
            cores = 4,
            iter =10)


ssei.path <- "C:/Users/gidon/OneDrive/Documents/ssEI"

smod.path <- paste0(ssei.path, "/inst/stan/")

mnb.path <- paste0(smod.path, "ssMNB.stan")
mnb2.path <- paste0(smod.path, "ssMNB2.stan")
mnd.path <- paste0(smod.path, "ssMND.stan")
co.path <- paste0(smod.path, "ssContextualOnion.stan")
co.new2.path <- paste0(smod.path, "ssContextual2.stan")
co.new3.path <- paste0(smod.path, "ssContextual3.stan")

set.seed(123)
test_set <- ssEI::mk_sim_tables(50, 3, 3)


test_for_stan <- list(
  n_areas = nrow(test_set$row_margins),
  R = ncol(test_set$row_margins),
  C = ncol(test_set$col_margins),
  row_margins = test_set$row_margins,
  col_margins = test_set$col_margins,
  lkj_param = 2
)

test_mnb <- rstan::stan(file = mnb.path,
                        data = test_for_stan,
                        cores = 4,
                        control = list(
                          "adapt_delta" = .9)
)

test_mnd <- rstan::stan(file = mnd.path,
                        data = test_for_stan,
                        cores = 4)

st1 <- mods_summary(list("mnb" = test_mnb,
                         "mnd" = test_mnd), test_set$row_rates)

st1$cv_eval$cv_eval_stats
st1$rr_eval$rr_eval_stats
st1$eval_plots$p_cv
st1$eval_plots$p_rr
st1$eval_plots$pf_rr

mnb_cv <- ei_cv_summary(test_mnb)
mnd_cv <- ei_cv_summary(test_mnd)
mnb_rr <- ei_row_rate_summary(test_mnb)
mnd_rr <- ei_row_rate_summary(test_mnd)

mnb_cv |> left_join(test_set$row_rates) |>
  summarise(cor(mean, actual_cell_value))

mnd_cv |> left_join(test_set$row_rates) |>
  summarise(cor(mean, actual_cell_value))

mnb_rr |> left_join(test_set$row_rates) |>
  summarise(cor(mean, actual_row_rate))

mnd_rr |> left_join(test_set$row_rates) |>
  summarise(cor(mean, actual_row_rate))

senc_for_stan <- list(
  n_areas = nrow(senc_rm),
  R = ncol(senc_rm),
  C = ncol(senc_cm),
  row_margins = senc_rm,
  col_margins = senc_cm,
  lkj_param = 2
)

senc_mnb <- rstan::stan(file = mnb.path,
                        data = senc_for_stan,
                        cores = 4)

senc_co_new3 <- rstan::stan(file = co.new3.path,
                            data = senc_for_stan,
                            cores = 4,
                            iter = 5000)

senc_co_new2 <- rstan::stan(file = co.new2.path,
                           data = senc_for_stan,
                           cores = 4)
senc_co_new <- rstan::stan(file = co.path,
                           data = senc_for_stan,
                           cores = 4)

ssum <- mods_summary(list("mnb2"=senc_mnb,
                          "co_new" = senc_co_new,
                          "co_new2"=senc_co_new2,
                          "co_new3"=senc_co_new3), senc2_flat |> rename(col_no=col, row_no=row, actual_cell_value = actual_value, actual_row_rate = actual_row_prop))

ssum <- mods_summary(list("mnb2"=senc_mnb,
                          "co new"=senc_co_new2), senc2_flat |> rename(col_no=col, row_no=row, actual_cell_value = actual_value, actual_row_rate = actual_row_prop))

ssum$cv_eval$cv_eval_stats
ssum$rr_eval$rr_eval_stats
ssum$eval_plots$p_cv
ssum$eval_plots$p_rr + ylim(0,1) + geom_point(aes(colour=paste(row_names, col_names), size=actual_cell_value))
ssum$eval_plots$pf_rr+ ylim(0,1)
ssum$rr_eval$rr_eval |> group_by(model, row_no, col_no) |> summarise(cor=cor(mean, actual_row_rate, use="pairwise.complete.obs")) |> arrange(row_no, col_no, model) |> View()

b_senc_cv <- ei_cv_summary(senc_mnb)
b_senc_rr <- ei_row_rate_summary(senc_mnb)

b_senc_cv |> left_join(senc2_flat, by=c("area_no", "row_no"="row", "col_no"="col"))|>
  summarise(cor(mean, actual_value))

b_senc_rr |> left_join(senc2_flat, by=c("area_no", "row_no"="row", "col_no"="col"))|>
  summarise(cor(mean, actual_row_prop, use="complete.obs"))

b_senc_rr |> left_join(senc2_flat, by=c("area_no", "row_no"="row", "col_no"="col"))|>
  filter(mean>-1) |>
  ggplot(aes(actual_row_prop, mean, colour=paste(row_names, col_names)))+
  geom_point(aes(size=actual_value))

b_senc_rr |> left_join(senc2_flat, by=c("area_no", "row_no"="row", "col_no"="col"))|>
  filter(mean>-1) |>
  ggplot(aes(actual_row_prop, mean, colour=paste(col_names, row_names)))+
  facet_grid(rows=vars(row_names), cols= vars(col_names)) +
  geom_point(aes(size=actual_value))

b_senc_rr |> left_join(senc2_flat, by=c("area_no", "row_no"="row", "col_no"="col"))|>
  filter(mean>-1) |>
  group_by(row_names, col_names) |>
  summarise(cor(mean, actual_row_prop, use="complete.obs"))

  ggplot(aes(actual_row_prop, mean, colour=paste(col_names, row_names)))+
  facet_grid(rows=vars(row_names), cols= vars(col_names)) +
  geom_point(aes(size=actual_value))


ei.path<-switch (Sys.info()["nodename"],
                   "GID-HOME-LENOVO"="C:/Users/Gidon/Dropbox/Analysis with R/EI",
                   "GID-LAPTOP-LENO" = "C:/Users/gidon/Dropbox/Analysis with R/EI",
                   "DM-GIA-099"="C:/Users/dgi0gc/Dropbox/Analysis with R/EI")
cacheloc <- paste0(ei.path, "/model_cache")

# read the education ethnicity crosstable from 2021 UK Census
c21_ed_eth <- read_csv(paste0(ei.path, "/data/uk_census_2021/ENG_LTLA_ED7_ETH6.csv"
)) |>  rename(la_code = 1,
              la_name = 2,
              eth_code = 3,
              eth_group = 4,
              ed_code = 5,
              ed_group =6,
              actual_cell_value = 7
) |>
  group_by(la_code) |>
  mutate(area_no = cur_group_id())

## Note: there are no missing values for ethnicity so remove the -8 cases for
c21_ed_eth <- c21_ed_eth |> filter(eth_code != -8)

set.seed(123)
c21_ed_eth_small <- xfun::cache_rds({
  c21_ed_eth |> filter(
    la_code %in% sample(unique(c21_ed_eth$la_code), 100)
  ) },
  rerun = FALSE,
  dir = cacheloc,
  file = "c21_ed_eth_small.rds")

prep_dat <- function(x, row_code, col_code){
  x <- x |>
    group_by(area_no) |>
    mutate(new_area_no = cur_group_id()) |>
    ungroup() |>
    select(-area_no) |>
    rename(area_no=new_area_no)
  row_nos <- x |>
    group_by({{row_code}}) |>
    summarise() |>
    mutate(row_no = dplyr::row_number(),
           row_name = paste0("row_no.", row_no))

  col_nos <- x |>
    group_by({{col_code}}) |>
    summarise() |>
    mutate(col_no = dplyr::row_number(),
           col_name = paste0("col_no.", col_no))

  row_margins_long <- x |>
    group_by(area_no, {{row_code}}) |>
    summarise(actual_row_margin = sum(actual_cell_value)) |>
    ungroup() |>
    left_join(row_nos)

  row_margins_wide <- row_margins_long |>
    select(area_no, row_name, actual_row_margin) |>
    pivot_wider(names_from = row_name,
                values_from = "actual_row_margin")|>
    select(-area_no)

  col_margins_long <- x |>
    group_by(area_no, {{col_code}}) |>
    summarise(actual_col_margin = sum(actual_cell_value)) |>
    ungroup() |>
    left_join(col_nos)

  col_margins_wide <- col_margins_long |>
    select(area_no, col_name, actual_col_margin) |>
    pivot_wider(names_from = "col_name",
                values_from = "actual_col_margin") |>
    select(-area_no)

  cv_rr_long <- x |>
    left_join(col_nos) |>
    left_join(row_nos) |>
    left_join(row_margins_long |>
                select(area_no, row_no, actual_row_margin)) |>
    mutate(actual_row_rate = actual_cell_value/actual_row_margin)


  list(
    row_margins_long = row_margins_long,
    row_margins_wide = row_margins_wide,
    col_margins_long = col_margins_long,
    col_margins_wide = col_margins_wide,
    cv_rr_long = cv_rr_long,
    row_nos = row_nos,
    col_nos = col_nos
  )

}
set.seed(1234)
ndists <- seq(from = 20, to = 200, by =20)
csubsets <- vector("list", length = length(ndists))
for(a in 1:length(ndists)){
  s <- ndists[a]
  this_dat <- c21_ed_eth |>
    filter(la_code %in% sample(unique(c21_ed_eth$la_code), s))
  csubsets[[a]] <- prep_dat(this_dat, eth_code, ed_code)
}

csubsets <- xfun::cache_rds({
  csubsets
},
rerun = FALSE,
dir = cacheloc,
file = "csubsets.rds"
)
ces1_av <- xfun::cache_rds(
  {
    ei_estimate(row_margins = csubsets[[1]]$row_margins_wide,
                col_margins = csubsets[[1]]$col_margins_wide,
                model="average"
    )

  },
  rerun = FALSE,
  dir = cacheloc,
  file = "cs1_av.rds"
)

prep_for_stan <- function(rm, cm, lkj_param = 2){
  list(
    n_areas = nrow(rm),
    R = ncol(rm),
    C = ncol(cm),
    row_margins = rm,
    col_margins = cm,
    lkj_param = lkj_param
  )
}
cs1_dat <- prep_for_stan(csubsets[[1]]$row_margins_wide,
                       csubsets[[1]]$col_margins_wide
)

ces1_mnb1 <- xfun::cache_rds(
  {
    rstan::stan(file = mnb.path,
                data = cs1_dat,
                cores = 4)
  },
  rerun = FALSEpo,
  dir = cacheloc,
  file = "cs1_mnb1.rds"
)

ces1_mnb2 <- xfun::cache_rds(
  {
    rstan::stan(file = mnb2.path,
                data = cs1_dat,
                cores = 4)
  },
  rerun = FALSE,
  dir = cacheloc,
  file = "cs1_mnb2.rds"
)

ces1_co <- xfun::cache_rds(
  {
    ei_estimate(row_margins = csubsets[[1]]$row_margins_wide,
                col_margins = csubsets[[1]]$col_margins_wide,
                model="contextualOnion"
    )

  },
  rerun = FALSE,
  dir = cacheloc,
  file = "cs1_co.rds"
)

ces1_co_new <- xfun::cache_rds(
  {
    rstan::stan(file = co.path,
                data = cs1_dat,
                cores = 4)

  },
  rerun = FALSE,
  dir = cacheloc,
  file = "cs1_co_new.rds"
)

ces1_co_new2 <- xfun::cache_rds(
  {
    rstan::stan(file = co.path,
                data = cs1_dat,
                cores = 4)

  },
  rerun = FALSE,
  dir = cacheloc,
  file = "cs1_co_new2.rds"
)

ces1_co_new3 <- xfun::cache_rds(
  {
    ei_estimate(row_margins = csubsets[[1]]$row_margins_wide,
                col_margins = csubsets[[1]]$col_margins_wide,
                model="contextualOnion"
    )

  },
  rerun = FALSE,
  dir = cacheloc,
  file = "cs1_co_new2.rds"
)

r1 <- mods_summary(list("mnd"=ces1_av,
                        "mnb1"=ces1_mnb1,
                        "mnb2"=ces1_mnb2),
                   actual_long = csubsets[[1]]$cv_rr_long)

r1 <- mods_summary(list("co"=ces1_co,
                        "co_new"=ces1_co_new,
                        "co_new2"=ces1_co_new2),
                   actual_long = csubsets[[1]]$cv_rr_long)

r1$cv_eval$cv_eval_stats
r1$rr_eval$rr_eval_stats

r1$eval_plots$p_rr
r1$eval_plots$pf_rr

r1$eval_plots$p_cv


prep.king <- function(cm, rm){
  formula <- as.formula(paste0("cbind(", paste(names(cm), collapse = ","), ") ~ cbind(", paste(names(rm), collapse=","), ")"))
  list(
    formula = formula,
    data = cbind(cm, rm)
  )
}


kc1_dat <- prep_king(csubsets[[1]]$col_margins_wide, csubsets[[1]]$row_margins_wide)




kc1.tune.nocov <- tuneMD(kc1_dat$formula, data = kc1_dat$data,
                         ntunes = 10, totaldraws = 100000)


kc1.out.nocov <- ei.MD.bayes(kc1_dat$formula,
                             covariate = NULL,
                             data = kc1_dat$data,
                             tune.list = kc1.tune.nocov)

gqc1_dat <- prep_gq(csubsets[[1]]$col_margins_wide, csubsets[[1]]$row_margins_wide)
gq.tune.c1 <- Tune(gqc1_dat$formula,
                  data = gqc1_dat$data,
                  num.iters = 20000,
                  num.runs = 15)
niter.gq <- 1500000
Chain1.c1 <- Analyze(gqc1_dat$formula,
                       rho.vec = gq.tune.c1$rhos,
                       data = gqc1_dat$data,
                       num.iters = niter.gq,
                       burnin = 150000,
                       save.every = 1000,
                       debug = 1,
                       keepNNinternals = 100,
                       keepTHETAS = 100)
Chain2.c1 <- Analyze(gqc1_dat$formula,
                     rho.vec = gq.tune.c1$rhos,
                     data = gqc1_dat$data,
                     num.iters = niter.gq,
                     burnin = 150000,
                     save.every = 1000,
                     debug = 1,
                     keepNNinternals = 100,
                     keepTHETAS = 100)
Chain3.c1 <- Analyze(gqc1_dat$formula,
                     rho.vec = gq.tune.c1$rhos,
                     data = gqc1_dat$data,
                     num.iters = niter.gq,
                     burnin = 150000,
                     save.every = 1000,
                     debug = 1,
                     keepNNinternals = 100,
                     keepTHETAS = 100)
c1.MCMClist <- mcmc.list(Chain1.c1, Chain2.c1,Chain3.c1)
gelman.diag(c1.MCMClist, multivariate = FALSE)

c1s <- mods_summary(list("king"= kc1.out.nocov, "average"=ces1_av, "co new"=ces1_co_new ), actual_long =csubsets[[1]]$cv_rr_long)
c1s <- mods_summary(list("king"= kc1.out.nocov, gq = c1.MCMClist, "average"=ces1_av, "co new"=ces1_co_new, "co orig" = ces1_co ), actual_long =csubsets[[1]]$cv_rr_long)
c1s <- mods_summary(list(gq = c1.MCMClist, "co new"=ces1_co_new), actual_long =csubsets[[1]]$cv_rr_long)
c1s <- mods_summary(list("av"=ces1_av), actual_long =csubsets[[1]]$cv_rr_long)
c1s$cv_eval$cv_eval_stats
c1s$rr_eval$rr_eval_stats
c1s$eval_plots$p_cv
c1s$eval_plots$p_rr
c1s$eval_plots$pf_rr


cedeth_co <- xfun::cache_rds(
  {
    ei_estimate(row_margins = row_margins_eth_wide,
                col_margins = col_margins_ed_wide,
                model="contextualOnion"
    )

  },
  rerun = FALSE,
  dir = cacheloc,
  file = "cdeth_co.rds"
)


library(parallel)

gq.Analyze <- function(formula, data, tune, niter.gq, burnin.iter = 150000, save.every=1000, debug = 1){
    RxCEcolInf::Analyze(formula,
                       rho.vec = tune$rhos,
                       data = data,
                       num.iters = niter.gq,
                       burnin = burnin.iter,
                       save.every = save.every,
                       debug = 1,
                       keepNNinternals = 100,
                       keepTHETAS = 100)
}

gq.Analyze <- function(x, ...){
  RxCEcolInf::Analyze(...)
}

res.list <- list(1,1,1)
c1 <-makeCluster(3, methods = FALSE)
clusterExport(c1, c("gqc1_dat", "gq.tune.c1"))
res <-  parLapply(c1, res.list,
               fun = gq.Analyze,
               fstring = gqc1_dat$formula,
               rho.vec = gq.tune.c1$rhos,
               data = gqc1_dat$data,
               num.iters = 1000,
               burnin = 100,
               save.every = 10,
               debug = 1,
               keepNNinternals = 100,
               keepTHETAS = 100
         )
stopCluster(c1)


##----simpleregmethods----

# This model runs a regression obsC ~ Poisson(predC)
# Where obsC[j,c] is the observed column total of column c in area j
# and predC[j,c] = beta[c] * row_margins[j,] where beta[c] is a column specific
reg1.path <- paste0(smod.path, "eiOreg.stan")
reg1 <- rstan::stan(file = reg1.path,
                data = cs1_dat,
                cores = 4)

reg2.path <- paste0(smod.path, "eiOreg2.stan")
reg2 <- rstan::stan(file = reg2.path,
                    data = cs1_dat,
                    cores = 4)


cprops <-cs1_dat$col_margins/rowSums(cs1_dat$col_margins)
rprops <- cs1_dat$row_margins/rowSums(cs1_dat$row_margins)
props <- cbind(cprops, rprops)
props$ref_row <- props$row_no.4
props_dm <- props - colSums(props)/nrow(props)

props_dm <- apply(props, 2, function(x){x-mean(x)}) |> data.frame()

lm(col_no.1 ~ log(row_no.1/ref_row) + log(row_no.2/ref_row) + log(row_no.3/ref_row) + log(row_no.4/ref_row) + log(row_no.5/ref_row), props) |> summary()

lm(col_no.1 ~ 1 +row_no.1  + row_no.3 + row_no.4 + row_no.5, props_dm) |> summary()
lm(col_no.1 ~ row_no.2, props_dm) |> summary()

cs1_overall <- csubsets[[1]]$cv_rr_long |> group_by(row_no) |> mutate(row_total=sum(actual_cell_value)) |> group_by(col_no) |> mutate(col_total = sum(actual_cell_value)) |> group_by(row_no, row_total, col_no, col_total) |> summarise(actual_cell_value = sum(actual_cell_value)) |>
  mutate(row_rate=actual_cell_value/row_total)

cs1_overall |> View()

lm(I(arm::logit(col_no.1)) ~ I(log(row_no.1/row_no.4)) +I(log(row_no.2/row_no.4)) + I(log(row_no.3/row_no.4))+ I(log(row_no.5/row_no.4)), props) |> summary()

lm(I(arm::logit(col_no.1)) ~ I(log(row_no.3/row_no.4)), props) |> summary()
lm(I(arm::logit(col_no.1)) ~ I(log(row_no.3/row_no.5)), props) |> summary()


ces2_av <- xfun::cache_rds(
  {
    ei_estimate(row_margins = csubsets[[2]]$row_margins_wide,
                col_margins = csubsets[[2]]$col_margins_wide,
                model="average"
    )

  },
  rerun = FALSE,
  dir = cacheloc,
  file = "cs2_av.rds"
)
c2s <- mods_summary(list("average"=ces2_av), actual_long =csubsets[[2]]$cv_rr_long)
c2s$eval_plots$pf_rr


ces3_av <- xfun::cache_rds(
  {
    ei_estimate(row_margins = csubsets[[3]]$row_margins_wide,
                col_margins = csubsets[[3]]$col_margins_wide,
                model="average"
    )

  },
  rerun = FALSE,
  dir = cacheloc,
  file = "cs3_av.rds"
)
c3s <- mods_summary(list("average"=ces3_av), actual_long =csubsets[[3]]$cv_rr_long)
c3s$eval_plots$pf_rr
ces4_av <- xfun::cache_rds(
  {
    ei_estimate(row_margins = csubsets[[4]]$row_margins_wide,
                col_margins = csubsets[[4]]$col_margins_wide,
                model="average"
    )

  },
  rerun = FALSE,
  dir = cacheloc,
  file = "cs4_av.rds"
)
c4s <- mods_summary(list("average"=ces4_av), actual_long =csubsets[[4]]$cv_rr_long)
c4s$eval_plots$pf_rr

ces5_av <- xfun::cache_rds(
  {
    ei_estimate(row_margins = csubsets[[5]]$row_margins_wide,
                col_margins = csubsets[[5]]$col_margins_wide,
                model="average"
    )

  },
  rerun = FALSE,
  dir = cacheloc,
  file = "cs5_av.rds"
)
c5s <- mods_summary(list("average"=ces5_av), actual_long =csubsets[[5]]$cv_rr_long)
c5s$eval_plots$pf_rr

full_prep <-prep_dat(c21_ed_eth, eth_code, ed_code)

cprops <- full_prep$col_margins_wide/rowSums(full_prep$col_margins_wide)
rprops <- full_prep$row_margins_wide/rowSums(full_prep$row_margins_wide)
props <- cbind(cprops, rprops)

hist(props$row_no.5)


##----censusedsoc----

c21_ed_soc <- read_csv(paste0(ei.path, "/data/uk_census_2021/ENG_LTLA_ED7_SOC10.csv"
)) |>  rename(la_code = 1,
              la_name = 2,
              ed_code = 3,
              ed_group = 4,
              soc_code = 5,
              soc_group =6,
              actual_cell_value = 7
) |>
  group_by(la_code) |>
  mutate(area_no = cur_group_id())

## Note: missing values for education and Socio-economic classification so remove the -8 cases for both
c21_ed_soc <- c21_ed_soc |> filter(ed_code != -8, soc_code!=-8)

c21es_full <- prep_dat(c21_ed_soc, {soc_code}, {ed_code})


set.seed(1234)
ndists <- seq(from = 20, to = 200, by =20)
csubsets2 <- vector("list", length = length(ndists))
for(a in 1:length(ndists)){
  s <- ndists[a]
  this_dat <- c21_ed_soc |>
    filter(la_code %in% sample(unique(c21_ed_soc$la_code), s))
  csubsets2[[a]] <- prep_dat(this_dat, soc_code, ed_code)
}

cse1_av <- ei_estimate(csubsets2[[1]]$row_margins_wide,
                       csubsets2[[1]]$col_margins_wide,
                       model="average")
cse1s <- mods_summary(list("average"=cse1_av), actual_long =csubsets2[[1]]$cv_rr_long)
cse1s$eval_plots$p_cv
cse1s$eval_plots$pf_cv
cse1s$eval_plots$pf_rr

co.path <- paste0(smod.path, "ssContextualOnion.stan")

cse_ss_co_res <- list(length=10)
cse_ss_av_res <- list(length=10)
for(sset in 1:5){
  this_dat <- prep_for_stan(csubsets2[[sset]]$row_margins_wide,
                           csubsets2[[sset]]$col_margins_wide
  )
  cse_ss_co_res[[sset]] <- xfun::cache_rds(
    {
      rstan::stan(file = co.path,
                  data = this_dat,
                  cores = 4)
    },
    rerun = FALSE,
    dir = paste0(cacheloc, "/"),
    file = paste0("cse_", sset, "_co_new.rds")
  )
}
mnd.path <- paste0(smod.path, "ssMND.stan")
for(sset in 2:5){
  this_dat <- prep_for_stan(csubsets2[[sset]]$row_margins_wide,
                            csubsets2[[sset]]$col_margins_wide
  )
  cse_ss_av_res[[sset]] <- xfun::cache_rds(
    {
      rstan::stan(file = mnd.path,
                  data = this_dat,
                  cores = 4)
    },
    rerun = FALSE,
    dir = paste0(cacheloc, "/"),
    file = paste0("cse_", sset, "_av.rds")
  )
}

cse1s <- mods_summary(list("average"=cse1_av,
                           "co" = cse_ss_co_res[[1]]), actual_long =csubsets2[[1]]$cv_rr_long)


cse1s$cv_eval$cv_eval_stats
cse1s$rr_eval$rr_eval_stats

cse1s$eval_plots$pf_cv
cse1s$eval_plots$pf_rr
cse1sa <- mods_summary(list(
                           "co" = cse_ss_co_res[[1]]), actual_long =csubsets2[[1]]$cv_rr_long)

cse2s <- mods_summary(list("co" = cse_ss_co_res[[2]]), actual_long =csubsets2[[2]]$cv_rr_long)
cse3s <- mods_summary(list("co" = cse_ss_co_res[[3]]), actual_long =csubsets2[[3]]$cv_rr_long)
cse2s$rr_eval$rr_eval_stats
cse2s$cv_eval$cv_eval_stats
cse3s$rr_eval$rr_eval_stats
cse3s$cv_eval$cv_eval_stats

cse1sa$eval_plots$pf_rr + ggtitle("subset 1")
cse2s$eval_plots$pf_rr + ggtitle("subset 2")
cse3s$eval_plots$pf_rr + ggtitle("subset 3")

set.seed(11)
sim_dat11 <- mk_sim_tables(100, 3, 3,
                           scale_sd = -1,
                           sd_sd = .2,
                           contextual_effects = TRUE,
                           lkj_eta = exp(-100))




sim_dat11$for_stan <- prep_for_stan(sim_dat11$row_margins, sim_dat11$col_margins)
sim11_mnb <- rstan::stan(file = mnb.path,
                                    data = sim_dat11$for_stan,
                                    cores = 4)

sim11_co_new <- rstan::stan(file = co.path,
                         data = sim_dat11$for_stan,
                         cores = 4)

s11s <- mods_summary(list("mnb"=sim11_mnb, "co new"=sim11_co_new), sim_dat11$row_rates)

s11s$cv_eval$cv_eval_stats
s11s$rr_eval$rr_eval_stats
s11s$eval_plots$p_cv
s11s$eval_plots$p_rr
s11s$eval_plots$pf_rr


sim_dat11$Sigma

s11rrwide <- sim_dat11$row_rates |> mutate(c = "c") |> unite(cell, c, row_no, col_no) |> select(area_no, cell, actual_row_rate) |> pivot_wider(names_from = cell, values_from = actual_row_rate)
corrplot::corrplot(cor(s11rrwide |> select(-area_no)))

s11rrwide|>
  ggplot(aes(c_2_2, c_3_3)) + geom_point() + stat_smooth(method="lm")

set.seed(12)
sim_dat12 <- mk_sim_tables(100, 3, 3,
                           scale_sd = -1,
                           sd_sd = .2,
                           contextual_effects = TRUE,
                           lkj_eta = exp(-300))

s12rrwide <- sim_dat12$row_rates |> mutate(c = "c") |>
  unite(cell, c, row_no, col_no) |>
  filter(actual_cell_value>10) |>
  select(area_no, cell, actual_row_rate) |>
  pivot_wider(names_from = cell, values_from = actual_row_rate)
corrplot::corrplot(cor(s12rrwide |> select(-area_no), use="pairwise.complete.obs"), method = "number")
s12rrwide|>
  ggplot(aes(c_1_1, c_2_1)) + geom_point() + stat_smooth(method="lm")


actualOmegal <- senc2_flat |> select(area_no, row, col, actual_row_prop) |>
  # mutate(col = paste0(",", col,"]"), row = paste0("Omega[", row)) |>
  pivot_wider(names_from="col", values_from ="actual_row_prop") |>
  mutate(across(c(3:5), function(x) log(x/`1` ))) |>
  select(-`1`)|>
  pivot_longer(cols = c(`2`,`3`), names_to = "col") |>
  mutate(cell = ((as.numeric(row) - 1) * 2) + as.numeric(col)-1) |>
  select(-row, -col) |>
  pivot_wider(names_from ="cell", values_from= value) |>
  drop_na() |>
  select(-area_no) |>
  rowwise() |>
  filter(!any(is.infinite(c_across(where(is.numeric)))))|>
  cor(use="pairwise.complete.obs") |>
  data.frame()|>
  setNames(1:6) |>
  rownames_to_column("idx1") |>
  pivot_longer(cols = - idx1, names_to = "idx2") |>
  mutate(cell = paste0("Omega[", idx1, ",", idx2, "]"))

actualOmega <- pull(actualOmegal, value)
names(actualOmega) <- pull(actualOmegal, cell)

plot(actualOmega, bb_o)
plot(actualOmega, aa_o)

library(ggrepel)
data.frame(actual=actualOmega,
           old = aa_o,
           new = bb_o,
           cell = names(aa_o),
           idx1 = pull(actualOmegal, idx1),
           idx2 = pull(actualOmegal, idx2)) |>
  pivot_longer(cols = c(old, new)) |>
  filter(actual < 1)|>
  filter(idx1 <= idx2) |>
  ggplot(aes(actual, value, colour = name))+
  geom_point()+
  geom_text_repel(aes(label = cell), position="jitter") +
  geom_abline(intercept =0, slope =1)

actualsigmal <- senc2_flat |> select(area_no, row, col, actual_row_prop) |>
  # mutate(col = paste0(",", col,"]"), row = paste0("Omega[", row)) |>
  pivot_wider(names_from="col", values_from ="actual_row_prop") |>
  mutate(across(c(3:5), function(x) log(x/`1` ))) |>
  select(-`1`)|>
  pivot_longer(cols = c(`2`,`3`), names_to = "col") |>
  mutate(cell = ((as.numeric(row) - 1) * 2) + as.numeric(col)-1)|>
  drop_na() |>
  select(-area_no) |>
  rowwise() |>
  filter(!any(is.infinite(c_across(where(is.numeric))))) |>
  group_by(cell) |>
  summarise(sd = sd(value))

data.frame(actual=actualsigmal$sd, old = aa_s, new=bb_s[1:6]) |>
  ggplot()+
  geom_point(colour="red", aes(actual, old))+
  geom_point(colour="blue", aes(actual, new))+
  geom_abline(intercept=0, slope =1)+
  ylim(0,1.5)+xlim(0,1.5)


## ----testnegbinom-----

data("st_louis_census", package = "bayesjackman")

negbin.dat<- within(list(), {
  y = st_louis_census$i8094
  N <- length(y)}))
nbtest.path <- paste0 (smod.path, "test_negbinom.stan")
nbt <- rstan::stan(nbtest.path,
                   data = negbin.dat)
library(rstanarm)
negbin_fit2 <- rstanarm::stan_glm.nb(i8094 ~ 1, data = st_louis_census)
