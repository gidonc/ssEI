library(ei)
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

senc_average <- ei_estimate(senc_rm, senc_cm, model="average")
senc_contextual <- ei_estimate(senc_rm, senc_cm)
senc_contextualonion <- ei_estimate(senc_rm, senc_cm, model="contextualOnion")
e_senc_average <- rstan::extract(senc_average)
e_senc_bb <- rstan::extract(bb)
e_senc_context <- rstan::extract(senc_contextual)
cv_draws <- ei_to_cvdraws(e_senc_average)
cv_draws <- ei_to_cvdraws(e_senc_context)
cv_draws <- ei_to_cvdraws(e_senc_bb)
cv_draws <-ei_cv_summary(senc_contextualonion)
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
