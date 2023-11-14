library(ei)
library(ggplot2)


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
e_senc_average <- rstan::extract(senc_average)
e_senc_bb <- rstan::extract(bb)
e_senc_context <- rstan::extract(senc_contextual)
cv_draws <- ei_to_cvdraws(e_senc_average)
cv_draws <- ei_to_cvdraws(e_senc_context)
cv_draws <- ei_to_cvdraws(e_senc_bb)
rp_draws <- cvdraws_to_rowprops(cv_draws)
rp_summary <- summarize_cvdraws(rp_draws)

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

rstan::Rhat(e_senc_context$L)

pcf <- c(0,4,8,12,16,20,24,28,32,36,40,44,48,52,56,60,64,68,72,76,80,84,88,92,96,100,104,108,112,116,120,124,128,132,136,140,144,148,152,156,160,164,168,172,176,180,184,188,192,196,200,204,208,212,216,220,224,228,232,236,240,244,248,252,256,260,264,268,272,276,280,284,288,292,296,300,304,308,312,316,320,324,328,332,336,340,344,348,352,356,360,364,368,372,376,380,384,388,392,396,400,404,408,412,416,420,424,428,432,436,440,444,448,452,456,460,464,468,472,476,480,484,488,492,496,500,504,508,512,516,520,524,528,532,536,540,544,548,552,556,560,564,568,572,576,580,584,588,592,596,600,604,608,612,616,620,624,628,632,636,640,644,648,652,656,660,664,668,672,676,680,684,688,692,696,700,704,708,712,716,720,724,728)

res <- vector(length = 0)
for (j in 1:n_areas){
     # lambda[j] = rep_array(0, R - 1, C - 1);
     for (r in 1:2){
       for (c in 1:2){
        res <- c(res, pcf[j] + (r - 1) * 2 + c)
     }
   }
 }
res - 1:(n_areas*2*2)

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
