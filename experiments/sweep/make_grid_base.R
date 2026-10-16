## Low-dimensional base on the main sweep's 108 datasets: Rscript experiments/sweep/make_grid_base.R -> experiments/sweep/grid_base.csv
## ssei_fixed (interior fixed: E_rc raked to each area's margins) and ssei_scaled (interior scaled by an estimated kappa).
## The other methods on the same datasets are in results/results_long.csv (identical data).
g <- expand.grid(seed = 1:3, seg = c(0.4, 0.8, 1.6), assim = c(0, 1), sort = c(0, 0.15, 0.3), assoc = c(0.3, 1), sigma = 0.15,
                 method = c("ssei_fixed", "ssei_scaled"), stringsAsFactors = FALSE)
g <- g[, c("assoc", "sort", "assim", "seg", "sigma", "seed", "method")]; g$task <- seq_len(nrow(g))
write.csv(g, "experiments/sweep/grid_base.csv", row.names = FALSE); cat(nrow(g), "tasks\n")
