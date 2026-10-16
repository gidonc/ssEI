## Reference methods only, sorting linear in the shares (not the log shares): Rscript experiments/sweep/make_grid_form.R -> grid_form.csv
g <- expand.grid(seed = 1:3, seg = c(0.4, 0.8, 1.6), assim = 0, sort = c(0.15, 0.3), assoc = c(0.3, 1), rank = c(0L, 1L), form = "share",
                 sigma = 0.15, method = "refs", stringsAsFactors = FALSE)
g <- g[, c("assoc", "sort", "assim", "seg", "sigma", "rank", "form", "seed", "method")]; g$task <- seq_len(nrow(g))
write.csv(g, "experiments/sweep/grid_form.csv", row.names = FALSE); cat(nrow(g), "tasks\n")
