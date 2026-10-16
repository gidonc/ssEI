## Writes experiments/sweep/grid.csv: one row per (dataset, method) task for the Slurm array.
## Datasets: assoc x sort x assim x seg x seed (sigma fixed at 0.15, 60 areas, 4 x 4). Edit the levels here to change the sweep.
lv <- expand.grid(seed = 1:3, seg = c(0.4, 0.8, 1.6), assim = c(0, 1), sort = c(0, 0.15, 0.3), assoc = c(0.3, 1), sigma = 0.15,
                  method = c("refs", "ssei_est", "ssei_s08", "ssei_kingqi", "king"), stringsAsFactors = FALSE)
lv <- lv[, c("assoc", "sort", "assim", "seg", "sigma", "seed", "method")]
lv$task <- seq_len(nrow(lv))
write.csv(lv, "experiments/sweep/grid.csv", row.names = FALSE)
cat(nrow(lv), "tasks;", nrow(unique(lv[, 1:6])), "datasets x", length(unique(lv$method)), "methods\n")
