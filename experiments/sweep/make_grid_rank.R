## Grid for rank-1 against full-rank sorting: Rscript experiments/sweep/make_grid_rank.R  (writes experiments/sweep/grid_rank.csv)
## Rank-1 datasets with sorting, all methods; full-rank datasets (and no-sorting ones) only need the reference methods, because their
## King and ssEI fits are already in results/results_long.csv (same seeds, identical data). assim 0 throughout; sigma 0.15.
g1 <- expand.grid(seed = 1:3, seg = c(0.4, 0.8, 1.6), assim = 0, sort = c(0.15, 0.3), assoc = c(0.3, 1), rank = 1L, sigma = 0.15,
                  method = c("refs", "ssei_est", "ssei_kingqi", "king"), stringsAsFactors = FALSE)
g0 <- expand.grid(seed = 1:3, seg = c(0.4, 0.8, 1.6), assim = 0, sort = c(0, 0.15, 0.3), assoc = c(0.3, 1), rank = 0L, sigma = 0.15,
                  method = "refs", stringsAsFactors = FALSE)
g <- rbind(g0, g1)[, c("assoc", "sort", "assim", "seg", "sigma", "rank", "seed", "method")]; g$task <- seq_len(nrow(g))
write.csv(g, "experiments/sweep/grid_rank.csv", row.names = FALSE); print(table(g$method, g$rank))
