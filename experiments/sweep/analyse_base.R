## Low-dimensional base (ssei_fixed, ssei_scaled) against the first sweep's methods on the same 108 datasets
d <- rbind(read.csv("experiments/sweep/results/results_long.csv"), read.csv("experiments/sweep/results/base_long.csv")[, names(read.csv("experiments/sweep/results/results_long.csv"))])
bad <- d[d$method == "ssei_fixed" & !is.na(d$rhat) & d$rhat > 1.1, ]; cat("ssei_fixed fits with Rhat > 1.1:\n"); print(bad[, c("assoc", "sort", "assim", "seg", "seed", "err", "rhat")], row.names = FALSE)
meth <- c("goodman_raked", "king", "ssei_est", "ssei_kingqi", "ssei_fixed", "ssei_scaled", "rake_true")
W <- reshape(d[, c("assoc", "sort", "assim", "seg", "seed", "method", "err")], idvar = c("assoc", "sort", "assim", "seg", "seed"), timevar = "method", direction = "wide"); names(W) <- sub("^err\\.", "", names(W))
Wn <- reshape(d[, c("assoc", "sort", "assim", "seg", "seed", "method", "nat")], idvar = c("assoc", "sort", "assim", "seg", "seed"), timevar = "method", direction = "wide"); names(Wn) <- sub("^nat\\.", "", names(Wn))
dc <- d[!is.na(d$cover), ]
for (fac in c("sort", "seg", "assoc", "assim")) { cat("\nMean area error (national) by", fac, "\n")
  m <- aggregate(W[, meth], list(lev = W[[fac]]), mean); mn <- aggregate(Wn[, meth], list(lev = Wn[[fac]]), mean)
  for (i in seq_len(nrow(m))) cat(sprintf("  %-5s", m$lev[i]), sprintf("%s %5.2f (%5.2f)", substr(meth, 1, 10), unlist(m[i, -1]), unlist(mn[i, -1])), "\n") }
cat("\nOverall:\n"); print(round(rbind(area = colMeans(W[, meth]), national = colMeans(Wn[, meth])), 2))
pc <- function(x, y) cat(sprintf("  %-12s vs %-13s mean diff %6.2f, first better %3.0f%% | no sorting %6.2f (%3.0f%%) | sorting %6.2f (%3.0f%%)\n", x, y,
  mean(W[[x]] - W[[y]]), 100 * mean(W[[x]] < W[[y]]), mean((W[[x]] - W[[y]])[W$sort == 0]), 100 * mean((W[[x]] < W[[y]])[W$sort == 0]),
  mean((W[[x]] - W[[y]])[W$sort > 0]), 100 * mean((W[[x]] < W[[y]])[W$sort > 0])))
cat("\nPaired comparisons (108 datasets):\n")
pc("ssei_fixed", "ssei_est"); pc("ssei_scaled", "ssei_est"); pc("ssei_fixed", "ssei_scaled"); pc("ssei_fixed", "goodman_raked"); pc("ssei_fixed", "king"); pc("ssei_scaled", "king")
cat("\nNo sorting, by segregation x association (area error):\n"); print(aggregate(W[W$sort == 0, meth], list(seg = W$seg[W$sort == 0], assoc = W$assoc[W$sort == 0]), function(v) round(mean(v), 2)))
cat("\nCoverage (area cells 50+ / national) by sorting:\n"); print(aggregate(cbind(cover, cover_nat) ~ method + sort, dc[dc$method %in% c("king", "ssei_est", "ssei_kingqi", "ssei_fixed", "ssei_scaled"), ], function(v) round(mean(v), 2)))
cat("\nMinutes per fit (median):\n"); print(aggregate(mins ~ method, d[d$method %in% meth, ], median))
