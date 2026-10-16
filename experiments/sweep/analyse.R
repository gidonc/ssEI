## Analysis of the simulation sweep: Rscript experiments/sweep/analyse.R [results_long.csv]
a <- commandArgs(TRUE); f <- if (length(a)) a[1] else "experiments/sweep/results/results_long.csv"
d <- read.csv(f); options(width = 200)
meth <- c("indep", "goodman_raked", "king", "ssei_est", "ssei_s08", "ssei_kingqi"); d$cell <- with(d, paste(assoc, sort, assim, seg))
W <- reshape(d[, c("assoc", "sort", "assim", "seg", "seed", "cell", "method", "err")], idvar = c("assoc", "sort", "assim", "seg", "seed", "cell"), timevar = "method", direction = "wide")
names(W) <- sub("^err\\.", "", names(W))
Wn <- reshape(d[, c("assoc", "sort", "assim", "seg", "seed", "method", "nat")], idvar = c("assoc", "sort", "assim", "seg", "seed"), timevar = "method", direction = "wide"); names(Wn) <- sub("^nat\\.", "", names(Wn))
cat("== Mean area error (national error) by factor level, averaged over everything else ==\n")
for (fac in c("assoc", "sort", "assim", "seg")) { cat("--", fac, "\n")
  m <- aggregate(W[, c(meth, "rake_true")], list(lev = W[[fac]]), mean); mn <- aggregate(Wn[, c(meth, "rake_true")], list(lev = Wn[[fac]]), mean)
  for (i in seq_len(nrow(m))) cat(sprintf("  %-5s", m$lev[i]), sprintf("%s %5.2f (%5.2f)", substr(c(meth, "rake_true"), 1, 9), unlist(m[i, -1]), unlist(mn[i, -1])), "\n") }
cat("\n== Overall mean area / national error ==\n"); print(round(rbind(area = colMeans(W[, c(meth, "rake_true")]), national = colMeans(Wn[, c(meth, "rake_true")])), 2))
cat("\n== Best feasible method per design cell (mean over 3 seeds); 36 cells ==\n")
M <- aggregate(W[, meth], list(assoc = W$assoc, sort = W$sort, assim = W$assim, seg = W$seg), mean)
best <- meth[apply(M[, meth], 1, which.min)]; print(table(best))
cat("best by sorting level:\n"); print(table(M$sort, best))
cat("\n== Paired comparisons (per dataset, n = 108): mean difference, share of datasets where the first is better ==\n")
pc <- function(x, y) cat(sprintf("  %-14s vs %-14s mean diff %6.2f | first better in %3.0f%% | within sort 0: %6.2f (%3.0f%%), sort>0: %6.2f (%3.0f%%)\n", x, y,
  mean(W[[x]] - W[[y]]), 100 * mean(W[[x]] < W[[y]]), mean((W[[x]] - W[[y]])[W$sort == 0]), 100 * mean((W[[x]] < W[[y]])[W$sort == 0]),
  mean((W[[x]] - W[[y]])[W$sort > 0]), 100 * mean((W[[x]] < W[[y]])[W$sort > 0])))
pc("ssei_est", "goodman_raked"); pc("ssei_est", "king"); pc("ssei_kingqi", "king"); pc("ssei_s08", "king"); pc("ssei_kingqi", "ssei_s08"); pc("ssei_kingqi", "ssei_est"); pc("king", "goodman_raked")
cat("\n== Seed-to-seed spread: sd over seeds within a cell, averaged over cells ==\n")
print(round(sapply(meth, function(m) mean(tapply(W[[m]], W$cell, sd))), 2))
cat("\n== Coverage (area cells 50+ / national table), by sorting and segregation ==\n")
dc <- d[!is.na(d$cover), ]
print(aggregate(cbind(cover, cover_nat) ~ method + sort, dc, function(v) round(mean(v), 2)))
print(aggregate(cbind(cover, cover_nat) ~ method + seg, dc, function(v) round(mean(v), 2)))
cat("\n== ssEI E_rc Rhat: share > 1.05, max ==\n"); print(aggregate(rhat ~ method, dc[!is.na(dc$rhat), ], function(v) round(c(gt105 = mean(v > 1.05), max = max(v)), 3)))
cat("\n== Where ssei_s08 / kingqi blow up: area error > 15 with no sorting ==\n"); print(d[d$sort == 0 & d$err > 15 & d$method %in% meth[-1], c("assoc", "assim", "seg", "seed", "method", "err", "nat", "cover", "rhat")], row.names = FALSE)
