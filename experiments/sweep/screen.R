## Margin-only diagnostics for each sweep dataset (regenerated from the seed), joined to the results
suppressMessages(source("experiments/sim_tables.R")); d <- read.csv("experiments/sweep/results/results_long.csv")
g <- unique(d[, c("assoc", "sort", "assim", "seg", "sigma", "seed")])
scr <- do.call(rbind, lapply(seq_len(nrow(g)), function(i) { p <- g[i, ]
  x <- sim_tables(J = 60, R = 4, C = 4, assoc = p$assoc, sort = p$sort, sigma = p$sigma, seg = p$seg, seed = p$seed, assim = p$assim)
  kc <- x$kc; rm <- apply(kc, 1:2, sum); cm <- apply(kc, c(1, 3), sum); N <- sum(kc); xs <- rm / rowSums(rm); ys <- cm / rowSums(cm); w <- rowSums(rm)
  B <- solve(crossprod(xs * sqrt(w)), crossprod(xs * w, ys))                     # unconstrained Goodman (weighted), R x C
  res <- ys - xs %*% B; r2 <- 1 - colSums(w * res^2) / colSums(w * sweep(ys, 2, colSums(w * ys) / sum(w))^2)
  ind <- array(0, dim(kc)); for (j in 1:60) ind[j, , ] <- outer(rm[j, ], cm[j, ]) / w[j]
  data.frame(p, out01 = mean(B < 0 | B > 1), out_mass = sum(pmax(-B, 0) + pmax(B - 1, 0)) / 4, r2 = mean(r2), seg_obs = mean(apply(log(xs), 2, sd)),
             indep_err = 50 * sum(abs(ind - kc)) / N) }))
W <- reshape(d[, c("assoc", "sort", "assim", "seg", "seed", "method", "err")], idvar = c("assoc", "sort", "assim", "seg", "seed"), timevar = "method", direction = "wide")
names(W) <- sub("^err\\.", "", names(W)); W <- merge(W, scr)
cat("check regeneration: max |indep error recomputed - Hamilton| =", max(abs(W$indep - W$indep_err)), "\n\n")
W$gain_kq <- W$ssei_est - W$ssei_kingqi                                          # > 0: the King-like set-up wins
cat("Mean diagnostics by sorting level:\n"); print(aggregate(cbind(out01, out_mass, r2, seg_obs, gain_kq) ~ sort, W, function(v) round(mean(v), 3)))
cat("\nBy sorting x segregation:\n"); print(aggregate(cbind(out01, out_mass, r2, gain_kq) ~ sort + seg, W, function(v) round(mean(v), 3)))
cat("\nCorrelation of gain (est - kingqi) with diagnostics:\n"); print(round(cor(W$gain_kq, W[, c("out01", "out_mass", "r2", "seg_obs")]), 2))
cat("\nRule: use kingqi when out_mass > t, else est. Mean error of the rule vs always-est / always-kingqi / King / best-of-two oracle:\n")
for (t in c(0, .02, .05, .1, .15, .2, .3)) { rule <- ifelse(W$out_mass > t, W$ssei_kingqi, W$ssei_est)
  cat(sprintf("  t %.2f: rule %.2f | est %.2f kingqi %.2f king %.2f | oracle %.2f | kingqi chosen %3.0f%%\n", t, mean(rule), mean(W$ssei_est), mean(W$ssei_kingqi), mean(W$king),
              mean(pmin(W$ssei_est, W$ssei_kingqi)), 100 * mean(W$out_mass > t))) }
saveRDS(W, "experiments/sweep/results/wide_with_screen.rds")
