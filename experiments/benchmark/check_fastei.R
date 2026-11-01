## fastei (Ubilla Pavez, Hermosilla & Thraves 2026, Statistics and Computing): EM for one probability matrix shared by all units,
## run on the same set-ups as ssEI (units = constituencies / electorates / precincts; truth known).
## Area error uses fastei's own per-unit expected counts (E-step under the fitted matrix, given each unit's margins), and also the
## fitted matrix raked to each unit; national error from the summed per-unit estimates.
## usage: Rscript experiments/benchmark/check_fastei.R [path to a library holding fastei]
a <- commandArgs(TRUE); if (length(a)) .libPaths(c(a[1], .libPaths()))
suppressMessages({ library(fastei); source("experiments/data_prep.R"); source("experiments/benchmark/bm_data.R") })
ipf <- function(X, rt, ct, it = 2000) { for (i in 1:it) { rs <- rowSums(X); X <- X * ifelse(rs > 0, rt / pmax(rs, 1e-300), 0); cs <- colSums(X)
  X <- sweep(X, 2, ifelse(cs > 0, ct / pmax(cs, 1e-300), 0), "*"); if (i %% 25 == 0 && max(abs(rowSums(X) - rt)) < 1e-7) break }; X }
dsets <- list(senc = bm_data("senc")$kc, scotland = prep_scotland(7, 73)$kc, nz2020 = prep_nz(2020, 6, "candidate", 1e6)$kc)
for (nm in names(dsets)) {
  kc <- dsets[[nm]]; J <- dim(kc)[1]; N <- sum(kc); kt <- apply(kc, 2:3, sum)
  W <- apply(kc, 1:2, sum); X <- apply(kc, c(1, 3), sum)          # W: units x groups (rows), X: units x candidates (columns)
  cat(sprintf("\n%s (%d units, %d x %d): area (own E-step) | area (matrix raked) | national\n", nm, J, dim(kc)[2], dim(kc)[3]))
  for (cfg in list(list("mult", FALSE), list("mvn_pdf", FALSE), list("mvn_cdf", FALSE), list("mvn_pdf", TRUE), list("mult", TRUE))) {
    t0 <- Sys.time()
    m <- tryCatch(run_em(X = X, W = W, method = cfg[[1]], symmetric = cfg[[2]], maxtime = 600), error = function(e) { cat("  ", cfg[[1]], "error:", conditionMessage(e), "\n"); NULL })
    if (is.null(m)) next
    eo <- aperm(m$expected_outcome, c(3, 1, 2))                       # units x rows x cols
    P <- m$prob; rk <- array(0, dim(kc)); for (j in 1:J) rk[j, , ] <- ipf(pmax(P, 1e-9) * W[j, ], W[j, ], X[j, ])
    cat(sprintf("  %-8s%-5s %6.2f | %6.2f | %5.2f   (margin misfit %.1f votes; %.1fs)\n", cfg[[1]], if (cfg[[2]]) "sym" else "", 50 * sum(abs(eo - kc)) / N,
                50 * sum(abs(rk - kc)) / N, 50 * sum(abs(apply(eo, 2:3, sum) - kt)) / N, max(abs(apply(eo, c(1, 3), sum) - X)), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
}
