## Does ei.MD.bayes do badly on the per-district benchmark because it is unconverged? Sampled districts, three settings:
##   default: eiPack defaults (no tuneMD; burn-in 1000, 1000 draws, thin 1)
##   manual:  Pavia & Romero (2023) settings (tuneMD 10 x 100000; burn-in 100000; 1000 draws, thin 100), as bm_eiMD defaults
##   long:    manual, then extended (more draws at thin 100) until converged or `max_ext` extensions
## 2 chains each. Convergence: split R-hat and bulk ESS of the district-table cells (summed over polling places) and the alphas.
## usage: Rscript experiments/pavia_bench/king_conv.R <dataset> <n districts> <out.rds> [max_ext]
a <- commandArgs(TRUE); ds <- a[1]; nd <- as.integer(a[2]); out <- a[3]; max_ext <- if (length(a) >= 4) as.integer(a[4]) else 0
.libPaths(c("/tmp/claude-0/-home-claude/e71a6b8d-b589-5b64-8f29-846a7938884f/scratchpad/rlib", .libPaths()))
suppressMessages({ source("experiments/benchmark/bm_lib.R"); source("experiments/pavia_bench/instances.R"); library(fastei) })
ins <- pavia_instances(ds); set.seed(20261101); pick <- sort(sample(length(ins), nd))
rh <- function(x) {   # split R-hat, x = iterations x chains
  n <- floor(nrow(x) / 2); y <- cbind(x[1:n, , drop = FALSE], x[(n + 1):(2 * n), , drop = FALSE])
  B <- n * var(colMeans(y)); W <- mean(apply(y, 2, var)); sqrt(((n - 1) / n * W + B / n) / W) }
tot_conv <- function(f) {   # R-hat of the district-table cells
  D <- f$draws; ch <- f$chain; tot <- apply(D, c(1, 3, 4), sum)
  r <- apply(tot, 2:3, function(v) { m <- sapply(sort(unique(ch)), function(k) v[ch == k]); if (sd(v) == 0) 1 else rh(m) }); max(r) }
res <- list()
for (p in pick) {
  i <- ins[[p]]; row <- list(district = i$name, units = nrow(i$rm), R = ncol(i$rm), C = ncol(i$cm))
  m <- run_em(X = i$cm, W = i$rm, method = "mvn_pdf", symmetric = TRUE); row$fastei <- ei_index(apply(m$expected_outcome, 1:2, sum), i$truth)
  ## default eiPack settings, 2 chains
  t0 <- Sys.time(); d <- bm_codes(i$rm, i$cm)
  runs <- parallel::mclapply(1:2, function(k) { set.seed(k); eiPack::ei.MD.bayes(d$formula, data = d$data, sample = 1000, thin = 1, burnin = 1000) }, mc.cores = 2)
  ests <- sapply(runs, function(r) { B <- as.matrix(r$draws$Beta); mb <- colMeans(B)
    m <- regmatches(names(mb), regexec("^beta\\.r([0-9]+)\\.c([0-9]+)\\.([0-9]+)$", names(mb))); ix <- t(vapply(m, function(x) as.integer(x[2:4]), integer(3)))
    T <- matrix(0, ncol(i$rm), ncol(i$cm)); for (q in seq_along(mb)) T[ix[q, 1], ix[q, 2]] <- T[ix[q, 1], ix[q, 2]] + mb[q] * d$rm[ix[q, 3], ix[q, 1]]; as.vector(T) })
  row$default <- ei_index(matrix(rowMeans(ests), ncol(i$rm)), i$truth); row$default_chain_gap <- ei_index(matrix(ests[, 1], ncol(i$rm)), matrix(ests[, 2], ncol(i$rm)))
  row$default_secs <- as.numeric(Sys.time() - t0, units = "secs")
  ## manual settings, then extensions
  f <- bm_eiMD(i$rm, i$cm, chains = 2, cores = 2, max_extensions = max_ext)
  row$manual_secs <- f$secs; row$segments <- nrow(f$attempts)
  if (nrow(f$attempts) > 1) {   # error after the first segment only (= the manual settings), from its draws
    n1 <- f$attempts$draws_per_chain[1]; keep <- unlist(lapply(split(seq_along(f$chain), f$chain), function(ix) ix[1:n1]))
    row$manual <- ei_index(apply(f$draws[keep, , , , drop = FALSE], 3:4, sum) / length(keep), i$truth)
  } else row$manual <- ei_index(apply(f$est, 2:3, sum), i$truth)
  row$manual_rhat_cells <- f$attempts$max_rhat[1]; row$final <- ei_index(apply(f$est, 2:3, sum), i$truth)
  row$final_rhat_all <- f$conv$max_rhat; row$final_rhat_table <- tot_conv(f); row$final_ess <- f$conv$min_ess_bulk
  res[[length(res) + 1]] <- as.data.frame(row); saveRDS(do.call(rbind, res), out)
  cat(sprintf("%-28s units %3d %dx%d | fastei %.2f | default %.2f (chain gap %.1f, %.0fs) | manual %.2f (R-hat %.2f) | final %.2f (R-hat all %.2f, table %.2f, segments %d, %.0f min)\n",
              i$name, row$units, row$R, row$C, row$fastei, row$default, row$default_chain_gap, row$default_secs, row$manual, row$manual_rhat_cells,
              row$final, row$final_rhat_all, row$final_rhat_table, row$segments, row$manual_secs / 60))
}
