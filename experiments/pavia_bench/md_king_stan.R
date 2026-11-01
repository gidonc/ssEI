## Converged answers for the ei.MD.bayes model on the per-district benchmark: the same model (eiPack default prior) sampled by HMC.
## usage: Rscript experiments/pavia_bench/md_king_stan.R <dataset> <n districts> <out.rds> [chains] [iter] [which of the sampled districts, e.g. 1,2] [adapt_delta]
a <- commandArgs(TRUE); ds <- a[1]; nd <- as.integer(a[2]); out <- a[3]; ch <- if (length(a) >= 4) as.integer(a[4]) else 2; it <- if (length(a) >= 5) as.integer(a[5]) else 2000
suppressMessages({ library(rstan); source("experiments/pavia_bench/instances.R") }); rstan_options(auto_write = TRUE)
sm <- stan_model("experiments/pavia_bench/md_king.stan")
ins <- pavia_instances(ds); set.seed(20261101); pick <- sort(sample(length(ins), nd))
only <- if (length(a) >= 6) as.integer(strsplit(a[6], ",")[[1]]) else seq_along(pick); pick <- pick[only]
ad <- if (length(a) >= 7) as.numeric(a[7]) else 0.9
res <- list()
for (p in pick) {
  i <- ins[[p]]; J <- nrow(i$rm); R <- ncol(i$rm); C <- ncol(i$cm)
  dat <- list(J = J, R = R, C = C, X = i$rm / rowSums(i$rm), Y = round(i$cm), lambda1 = 4, lambda2 = 2)
  t0 <- Sys.time()
  fit <- sampling(sm, data = dat, chains = ch, cores = ch, iter = it, warmup = it / 2, seed = 1, refresh = 0, control = list(adapt_delta = ad))
  B <- rstan::extract(fit, "beta")$beta                                    # draws x J x R x C
  tab <- apply(sweep(B, c(2, 3), i$rm, "*"), c(1, 3, 4), sum)             # district table per draw
  est <- apply(tab, 2:3, mean)
  ## R-hat / ESS of the district-table cells and the alphas
  arr <- as.array(fit, pars = "alpha"); tabc <- array(NA, c(dim(arr)[1], ch, R * C))
  perm <- rstan::extract(fit, "beta", permuted = FALSE)                     # iter x chain x params
  for (k in 1:ch) { Bk <- array(perm[, k, ], c(dim(perm)[1], J, R, C)); tabc[, k, ] <- matrix(apply(sweep(Bk, c(2, 3), i$rm, "*"), c(1, 3, 4), sum), dim(perm)[1]) }
  ei_chain <- sapply(1:ch, function(k) ei_index(matrix(colMeans(tabc[, k, ]), R, C), i$truth))
  saveRDS(list(tabc = tabc, truth = i$truth, alpha = arr), sub("\\.rds$", paste0("_", p, "_draws.rds"), out))
  rh <- max(apply(tabc, 3, function(x) rstan::Rhat(x)), na.rm = TRUE); es <- min(apply(tabc, 3, function(x) rstan::ess_bulk(x)), na.rm = TRUE)
  rha <- max(apply(arr, 3, rstan::Rhat))
  sp <- do.call(rbind, get_sampler_params(fit, inc_warmup = FALSE))
  row <- data.frame(district = i$name, units = J, R = R, C = C, ei = ei_index(est, i$truth), rhat_table = rh, ess_table = es, rhat_alpha = rha,
                    div = sum(sp[, "divergent__"]), ei_chains = paste(sprintf("%.2f", ei_chain), collapse = " "), mins = as.numeric(Sys.time() - t0, units = "mins"))
  res[[length(res) + 1]] <- row; saveRDS(do.call(rbind, res), out)
  cat(sprintf("%-28s units %3d %dx%d | MD (HMC) EI %.2f, by chain %s | R-hat table %.3f (ESS %.0f), alpha %.3f | div %d | %.1f min\n",
              i$name, J, R, C, row$ei, row$ei_chains, rh, es, rha, row$div, row$mins))
}
