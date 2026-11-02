## ssEI on the per-district benchmark with the generic set-up of the simulation sweep (clustered bisection basis, no loyalty):
##   est:    sigma estimated, QI prior centre
##   king:   King-like: Dirichlet(4) prior on the average table's row splits (E_alpha = 4, uniform centre 1/C, as eiPack's
##           Gamma(4, 2) gives), area spread fixed at sqrt(trigamma(2)) = 0.80 (an area Dirichlet at King's prior-mean alpha of about 2)
##   kingqi: as king with the prior centred on quasi-independence (the sweep's ssei_kingqi)
## usage: Rscript experiments/pavia_bench/ssei_kinglike.R <dataset> <n districts> <out.rds> [which of the sampled districts, e.g. 1,2] [variants, e.g. king,kingqi,est]
a <- commandArgs(TRUE); ds <- a[1]; nd <- as.integer(a[2]); out <- a[3]
suppressMessages({ source("experiments/rowwise_helpers.R"); source("experiments/pavia_bench/instances.R") })
ins <- pavia_instances(ds); set.seed(20261101); pick <- sort(sample(length(ins), nd))
if (length(a) >= 4) pick <- pick[as.integer(strsplit(a[4], ",")[[1]])]
vars <- if (length(a) >= 5) strsplit(a[5], ",")[[1]] else c("king", "kingqi", "est")
res <- list()
for (p in pick) {
  i <- ins[[p]]; rm <- i$rm; cm <- i$cm; J <- nrow(rm); R <- ncol(rm); C <- ncol(cm)
  kc <- array(0, c(J, R, C)); for (j in 1:J) kc[j, , ] <- outer(rm[j, ], cm[j, ]) / sum(rm[j, ])   # only the margins are used
  Vr <- rw_sbp_to_v(rw_cluster_sbp(rm + 0.5)); Vc <- rw_sbp_to_v(rw_cluster_sbp(cm + 0.5)); V <- cbind(kronecker(Vr, matrix(1 / sqrt(C), C, 1)), kronecker(diag(R), Vc))
  for (v in vars) {
    args <- list(rem_col = rep(which.max(colSums(cm)), R), tiers = rep(1, C - 1), loyal_rows = integer(0))
    args <- switch(v, est = c(args, list(E_centre = "qi")),
                   king = c(args, list(E_centre = "uniform", E_alpha = 4, sigma_fixed = sqrt(trigamma(2)))),
                   kingqi = c(args, list(E_centre = "qi", E_alpha = 4, sigma_fixed = sqrt(trigamma(2)))))
    t0 <- Sys.time()
    fit <- tryCatch(suppressMessages(do.call(rw_fit, c(list(kc, V, alloc = 3L), args, list(chains = 4, cores = 2, iter = 1000, warmup = 500, seed = 1234, refresh = 0)))),
                    error = function(e) { cat("  error:", conditionMessage(e), "\n"); NULL })
    if (is.null(fit)) next
    cv <- rstan::extract(fit, "cell_values")$cell_values; tab <- apply(cv, c(1, 3, 4), sum)
    perm <- rstan::extract(fit, "cell_values", permuted = FALSE); nch <- dim(perm)[2]
    tabc <- sapply(1:nch, function(k) { x <- array(colMeans(perm[, k, ]), c(J, R, C)); as.vector(apply(x, 2:3, sum)) })
    ei_ch <- apply(tabc, 2, function(t) ei_index(matrix(t, R, C), i$truth))
    s <- rstan::summary(fit, pars = "E_rc")$summary; sp <- do.call(rbind, rstan::get_sampler_params(fit, inc_warmup = FALSE))
    row <- data.frame(district = i$name, variant = v, units = J, R = R, C = C, ei = ei_index(apply(tab, 2:3, mean), i$truth),
                      ei_chains = paste(sprintf("%.2f", ei_ch), collapse = " "), rhat = max(s[, "Rhat"], na.rm = TRUE), div = sum(sp[, "divergent__"]),
                      mins = as.numeric(Sys.time() - t0, units = "mins"))
    res[[length(res) + 1]] <- row; saveRDS(do.call(rbind, res), out)
    cat(sprintf("%-22s %-7s units %3d %dx%d | EI %.2f (chains %s) | E R-hat %.3f | div %d | %.1f min\n", i$name, v, J, R, C, row$ei, row$ei_chains, row$rhat, row$div, row$mins))
  }
}
