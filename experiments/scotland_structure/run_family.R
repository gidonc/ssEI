## Scotland 7x7, 73 areas, average table restricted to a family (E_family): the default row-wise set-up (as the 'sco_base' run: loyalty for the four major parties, two sigma tiers,
## off-diagonal prior scale 0.5, free interior, 2 chains x 1000), with the average table restricted by E_family.
## usage: Rscript experiments/scotland_structure/run_family.R <tag> <out dir> [family, default qi] [prior sd of the free coordinates, default 1.5] [sigma fixed | none] [interior free | fixed | scaled] [held, e.g. spoilt:0.10 | none]
a <- commandArgs(TRUE); tag <- a[1]; out <- a[2]; fam <- if (length(a) >= 3) a[3] else "qi"; fsd <- if (length(a) >= 4) as.numeric(a[4]) else 1.5; offs <- "none"; osd <- NULL
sfix <- if (length(a) >= 5 && a[5] != "none") as.numeric(a[5]) else NULL   # 5th argument: sigma fixed (weaker or stronger pooling), default estimated
intr <- if (length(a) >= 6) a[6] else "free"                                   # 6th argument: interior free | fixed | scaled
## 7th argument: held directions, e.g. spoilt:0.10 (the row's own-column share held, not learned); several separated by commas
hold <- if (length(a) >= 7 && a[7] != "none") lapply(strsplit(strsplit(a[7], ",")[[1]], ":"), function(x) list(row = x[1], share = as.numeric(x[2]))) else NULL
source("experiments/rowwise_helpers.R"); suppressMessages(source("experiments/data_prep.R"))
d <- prep_scotland(7, 73); loyal <- d$row_names[1:4]; if (!is.null(osd)) loyal <- union(loyal, strsplit(offs, ",")[[1]]); R <- 7; C <- 7; kc <- d$kc; dimnames(kc) <- list(NULL, d$row_names, d$col_names)
V <- rw_row_priority_basis(d$row_names, d$col_names, NULL, bisect_order = colSums(d$cm), loyal_rows = loyal)
args <- list(rem_col = seq_len(R), tiers = c(1, 1, rep(2, C - 3)), loyalty_mean = qlogis(0.75), loyal_rows = loyal, E_sd_scale_offdiag = 0.5, E_centre = "qi")
t0 <- Sys.time()
fit <- do.call(rw_fit, c(list(kc, V), args, list(E_family = if (fam != "none") fam, E_family_sd = fsd, E_family_fix = hold, sigma_fixed = sfix,
                                                            interior = intr, sigma_prior = if (intr == "free") c(log(0.3), 0.5) else c(log(0.7), 0.7), alloc = if (intr == "fixed") 4L else 3L, chains = 2, cores = 2, iter = 1000, warmup = 500, seed = 1234)))
mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cv <- rstan::extract(fit, "cell_values")$cell_values; N <- sum(kc); kt <- apply(kc, 2:3, sum)
est <- apply(cv, 2:4, mean); lo <- apply(cv, 2:4, quantile, 0.05); hi <- apply(cv, 2:4, quantile, 0.95); tot <- apply(cv, c(1, 3, 4), sum)
s <- rstan::summary(fit, pars = c("E_rc", "sigma_group_mu"))$summary; sp <- do.call(rbind, rstan::get_sampler_params(fit, inc_warmup = FALSE))
res <- list(tag = tag, offsets = offs, mins = mins, est = est, ei = 50 * sum(abs(est - kc)) / N, ei_sum = 50 * sum(abs(apply(est, 2:3, sum) - kt)) / N,
            ei_row = 50 * apply(abs(est - kc), 2, sum) / N, row_err = 50 * apply(abs(est - kc), 2, sum) / apply(kc, 2, sum),
            cover = mean((kc >= lo & kc <= hi)[kc >= 50]), cover_sum = mean(kt >= apply(tot, 2:3, quantile, 0.05) & kt <= apply(tot, 2:3, quantile, 0.95)),
            nat_rates = apply(est, 2:3, sum) / rowSums(apply(est, 2:3, sum)), div = sum(sp[, "divergent__"]), rhat = max(s[, "Rhat"]), neff = min(s[, "n_eff"]))
dir.create(out, showWarnings = FALSE, recursive = TRUE); saveRDS(res, file.path(out, paste0(tag, ".rds")))
cat(sprintf("DONE %s (family %s): area %.2f | national %.2f | row errors %s | cover %.2f / %.2f | div %d Rhat %.3f n_eff %.0f | %.1f min\n", tag, fam, res$ei, res$ei_sum,
            paste(sprintf("%.1f", res$row_err), collapse = " "), res$cover, res$cover_sum, res$div, res$rhat, res$neff, mins))
cat("spoilt row, national split (%): estimate", paste(sprintf("%.1f", 100 * res$nat_rates[5, ]), collapse = " "), "| truth", paste(sprintf("%.1f", 100 * kt[5, ] / sum(kt[5, ])), collapse = " "), "\n")
