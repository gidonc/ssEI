## usage: Rscript run_sim.R <assoc> <sort> <sigma> <seed> <method: refs | king | ssei_est | ssei_s08 | ssei_s04 | ssei_kingqi>
a <- commandArgs(TRUE); S <- Sys.getenv("S"); assoc <- as.numeric(a[1]); srt <- as.numeric(a[2]); sg <- as.numeric(a[3]); seed <- as.integer(a[4]); method <- a[5]; seg <- if (length(a) >= 6) as.numeric(a[6]) else 0.8
options(max.print = 5); `%||%` <- function(x, y) if (is.null(x)) y else x
suppressMessages({ source("experiments/rowwise_helpers.R"); source("experiments/sim_tables.R") })
x <- sim_tables(J = 60, R = 4, C = 4, assoc = assoc, sort = srt, sigma = sg, seg = seg, seed = seed); kc <- x$kc; J <- dim(kc)[1]; R <- dim(kc)[2]; C <- dim(kc)[3]
dimnames(kc) <- list(NULL, paste0("r", 1:R), paste0("c", 1:C)); N <- sum(kc); kt <- apply(kc, 2:3, sum); rm <- apply(kc, 1:2, sum); cm <- apply(kc, c(1, 3), sum)
tag <- sprintf("sim_a%.1f_s%.1f_g%.2f_e%.1f_%d_%s", assoc, srt, sg, seg, seed, method); t0 <- Sys.time(); cover <- c(NA, NA)
if (method == "refs") {
  library(quadprog); sc <- readLines(file.path(S, "screen.R")); eval(parse(text = sc[grep("^rake <-", sc):(grep("^dll <-", sc) - 1)]))
  ind <- array(0, dim(kc)); for (j in 1:J) ind[j, , ] <- outer(rm[j, ], cm[j, ]) / sum(rm[j, ])
  gb <- goodman(kc) * rowSums(kt)
  out <- list(indep = ind, goodman_raked = rake(kc, gb), rake_true = rake(kc, kt))
} else if (method == "king") {
  library(eiPack); df <- data.frame(rm, cm); names(df) <- c(paste0("R", 1:R), paste0("C", 1:C))
  f <- as.formula(paste0("cbind(", paste0("C", 1:C, collapse = ","), ") ~ cbind(", paste0("R", 1:R, collapse = ","), ")"))
  set.seed(1); tl <- tuneMD(f, data = df, ntunes = 10, totaldraws = 50000)
  md <- ei.MD.bayes(f, data = df, sample = 1000, thin = 100, burnin = 200000, ret.beta = "r", verbose = 0, tune.list = tl)
  bd <- md$draws$Beta; est <- array(0, dim(kc)); lo <- est; hi <- est
  for (r in 1:R) for (c in 1:C) for (j in 1:J) { v <- bd[, paste0("beta.R", r, ".C", c, ".", j)] * rm[j, r]; est[j, r, c] <- mean(v); lo[j, r, c] <- quantile(v, 0.05); hi[j, r, c] <- quantile(v, 0.95) }
  out <- list(king = est); cover <- c(mean((kc >= lo & kc <= hi)[kc >= 50]), NA)
} else {
  Vr <- rw_sbp_to_v(rw_cluster_sbp(rm + 0.5)); Vc <- rw_sbp_to_v(rw_cluster_sbp(cm + 0.5)); V <- cbind(kronecker(Vr, matrix(1 / sqrt(C), C, 1)), kronecker(diag(R), Vc))
  args <- list(rem_col = rep(which.max(colSums(cm)), R), tiers = rep(1, C - 1), loyal_rows = integer(0), E_centre = "qi")
  args <- switch(method, ssei_est = args, ssei_s08 = c(args, list(sigma_fixed = sqrt(trigamma(2)))), ssei_s04 = c(args, list(sigma_fixed = 0.4)),
                 ssei_kingqi = c(args[-4], list(E_alpha = 4, E_centre = "qi", sigma_fixed = sqrt(trigamma(2)))))
  fit <- do.call(rw_fit, c(list(kc, V, alloc = 3L), args, list(chains = 2, cores = 2, iter = 800, warmup = 400, seed = 1234, refresh = 0)))
  cv <- rstan::extract(fit, "cell_values")$cell_values; est <- apply(cv, 2:4, mean); lo <- apply(cv, 2:4, quantile, 0.05); hi <- apply(cv, 2:4, quantile, 0.95)
  tot <- apply(cv, c(1, 3, 4), sum); out <- setNames(list(est), method)
  cover <- c(mean((kc >= lo & kc <= hi)[kc >= 50]), mean(kt >= apply(tot, 2:3, quantile, 0.05) & kt <= apply(tot, 2:3, quantile, 0.95)))
  sp <- do.call(rbind, rstan::get_sampler_params(fit, inc_warmup = FALSE)); s <- rstan::summary(fit, pars = "E_rc")$summary
  attr(out, "diag") <- c(rhat = max(s[, "Rhat"]), div = sum(sp[, "divergent__"]))
}
res <- lapply(out, function(e) c(err = 50 * sum(abs(e - kc)) / N, nat = 50 * sum(abs(apply(e, 2:3, sum) - kt)) / N))
saveRDS(list(settings = x$settings, method = method, res = res, cover = cover, diag = attr(out, "diag"), mins = as.numeric(difftime(Sys.time(), t0, units = "mins"))), file.path(S, "sim", paste0(tag, ".rds")))
for (m in names(res)) cat(sprintf("DONE %s %s: error %.2f national %.2f | cover %s | %.1f min\n", tag, m, res[[m]]["err"], res[[m]]["nat"], paste(round(cover, 2), collapse = "/"), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
