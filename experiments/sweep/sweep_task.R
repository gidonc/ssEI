## One task of the simulation sweep: simulate a dataset and fit one method. Run from the repository root:
##   Rscript experiments/sweep/sweep_task.R <task number> [output directory, default sweep_out] [grid file]
## Methods: refs (independence, Goodman raked, true pooled table raked) | king (eiPack ei.MD.bayes) |
##          ssei_est (sigma estimated) | ssei_s08 (sigma fixed 0.80) | ssei_kingqi (Dirichlet(4) mean prior on quasi-independence + sigma 0.80)
a <- commandArgs(TRUE); task <- as.integer(a[1]); out <- if (length(a) >= 2) a[2] else "sweep_out"; gridf <- if (length(a) >= 3) a[3] else "experiments/sweep/grid.csv"
g <- read.csv(gridf); p <- g[g$task == task, ]; stopifnot(nrow(p) == 1); dir.create(out, showWarnings = FALSE, recursive = TRUE)
tag <- sprintf("t%04d_a%.1f_s%.2f_m%.0f_e%.1f_g%.2f_%d_%s", task, p$assoc, p$sort, p$assim, p$seg, p$sigma, p$seed, p$method)
outf <- file.path(out, paste0(tag, ".rds")); if (file.exists(outf)) { cat("exists:", outf, "\n"); quit(save = "no") }
`%||%` <- function(x, y) if (is.null(x)) y else x; options(max.print = 5)
suppressMessages({ source("experiments/rowwise_helpers.R"); source("experiments/sim_tables.R") })
x <- sim_tables(J = 60, R = 4, C = 4, assoc = p$assoc, sort = p$sort, sigma = p$sigma, seg = p$seg, seed = p$seed, assim = p$assim)
kc <- x$kc; J <- dim(kc)[1]; R <- dim(kc)[2]; C <- dim(kc)[3]; dimnames(kc) <- list(NULL, paste0("r", 1:R), paste0("c", 1:C))
N <- sum(kc); kt <- apply(kc, 2:3, sum); rm <- apply(kc, 1:2, sum); cm <- apply(kc, c(1, 3), sum)
rake <- function(seed) { o <- array(0, dim(kc)); for (j in 1:J) { X <- seed / sum(seed) + 1e-12
  for (i in 1:500) { X <- X * (rm[j, ] / pmax(rowSums(X), 1e-300)); X <- sweep(X, 2, cm[j, ] / pmax(colSums(X), 1e-300), "*") }; o[j, , ] <- X }; o }
t0 <- Sys.time(); cover <- c(area = NA, national = NA); diag <- NULL
if (p$method == "refs") {
  library(quadprog); xs <- rm / rowSums(rm); ys <- cm / rowSums(cm); w <- rowSums(rm); D <- matrix(0, R * C, R * C); dv <- numeric(R * C)
  for (j in 1:J) { D <- D + w[j] * kronecker(tcrossprod(xs[j, ]), diag(C)); dv <- dv + w[j] * kronecker(xs[j, ], ys[j, ]) }
  D <- 2 * D / sum(w) + diag(1e-9, R * C); dv <- 2 * dv / sum(w); A <- cbind(t(kronecker(diag(R), matrix(1, 1, C))), diag(R * C))
  gb <- matrix(solve.QP(D, dv, A, c(rep(1, R), rep(0, R * C)), meq = R)$solution, R, C, byrow = TRUE); gb[gb < 0] <- 0; gb <- gb / rowSums(gb)
  ind <- array(0, dim(kc)); for (j in 1:J) ind[j, , ] <- outer(rm[j, ], cm[j, ]) / sum(rm[j, ])
  est <- list(indep = ind, goodman_raked = rake(gb * rowSums(kt)), rake_true = rake(kt))
} else if (p$method == "king") {
  library(eiPack); df <- data.frame(rm, cm); names(df) <- c(paste0("R", 1:R), paste0("C", 1:C))
  f <- as.formula(paste0("cbind(", paste0("C", 1:C, collapse = ","), ") ~ cbind(", paste0("R", 1:R, collapse = ","), ")"))
  set.seed(1); tl <- tuneMD(f, data = df, ntunes = 10, totaldraws = 50000)
  md <- ei.MD.bayes(f, data = df, sample = 1000, thin = 100, burnin = 200000, ret.beta = "r", verbose = 0, tune.list = tl)
  bd <- md$draws$Beta; e <- array(0, dim(kc)); lo <- e; hi <- e; dr <- array(0, c(nrow(bd), J, R, C))
  for (r in 1:R) for (c in 1:C) for (j in 1:J) dr[, j, r, c] <- bd[, paste0("beta.R", r, ".C", c, ".", j)] * rm[j, r]
  est <- list(king = apply(dr, 2:4, mean)); lo <- apply(dr, 2:4, quantile, 0.05); hi <- apply(dr, 2:4, quantile, 0.95); tot <- apply(dr, c(1, 3, 4), sum)
  cover <- c(area = mean((kc >= lo & kc <= hi)[kc >= 50]), national = mean(kt >= apply(tot, 2:3, quantile, 0.05) & kt <= apply(tot, 2:3, quantile, 0.95)))
} else {
  Vr <- rw_sbp_to_v(rw_cluster_sbp(rm + 0.5)); Vc <- rw_sbp_to_v(rw_cluster_sbp(cm + 0.5)); V <- cbind(kronecker(Vr, matrix(1 / sqrt(C), C, 1)), kronecker(diag(R), Vc))
  args <- list(rem_col = rep(which.max(colSums(cm)), R), tiers = rep(1, C - 1), loyal_rows = integer(0), E_centre = "qi")
  args <- switch(p$method, ssei_est = args, ssei_s08 = c(args, list(sigma_fixed = sqrt(trigamma(2)))),
                 ssei_kingqi = c(args, list(E_alpha = 4, sigma_fixed = sqrt(trigamma(2)))), stop("unknown method"))
  fit <- do.call(rw_fit, c(list(kc, V, alloc = 3L), args, list(chains = 2, cores = 2, iter = 800, warmup = 400, seed = 1234, refresh = 0)))
  cv <- rstan::extract(fit, "cell_values")$cell_values; est <- setNames(list(apply(cv, 2:4, mean)), p$method)
  lo <- apply(cv, 2:4, quantile, 0.05); hi <- apply(cv, 2:4, quantile, 0.95); tot <- apply(cv, c(1, 3, 4), sum)
  cover <- c(area = mean((kc >= lo & kc <= hi)[kc >= 50]), national = mean(kt >= apply(tot, 2:3, quantile, 0.05) & kt <= apply(tot, 2:3, quantile, 0.95)))
  sp <- do.call(rbind, rstan::get_sampler_params(fit, inc_warmup = FALSE)); s <- rstan::summary(fit, pars = "E_rc")$summary
  diag <- c(rhat = max(s[, "Rhat"]), n_eff_min = min(s[, "n_eff"]), divergent = sum(sp[, "divergent__"]))
}
res <- lapply(est, function(e) c(err = 50 * sum(abs(e - kc)) / N, nat = 50 * sum(abs(apply(e, 2:3, sum) - kt)) / N))
## screening measures available from the margins alone (for relating results to observable features)
xs <- rm / rowSums(rm); screen <- c(seg_obs = mean(apply(log(xs + 1e-6), 2, sd)), min_row_share = mean(apply(xs, 2, min)))
saveRDS(list(params = p, res = res, cover = cover, diag = diag, screen = screen, mins = as.numeric(difftime(Sys.time(), t0, units = "mins"))), outf)
for (m in names(res)) cat(sprintf("DONE %s %s: error %.2f national %.2f | cover %s\n", tag, m, res[[m]]["err"], res[[m]]["nat"], paste(round(cover, 2), collapse = "/")))
