## ssEI on the per-district benchmark instances (same sampled districts as king_conv.R).
##   generic: no loyalty structure, uniform prior centre (nothing about which candidate belongs to which party)
##   aligned: each party row paired with the candidate column whose polling-place shares correlate best with it (mutual best match,
##            correlation > 0.5; margins only), then the Scotland default set-up: loyalty for paired rows, QI centre, off-diagonal scale 0.5
## usage: Rscript experiments/pavia_bench/ssei_bench.R <dataset> <n districts> <out.rds> [generic|aligned|both]
a <- commandArgs(TRUE); ds <- a[1]; nd <- as.integer(a[2]); out <- a[3]; which_v <- if (length(a) >= 4) a[4] else "both"
suppressMessages({ source("experiments/rowwise_helpers.R"); source("experiments/pavia_bench/instances.R") })
ins <- pavia_instances(ds); set.seed(20261101); pick <- sort(sample(length(ins), nd))
align <- function(rm, cm) {   # reorder rows and columns so paired party r and candidate sit at the same index, pairs first
  xs <- rm / rowSums(rm); ys <- cm / rowSums(cm); K <- suppressWarnings(cor(xs, ys)); K[is.na(K)] <- 0
  pr <- integer(0); pc <- integer(0)
  for (r in seq_len(nrow(K))) { c <- which.max(K[r, ]); if (K[r, c] > 0.5 && which.max(K[, c]) == r) { pr <- c(pr, r); pc <- c(pc, c) } }
  ro <- c(pr, setdiff(seq_len(ncol(rm)), pr)); co <- c(pc, setdiff(seq_len(ncol(cm)), pc))
  list(ro = ro, co = co, n = length(pr))
}
res <- list()
for (p in pick) {
  i <- ins[[p]]; R <- ncol(i$rm); C <- ncol(i$cm); J <- nrow(i$rm)
  for (v in if (which_v == "both") c("generic", "aligned") else which_v) {
    rm <- i$rm; cm <- i$cm; ro <- seq_len(R); co <- seq_len(C); np <- 0
    if (v == "aligned") { al <- align(rm, cm); ro <- al$ro; co <- al$co; np <- al$n; rm <- rm[, ro, drop = FALSE]; cm <- cm[, co, drop = FALSE] }
    rn <- paste0("r", seq_len(R)); cn <- paste0("c", seq_len(C))
    kc <- array(0, c(J, R, C), dimnames = list(NULL, rn, cn))       # only the margins are used; cells set to independence
    for (j in 1:J) kc[j, , ] <- outer(rm[j, ], cm[j, ]) / sum(rm[j, ])
    loyal <- if (v == "aligned" && np > 0) rn[seq_len(min(np, C))] else integer(0)
    V <- rw_row_priority_basis(rn, cn, NULL, bisect_order = colSums(cm), loyal_rows = loyal)
    args <- list(kc, V, rem_col = pmin(seq_len(R), C), loyal_rows = loyal, chains = 2, cores = 2, iter = 1000, warmup = 500, seed = 1234, refresh = 0)
    if (v == "aligned" && np > 0) args <- c(args, list(loyalty_mean = qlogis(0.75), E_sd_scale_offdiag = 0.5, E_centre = "qi"))
    t0 <- Sys.time()
    fit <- tryCatch(suppressMessages(do.call(rw_fit, args)), error = function(e) { cat("  error:", conditionMessage(e), "\n"); NULL })
    if (is.null(fit)) next
    cv <- rstan::extract(fit, "cell_values")$cell_values; tab <- apply(cv, c(1, 3, 4), sum); est <- apply(tab, 2:3, mean)
    est0 <- matrix(0, R, C); est0[ro, co] <- est                     # back to the original order
    s <- rstan::summary(fit, pars = "E_rc")$summary
    row <- data.frame(district = i$name, variant = v, units = J, R = R, C = C, pairs = np, ei = ei_index(est0, i$truth),
                      rhat = max(s[, "Rhat"], na.rm = TRUE), mins = as.numeric(Sys.time() - t0, units = "mins"))
    res[[length(res) + 1]] <- row; saveRDS(do.call(rbind, res), out)
    cat(sprintf("%-28s %-8s units %3d %dx%d pairs %d | EI %.2f | Rhat %.3f | %.1f min\n", i$name, v, J, R, C, np, row$ei, row$rhat, row$mins))
  }
}
