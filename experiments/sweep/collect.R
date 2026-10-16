## Collect sweep results: Rscript experiments/sweep/collect.R [output directory, default sweep_out]
a <- commandArgs(TRUE); out <- if (length(a) >= 1) a[1] else "sweep_out"
fs <- list.files(out, pattern = "\\.rds$", full.names = TRUE); rows <- list()
for (f in fs) { x <- readRDS(f); p <- x$params
  for (m in names(x$res)) rows[[length(rows) + 1]] <- data.frame(p[, c("assoc", "sort", "assim", "seg", "sigma", "seed")], method = m,
    err = x$res[[m]]["err"], nat = x$res[[m]]["nat"], cover = if (m %in% c("indep", "goodman_raked", "rake_true")) NA else x$cover[["area"]],
    cover_nat = if (m %in% c("indep", "goodman_raked", "rake_true")) NA else x$cover[["national"]],
    rhat = if (is.null(x$diag)) NA else x$diag[["rhat"]], mins = x$mins, row.names = NULL) }
d <- do.call(rbind, rows); write.csv(d, file.path(out, "results_long.csv"), row.names = FALSE)
cat(nrow(d), "rows from", length(fs), "files; written to", file.path(out, "results_long.csv"), "\n")
## mean over seeds, wide by method
m <- aggregate(err ~ assoc + sort + assim + seg + method, d, mean)
w <- reshape(m, idvar = c("assoc", "sort", "assim", "seg"), timevar = "method", direction = "wide"); names(w) <- sub("^err\\.", "", names(w))
w <- w[order(w$assoc, w$sort, w$assim, w$seg), ]; options(width = 220); print(w, row.names = FALSE, digits = 3)
cat("\nmean coverage (area / national) and worst Rhat by method:\n")
print(aggregate(cbind(cover, cover_nat, rhat) ~ method, d[!is.na(d$cover), ], function(v) round(c(mean = mean(v)), 2), na.action = na.pass), row.names = FALSE)
