S <- Sys.getenv("S"); fs <- list.files(file.path(S, "sim"), pattern = "^sim_.*rds$", full.names = TRUE)
rows <- list()
for (f in fs) { x <- readRDS(f); st <- x$settings
  for (m in names(x$res)) rows[[length(rows) + 1]] <- data.frame(assoc = st$assoc, sort = st$sort, sigma = st$sigma, seg = st$seg, seed = st$seed, method = m,
    err = x$res[[m]]["err"], nat = x$res[[m]]["nat"], cover = if (m %in% c("indep", "goodman_raked", "rake_true")) NA else x$cover[1],
    cover_nat = if (m %in% c("indep", "goodman_raked", "rake_true")) NA else x$cover[2], rhat = if (is.null(x$diag)) NA else x$diag["rhat"]) }
d <- do.call(rbind, rows); rownames(d) <- NULL
ord <- c("rake_true", "indep", "goodman_raked", "ssei_est", "ssei_s08", "ssei_kingqi", "king"); d$method <- factor(d$method, levels = intersect(ord, unique(d$method)))
w <- reshape(d[, c("assoc", "sort", "sigma", "seg", "method", "err")], idvar = c("assoc", "sort", "sigma", "seg"), timevar = "method", direction = "wide")
names(w) <- sub("^err\\.", "", names(w)); w <- w[order(w$seg, w$assoc, w$sort, w$sigma), c("seg", "assoc", "sort", "sigma", levels(d$method))]
options(width = 200); cat("Area error by cell\n"); print(w, row.names = FALSE, digits = 3)
wn <- reshape(d[, c("assoc", "sort", "sigma", "seg", "method", "nat")], idvar = c("assoc", "sort", "sigma", "seg"), timevar = "method", direction = "wide")
names(wn) <- sub("^nat\\.", "", names(wn)); wn <- wn[order(wn$seg, wn$assoc, wn$sort, wn$sigma), c("seg", "assoc", "sort", "sigma", levels(d$method))]
cat("\nNational-table error by cell\n"); print(wn, row.names = FALSE, digits = 3)
cv <- aggregate(cbind(cover, cover_nat, rhat) ~ method, d[!is.na(d$cover), ], function(v) round(mean(v), 2), na.action = na.pass); cat("\nMean coverage (area cells 50+, national) and max Rhat by method\n"); print(cv, row.names = FALSE)
saveRDS(d, file.path(S, "sim_pilot_results.rds"))
