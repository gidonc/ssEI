## Saving, loading and tabulating benchmark results. Results are RDS files in experiments/benchmark/results/ named
## <election>_<method>.rds, so runs made on different machines can be copied into one folder and compared together.

if (!exists("BM_DIR")) BM_DIR <- "experiments/benchmark"
## where results are written and read: $BM_RESULTS if set (e.g. a scratch directory on a cluster), else experiments/benchmark/results
bm_results_dir <- function() { d <- Sys.getenv("BM_RESULTS", ""); if (nzchar(d)) d else file.path(BM_DIR, "results") }

bm_save <- function(res, election) {
  dir.create(bm_results_dir(), recursive = TRUE, showWarnings = FALSE)
  saveRDS(res, file.path(bm_results_dir(), sprintf("%s_%s.rds", election, gsub("[^A-Za-z0-9]", "", res$method))))
}

## one line per method: scores, convergence over every cell, time
bm_report <- function(results, kc) {
  rm <- apply(kc, c(1, 2), sum); cm <- apply(kc, c(1, 3), sum)
  rows <- lapply(results, function(res) {
    s <- bm_score(res$est, kc, res$draws); cv <- res$conv; g <- bm_margin_gap(res$est, rm, cm)
    data.frame(method = res$method, ei_area = s$ei_area, ei_area_median = s$ei_area_median, ei_total = s$ei_total, wpe_total = s$wpe_total,
               cover_area = s$cover_area, cover_total = s$cover_total, gap_row = g[["gap_row"]], gap_col = g[["gap_col"]],
               max_rhat = if (is.null(cv)) NA else cv$max_rhat, min_ess_bulk = if (is.null(cv)) NA else cv$min_ess_bulk,
               converged = if (is.null(cv)) NA else cv$converged, mins = res$secs / 60)
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL
  out[sapply(out, is.numeric)] <- lapply(out[sapply(out, is.numeric)], round, 3)
  out
}


## all saved results for an election, as a named list
bm_load <- function(election) {
  f <- list.files(bm_results_dir(), pattern = sprintf("^%s_.*\\.rds$", election), full.names = TRUE)
  stopifnot(length(f) > 0)
  res <- lapply(f, readRDS); names(res) <- vapply(res, `[[`, "", "method"); res
}

## Save an ssEI fit (the stanfit that rw_fit returns) in the same form as the other methods, e.g. at the end of an experiment script:
##   if (isTRUE(BM_SAVE)) bm_ssei_save(fit, kc, "scotland7")
bm_ssei_save <- function(fit, kc, election, label = "ssEI") {
  res <- bm_ssei(fit, apply(kc, c(1, 2), sum), label = label); bm_save(res, election)
  print(bm_report(list(res), kc)); invisible(res)
}
