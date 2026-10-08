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


## ---------------------------------------------------------------------------
## Into the package's own evaluation, mods_summary() in R/eval_ei.R
## ---------------------------------------------------------------------------
## A saved benchmark result (from any method) as the two tables mods_summary() wants for each model: cell values (cv) and row rates (rr), with
## columns area_no, row_no, col_no, mean, sd and the 2.5/25/50/75/97.5% quantiles when the method has draws (mean only otherwise). Area, row and
## column numbers are positions in the data, as in kc. Row rates of an empty row are -1, as in ssEI.
bm_as_mod <- function(res, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
  est <- res$est; d <- dim(est); J <- d[1]; R <- d[2]; C <- d[3]
  grid <- expand.grid(area_no = seq_len(J), row_no = seq_len(R), col_no = seq_len(C))        # same order as a J x R x C array
  rmh <- apply(est, c(1, 2), sum)                                                             # row totals (every method reproduces them)
  rmv <- rmh[cbind(grid$area_no, grid$row_no)]
  stats <- function(D) {                                                                      # D: draws x cells
    q <- t(apply(D, 2, stats::quantile, probs = probs, names = FALSE)); colnames(q) <- paste0(probs * 100, "%")
    cbind(mean = colMeans(D), sd = apply(D, 2, stats::sd), q)
  }
  if (is.null(res$draws)) {
    cvs <- cbind(mean = as.numeric(est)); rrs <- cbind(mean = as.numeric(est) / rmv)
  } else {
    D <- res$draws; dim(D) <- c(dim(D)[1], J * R * C)
    cvs <- stats(D); rrs <- stats(sweep(D, 2, ifelse(rmv > 0, rmv, 1), "/"))
  }
  rrs[rmv == 0, ] <- -1
  mk <- function(m) tibble::as_tibble(cbind(grid, as.data.frame(m, check.names = FALSE)))
  structure(list(cv = mk(cvs), rr = mk(rrs), method = res$method), class = "bm_mod")
}

## the truth in the long form mods_summary() wants: area_no, row_no, col_no, actual_cell_value, actual_row_rate, actual_row_margin
bm_actual_long <- function(kc) {
  d <- dim(kc); grid <- expand.grid(area_no = seq_len(d[1]), row_no = seq_len(d[2]), col_no = seq_len(d[3]))
  rm <- apply(kc, c(1, 2), sum)
  out <- cbind(grid, actual_cell_value = as.numeric(kc), actual_row_margin = rm[cbind(grid$area_no, grid$row_no)])
  out$actual_row_rate <- ifelse(out$actual_row_margin > 0, out$actual_cell_value / out$actual_row_margin, NA_real_)
  tibble::as_tibble(out)
}

## THE SIMPLE ROUTE from saved results: results = a list of saved results (bm_load("scotland7"), or list(readRDS(file))), kc = true cells.
##   res <- readRDS("experiments/benchmark/results/scotland7_eiMDbayes.rds")
##   s <- bm_mods_summary(list(res), D$kc)      # D <- bm_data("scotland7"); s$cv_eval$cv_eval_stats, s$eval_plots$p_cv, ...
## Models can be mixed with fits mods_summary() already knows (stanfit, eiMD, ...) by calling mods_summary() yourself with bm_as_mod(res).
bm_mods_summary <- function(results, kc) {
  if (!exists("mods_summary")) {
    suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })
    source(file.path(if (exists("BM_DIR")) dirname(dirname(BM_DIR)) else ".", "R", "eval_ei.R"))
  }
  mods <- lapply(results, bm_as_mod); names(mods) <- vapply(results, function(r) r$method, "")
  mods_summary(mods, bm_actual_long(kc))
}
