## Scoring and convergence summaries shared by every method in the benchmark.
##
## Everything works on an areas x R x C array of cell counts, so any method that can produce one can be scored:
##   est    J x R x C point estimate (posterior mean for the Bayesian methods)
##   draws  S x J x R x C posterior draws, or NULL for methods without them
##   kc     J x R x C true cells
##
## Scores
##   ei_area         error index over all area cells: 50 * sum|est - kc| / N        (the same number rw_report prints)
##   ei_area_median  median over areas of each area's own error index
##   ei_total        error index of the table summed over the areas                  (Pavia and Romero's EI)
##   wpe_total       weighted proportion error of the summed table, within-row proportions weighted by true counts (their WPE)
##   cover_area      share of area cells with 50+ true voters inside their central 90% interval
##   cover_total     share of cells of the summed table inside the central 90% interval of the summed draws
##   width_total     mean width of those intervals as a share of the cell's true count (only cells with 50+ voters)

bm_score <- function(est, kc, draws = NULL, level = 0.90, big = 50) {
  stopifnot(identical(unname(dim(est)), unname(dim(kc))), all(is.finite(est)))
  N <- sum(kc); a <- (1 - level) / 2
  nj <- apply(kc, 1, sum)
  ei_j <- 50 * apply(abs(est - kc), 1, sum) / nj
  est_t <- apply(est, 2:3, sum); kc_t <- apply(kc, 2:3, sum)
  rowprop <- function(m) { rs <- rowSums(m); m / ifelse(rs > 0, rs, 1) }      # an empty row has proportions 0, not NaN
  p_t <- rowprop(kc_t); phat_t <- rowprop(est_t)
  out <- list(
    ei_area = 50 * sum(abs(est - kc)) / N,
    ei_area_median = median(ei_j[nj > 0]),
    ei_area_max = max(ei_j[nj > 0]),
    ei_total = 50 * sum(abs(est_t - kc_t)) / N,
    wpe_total = 100 * sum(kc_t * abs(p_t - phat_t)) / N,
    cover_area = NA_real_, cover_total = NA_real_, width_total = NA_real_,
    n_big_area = sum(kc >= big), n_big_total = sum(kc_t >= big))
  if (!is.null(draws)) {
    stopifnot(identical(unname(dim(draws)[-1]), unname(dim(kc))))
    lo <- apply(draws, 2:4, quantile, a, names = FALSE); hi <- apply(draws, 2:4, quantile, 1 - a, names = FALSE)
    out$cover_area <- mean((kc >= lo & kc <= hi)[kc >= big])
    tot <- apply(draws, c(1, 3, 4), sum)                       # draw by draw, so the interval is that of the total
    lo_t <- apply(tot, 2:3, quantile, a, names = FALSE); hi_t <- apply(tot, 2:3, quantile, 1 - a, names = FALSE)
    out$cover_total <- mean(kc_t >= lo_t & kc_t <= hi_t)
    out$width_total <- mean(((hi_t - lo_t) / kc_t)[kc_t >= big])
  }
  out$ei_by_area <- ei_j
  out
}

## Convergence over everything the method reports, not a few chosen parameters.
##   x    iterations x chains x variables array (names on the third dimension)
## Variables that never change (an empty row's cells are always 0) carry no information and are dropped.
bm_conv <- function(x, rhat_ok = 1.01, ess_ok = 400) {
  stopifnot(length(dim(x)) == 3)
  sdv <- apply(x, 3, sd); keep <- is.finite(sdv) & sdv > 0
  s <- posterior::summarise_draws(posterior::as_draws_array(x[, , keep, drop = FALSE]),
                                  "rhat", "ess_bulk", "ess_tail")
  list(n_chains = dim(x)[2], n_iter = dim(x)[1], n_vars = sum(keep), n_constant = sum(!keep),
       max_rhat = max(s$rhat, na.rm = TRUE), frac_rhat_gt_1.01 = mean(s$rhat > 1.01, na.rm = TRUE),
       frac_rhat_gt_1.05 = mean(s$rhat > 1.05, na.rm = TRUE),
       min_ess_bulk = min(s$ess_bulk, na.rm = TRUE), min_ess_tail = min(s$ess_tail, na.rm = TRUE),
       converged = isTRUE(max(s$rhat, na.rm = TRUE) <= rhat_ok && min(s$ess_bulk, na.rm = TRUE) >= ess_ok),
       worst = s$variable[order(-s$rhat)][seq_len(min(5, nrow(s)))])
}

## S x J x R x C draws with the chain of each draw  ->  iterations x chains x (J*R*C) array for bm_conv
bm_chain_array <- function(draws, chain) {
  S <- dim(draws)[1]; nc <- length(unique(chain)); stopifnot(S %% nc == 0, all(table(chain) == S / nc))
  d <- dim(draws)[-1]; flat <- matrix(draws, S, prod(d))
  idx <- arrayInd(seq_len(prod(d)), d)
  colnames(flat) <- sprintf("cell[%d,%d,%d]", idx[, 1], idx[, 2], idx[, 3])
  out <- array(NA_real_, c(S / nc, nc, ncol(flat)), dimnames = list(NULL, NULL, colnames(flat)))
  for (k in seq_len(nc)) out[, k, ] <- flat[chain == sort(unique(chain))[k], , drop = FALSE]
  out
}

## Draws of an ssEI fit in the same form: S x J x R x C counts plus the chain of each draw.
bm_from_stanfit <- function(fit, pars = "cell_values") {
  a <- rstan::extract(fit, pars, permuted = FALSE)             # iterations x chains x parameters
  dn <- dimnames(a)[[3]]
  ij <- do.call(rbind, lapply(regmatches(dn, gregexpr("[0-9]+", dn)), as.integer))
  nit <- dim(a)[1]; nch <- dim(a)[2]; S <- nit * nch
  D <- array(0, c(S, max(ij[, 1]), max(ij[, 2]), max(ij[, 3])))
  flat <- matrix(a, S, dim(a)[3])                                 # iteration fastest, then chain
  for (k in seq_len(ncol(flat))) D[, ij[k, 1], ij[k, 2], ij[k, 3]] <- flat[, k]
  list(draws = D, chain = rep(seq_len(nch), each = nit))
}

bm_summarise <- function(draws) apply(draws, 2:4, mean)

## How far the estimated tables are from reproducing the margins they were fitted to, as an error index on the row and on the column
## margins (50 * sum|margin of est - margin| / N). Zero for methods that fit the margins exactly. A large gap on one side only
## means rows and columns have been swapped somewhere.
bm_margin_gap <- function(est, rm, cm) {
  N <- sum(rm)
  c(gap_row = 50 * sum(abs(apply(est, c(1, 2), sum) - as.matrix(rm))) / N, gap_col = 50 * sum(abs(apply(est, c(1, 3), sum) - as.matrix(cm))) / N)
}
