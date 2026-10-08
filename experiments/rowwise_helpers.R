## Helpers for the row-by-row margin model, shared by the senc, Scotland and New Zealand examples.
##
## They build the extra model inputs (row order, reference columns, whitening matrices, anchor matrices, starting values), run the fit
## and report it. They are here so the examples run against the package as it stands; they are meant to move into the package.
##
## Allocation functions (Beyond RxC names -> lflag_seq_or_expected): adjusted row = 3 (default), odds-ratio logit = 2.

library(ssEI)

## ---------------------------------------------------------------------------
## 1. Helpers for the row-by-row model
## ---------------------------------------------------------------------------

## within-row basis from a sign matrix (one column per split: +1 / -1 / 0)
rw_sbp_to_v <- function(S) {
  V <- matrix(0, nrow(S), ncol(S))
  for (k in seq_len(ncol(S))) {
    p <- which(S[, k] == 1); n <- which(S[, k] == -1)
    nc <- sqrt(length(p) * length(n) / (length(p) + length(n)))
    V[p, k] <- nc / length(p); V[n, k] <- -nc / length(n)
  }
  V
}

## Goodman regression pattern (R x C, rows sum to 1): reference table and starting point only
rw_goodman <- function(rm, cm, floor = 0.02) {
  X <- as.matrix(rm) / rowSums(rm); Y <- as.matrix(cm) / rowSums(cm)
  B <- solve(crossprod(X), crossprod(X, Y))
  B[] <- pmax(B, floor)
  B / rowSums(B)
}

## row-by-row allocation, as in the Stan function seq_alloc_ordered (mode 3).
##   w, m      row and column shares (same total)
##   lam       (R - 1) * (C - 1) logits, laid out row by row in allocation order
##   row_order allocation order of the rows (last = remainder row)
##   rem_col   for each row, its reference column
rw_alloc <- function(w, m, lam, row_order, rem_col) {
  R <- length(w); C <- length(m); sc <- m; X <- matrix(0, R, C)
  for (k in seq_len(R - 1)) {
    r <- row_order[k]; eta <- numeric(C)
    eta[-rem_col[r]] <- lam[(k - 1) * (C - 1) + seq_len(C - 1)]
    lo <- -80; hi <- 80                                  # shift so that the row takes exactly w[r]
    for (i in 1:120) { t <- (lo + hi) / 2; if (sum(sc * plogis(eta + t)) > w[r]) hi <- t else lo <- t }
    X[r, ] <- sc * plogis(eta + (lo + hi) / 2); sc <- sc - X[r, ]
  }
  X[row_order[R], ] <- sc
  X
}

## the logits that reproduce table X under rw_alloc
rw_inv <- function(X, row_order, rem_col) {
  C <- ncol(X); sc <- colSums(X); lam <- numeric(0)
  for (r in row_order[-length(row_order)]) {
    e <- log(X[r, ] / (sc - X[r, ]))
    lam <- c(lam, (e - e[rem_col[r]])[-rem_col[r]]); sc <- sc - X[r, ]
  }
  lam
}

## cell-by-cell allocation with odds-ratio placement, as in the Stan function seq_alloc_ordered (mode 2).
## Same layout of lam, row order and reference columns as rw_alloc; each cell is placed by the log odds ratio of its 2 x 2 table.
rw_or_cell <- function(sr, sc, rest, lam) {
  psi <- exp(lam); A <- psi - 1; B <- -(psi * (sr + sc) + rest - sr); Cq <- psi * sr * sc
  disc <- sqrt(max(B * B - 4 * A * Cq, 0))
  if (B <= 0) 2 * Cq / (-B + disc) else (-B - disc) / (2 * A)
}
rw_alloc_or <- function(w, m, lam, row_order, rem_col) {
  R <- length(w); C <- length(m); sr <- w; sc <- m; X <- matrix(0, R, C)
  for (k in seq_len(R - 1)) {
    r <- row_order[k]; last <- rem_col[r]; rest <- sum(sc); i <- 0
    for (c in seq_len(C)) if (c != last) {
      i <- i + 1; rest <- rest - sc[c]
      x <- rw_or_cell(sr[r], sc[c], rest, lam[(k - 1) * (C - 1) + i])
      X[r, c] <- x; sc[c] <- sc[c] - x; sr[r] <- sr[r] - x
    }
    X[r, last] <- sr[r]; sc[last] <- sc[last] - sr[r]
  }
  X[row_order[R], ] <- sc
  X
}
rw_inv_or <- function(X, row_order, rem_col) {
  C <- ncol(X); sc <- colSums(X); sr <- rowSums(X); lam <- numeric(0)
  for (r in row_order[-length(row_order)]) {
    last <- rem_col[r]; rest <- sum(sc)
    for (c in seq_len(C)) if (c != last) {
      x <- X[r, c]; rest <- rest - sc[c]
      lam <- c(lam, log(x) + log(rest - sr[r] + x) - log(sr[r] - x) - log(sc[c] - x))
      sc[c] <- sc[c] - x; sr[r] <- sr[r] - x
    }
    sc[last] <- sc[last] - sr[r]
  }
  lam
}

## Inputs that switch on the row-by-row model with rows conditioned on the data.
##   mr        list from build_margin_reduction()
##   rm, cm    row and column margins (areas x R, areas x C)
##   row_order allocation order; default smallest row first, largest row as remainder
##   rem_col   reference column of each row
##   alloc     3 = row by row, one adjustment per row ("adjusted row"); 2 = cell by cell, odds-ratio placement ("odds-ratio logit")
rw_setup <- function(mr, rm, cm, row_order = order(colSums(rm)), rem_col = rep(1L, ncol(rm)),
                     q_bar = rw_goodman(rm, cm), alloc = 3L) {
  stopifnot(alloc %in% c(2L, 3L))
  alloc_fun <- if (alloc == 3L) rw_alloc else rw_alloc_or
  inv_fun   <- if (alloc == 3L) rw_inv else rw_inv_or
  rm <- as.matrix(rm); cm <- as.matrix(cm)
  J <- nrow(rm); R <- ncol(rm); C <- ncol(cm); Dm1 <- mr$Dm1_model; K <- (R - 1) * (C - 1)
  V <- mr$V_ilr_model
  stopifnot(Dm1 == R * (C - 1))                           # rows conditioned: within-row coordinates only
  W <- (rm + 0.01) / (rowSums(rm) + R * 0.01)             # the model's smoothed shares (empty rows become very small rows)
  M <- (cm + 0.01) / (rowSums(cm) + C * 0.01)
  b_of <- function(X, w) as.vector(crossprod(V, as.vector(t(log(X / w)))))
  ## reference matrix: d (hierarchy coordinates) / d (interior logits) at the Goodman pattern raked to each area's margins
  mp_B <- array(0, c(J, Dm1, K)); h <- 1e-4
  for (j in seq_len(J)) {
    w <- W[j, ]; m <- M[j, ]; X <- w * q_bar
    for (i in 1:1000) { X <- sweep(X, 2, m / colSums(X), "*"); X <- X * (w / rowSums(X)) }
    l0 <- inv_fun(X, row_order, rem_col)
    for (k in seq_len(K)) {
      e <- numeric(K); e[k] <- h
      mp_B[j, , k] <- (b_of(alloc_fun(w, m, l0 + e, row_order, rem_col), w) -
                       b_of(alloc_fun(w, m, l0 - e, row_order, rem_col), w)) / (2 * h)
    }
  }
  ## whitening of the column-margin parameters: factor of the covariance of the last-column log-ratios under the counts
  mp_T <- array(0, c(J, C - 1, C - 1))
  for (j in seq_len(J)) { n <- cm[j, ] + 0.01; mp_T[j, , ] <- t(chol(diag(1 / n[-C], C - 1) + 1 / n[C])) }
  modifyList(mr, list(
    lflag_margin_param = 1L, lflag_mp_seq = 1L, lflag_mp_seq_anchor = 1L, n_pin = C - 1L,
    mp_b_ref = matrix(0, J, Dm1), mp_G = array(0, c(J, C - 1, Dm1)),      # not used by the sequential model
    mp_P = array(0, c(J, Dm1, C - 1)), mp_N = array(0, c(J, Dm1, K)),
    lflag_seq_or_expected = alloc, mp_row_order = as.array(as.integer(row_order)), mp_rem_col = as.array(as.integer(rem_col)),
    mp_newton_iters = 10L,
    lflag_mp_beta_whiten = 1L, mp_T = mp_T,
    lflag_mp_vol_scale = 1L,
    lflag_mp_seq_scale = 2L, mp_B = mp_B, mp_nc_w = rep(1, Dm1), mp_sigma0 = rep(0.5, Dm1),
    lflag_mp_scale_fast = 0L
  ))
}

## starting values: every area at E_rc's table, volume and margin parameters at 0
rw_init <- function(J, n_pin, K, E_mean, n_sigma) {
  function() list(log_volume_raw = 0, log_volume_rest = as.array(rep(0, J - 1)),
                  mp_beta = matrix(0, J, n_pin), mp_z = matrix(0, J, K),
                  E_rc_group_mu = as.array(E_mean), sigma_group_mu = as.array(rep(log(0.3), n_sigma)))
}


## ---------------------------------------------------------------------------
## Basis, priors and the fit itself
## ---------------------------------------------------------------------------

## Row-priority basis for a table whose rows and columns are the same categories (parties), orthonormal, rows first:
## row-margin coordinates (Helmert on the log row shares), then within each row: loyalty vs the rest, congruent (same bloc) vs cross,
## then bisections of what is left. `bloc` is a named vector of bloc labels for the categories (NULL: no blocs, so the loyalty split is
## followed by bisections of all the other columns). `bisect_order` (e.g. column totals) puts the largest columns first within a bisection.
rw_bisect <- function(cells, add_split) {
  if (length(cells) <= 1) return(invisible(NULL))
  mid <- ceiling(length(cells) / 2)
  add_split(cells[1:mid], cells[(mid + 1):length(cells)])
  rw_bisect(cells[1:mid], add_split); rw_bisect(cells[(mid + 1):length(cells)], add_split)
  invisible(NULL)
}
rw_row_priority_basis <- function(row_names, col_names, bloc = NULL, bisect_order = NULL) {
  R <- length(row_names); C <- length(col_names)
  if (is.null(bloc)) bloc <- setNames(seq_along(col_names), col_names)     # no blocs: every category its own bloc
  srt <- function(ix) if (is.null(bisect_order)) ix else ix[order(-bisect_order[ix], ix)]   # largest first, so bisections group like with like
  helmert <- function(N) { V <- matrix(0, N, N - 1); for (k in 1:(N - 1)) { nc <- sqrt(k * (k + 1)); V[1:k, k] <- 1 / nc; V[k + 1, k] <- -k / nc }; V }
  V_within <- matrix(0, R * C, R * (C - 1))
  for (r in seq_len(R)) {
    S <- matrix(0, C, C - 1); col <- 1
    add_split <- function(pos, neg) { S[pos, col] <<- 1; S[neg, col] <<- -1; col <<- col + 1 }
    loyal <- match(row_names[r], col_names)
    congruent <- srt(setdiff(which(bloc[col_names] == bloc[row_names[r]]), loyal))
    cross <- srt(setdiff(seq_len(C), c(loyal, congruent)))
    if (length(congruent) + length(cross) > 0) add_split(loyal, c(congruent, cross))
    if (length(congruent) > 0) { if (length(cross) > 0) add_split(congruent, cross); rw_bisect(congruent, add_split) }
    rw_bisect(cross, add_split)
    stopifnot(col - 1 == C - 1)
    V_within[(r - 1) * C + 1:C, (r - 1) * (C - 1) + 1:(C - 1)] <- rw_sbp_to_v(S)
  }
  V <- cbind(kronecker(helmert(R), matrix(1 / sqrt(C), C, 1)), V_within)
  stopifnot(all.equal(crossprod(V), diag(R * C - 1), tolerance = 1e-8))
  V
}

## Centre for the non-loyalty coordinates: quasi-independence. Loyalty is estimated from the margins alone (Goodman regression table raked to
## each area's margins, as in rw_setup); the voters who do not stay are then spread over the other columns in proportion to what each
## column has left (IPF on the off-diagonal). Returns, for each coordinate, log(sum expected on + side) - log(sum on - side), NA if undefined.
rw_qi_centre <- function(V_model, rm, cm) {
  rm <- as.matrix(rm); cm <- as.matrix(cm); J <- nrow(rm); R <- ncol(rm); C <- ncol(cm)
  q_bar <- rw_goodman(rm, cm)
  W <- (rm + 0.01) / (rowSums(rm) + R * 0.01); M <- (cm + 0.01) / (rowSums(cm) + C * 0.01)
  T <- matrix(0, R, C)
  for (j in seq_len(J)) {
    X <- W[j, ] * q_bar
    for (i in 1:200) { X <- sweep(X, 2, M[j, ] / colSums(X), "*"); X <- X * (W[j, ] / rowSums(X)) }
    T <- T + sum(rm[j, ]) * X
  }
  off <- T; diag(off) <- 0; rs <- rowSums(off); cs <- colSums(off)
  E <- outer(rs, cs) * (1 - diag(R))
  for (i in 1:500) { E <- E * (rs / rowSums(E)); E[is.nan(E)] <- 0; E <- sweep(E, 2, cs / colSums(E), "*"); E[is.nan(E)] <- 0 }
  Ev <- as.vector(t(E))
  apply(V_model, 2, function(v) {
    sp <- sum(Ev[v > 1e-10]); sn <- sum(Ev[v < -1e-10])
    if (sp > 0 && sn > 0) log(sp) - log(sn) else NA_real_
  })
}

## Coordinates that split a row's bloc partners (congruent columns) from the rest; +1 if the partners are on the + side, -1 if on the - side.
rw_bloc_coords <- function(V_model, bloc, row_names, col_names) {
  C <- length(col_names); out <- numeric(ncol(V_model))
  for (k in seq_len(ncol(V_model))) {
    v <- V_model[, k]; nz <- which(abs(v) > 1e-10); rr <- unique(((nz - 1) %/% C) + 1)
    if (length(rr) != 1) next
    G <- setdiff(which(bloc[col_names] == bloc[row_names[rr]]), rr)
    P <- ((which(v > 1e-10) - 1) %% C) + 1; N <- ((which(v < -1e-10) - 1) %% C) + 1
    if (length(G) > 0 && length(G) < C - 1) { if (setequal(P, G)) out[k] <- 1 else if (setequal(N, G)) out[k] <- -1 }
  }
  out
}

## Fit one table and report. All the choices are arguments; the defaults are the senc configuration.
##   kc          true cells (areas x R x C), used only to score the fit
##   V_ilr       orthonormal basis, row-margin coordinates first
##   alloc       3 = adjusted row, 2 = odds-ratio logit
##   row_order   allocation order of the rows (last = remainder); default smallest row first
##   rem_col     reference column of each row; default column 1. For tables with the same categories on both margins use
##               seq_len(R): each row's own (loyal) column.
##   tiers       NULL: one sigma per coordinate. Otherwise the sigma group of each within-row coordinate (length C - 1), shared by all rows.
##   cores, refresh  passed to rstan::sampling; with cores > 1 rstan shows no progress in RStudio, so for a quick look use
##               chains = 1, iter = 100, warmup = 50, refresh = 10 (prints the gradient time)
##   E_sd_scale, E_sd_scale_small  multiply the E_rc prior sd (all coordinates / coordinates on small columns only)
##   E_centre    "uniform" (the logit of a uniform split) or "qi" (quasi-independence from loyalty estimated from the margins) for the non-loyalty coordinates
##   bloc, bloc_affinity  with bloc (named labels of the categories), adds bloc_affinity (a log odds ratio) to the coordinates that split a row's bloc partners from the rest
##   E_centre_shift  (experiments) moves the non-loyalty centres by this many prior sds, alternating in sign
##   E_sd_scale_offdiag  multiplies the E_rc prior sd of every coordinate except each row's loyalty split
##   ...         passed to ei_estimate and on to rstan::sampling, e.g. pars = ..., include = FALSE, or return_data = TRUE
##   loyalty_mean  NULL: every coordinate's prior centred on the logit of a uniform split. Otherwise the logit mean for the first (loyalty)
##               split of each row, e.g. qlogis(0.75).
rw_fit <- function(kc, V_ilr, alloc = 3L, row_order = NULL, rem_col = NULL, tiers = NULL, loyalty_mean = NULL,
                   chains = 4, iter = 1000, warmup = 500, seed = 1234,
                   cores = chains, refresh = max(iter %/% 10, 1), E_sd_scale = 1, E_sd_scale_small = E_sd_scale, small_frac = 0.05, E_sd_scale_offdiag = 1,
                   E_centre = c("uniform", "qi"), bloc = NULL, bloc_affinity = 0, E_centre_shift = 0, ...) {
  J <- dim(kc)[1]; R <- dim(kc)[2]; C <- dim(kc)[3]
  rm <- apply(kc, c(1, 2), sum); cm <- apply(kc, c(1, 3), sum)      # keep the category names when kc has them
  if (is.null(row_order)) row_order <- order(colSums(rm))
  if (is.null(rem_col)) rem_col <- rep(1L, R)
  mr  <- build_margin_reduction(V_ilr, rm, R = R, C = C)    # row-margin coordinates removed: rows conditioned on the data
  mr  <- rw_setup(mr, rm, cm, row_order = row_order, rem_col = rem_col, alloc = alloc)
  Dm1 <- mr$Dm1_model
  ## with sigma shared in a few groups the scaling has a low-rank form: same density, cheaper gradient (about 1.5x at 7 x 7, 73 areas)
  if (!is.null(tiers)) mr$lflag_mp_scale_fast <- 2L
  ## priors on the average table, one logit per split: the logit implied by a uniform split, with its sd
  m_pos <- colSums(mr$V_ilr_model > 1e-10); n_neg <- colSums(mr$V_ilr_model < -1e-10)
  E_mean <- digamma(m_pos) - digamma(n_neg); E_sd <- sqrt(trigamma(m_pos) + trigamma(n_neg))
  ## optional tightening of the E_rc prior sd: E_sd_scale for every coordinate, E_sd_scale_small for coordinates whose loadings
  ## have a side made only of small columns (column share of the total below small_frac); both default 1 = the uniform-split prior
  if (E_sd_scale != 1 || E_sd_scale_small != 1) {
    small_col <- which(colSums(cm) / sum(cm) < small_frac)
    cell_col <- ((seq_len(nrow(mr$V_ilr_model)) - 1) %% C) + 1
    ## small coordinate: one side of the split consists only of small columns (e.g. a spoilt-ballot or minor-party share)
    is_small <- apply(mr$V_ilr_model, 2, function(v) all(cell_col[v > 1e-10] %in% small_col) || all(cell_col[v < -1e-10] %in% small_col))
    E_sd <- E_sd * ifelse(is_small, E_sd_scale_small, E_sd_scale)
    message(sprintf("E_rc prior sd scaled: %d of %d coordinates small-cell (scale %.2f), the rest %.2f",
                    sum(is_small), Dm1, E_sd_scale_small, E_sd_scale))
  }
  ## E_sd_scale_offdiag: scale the sd of every coordinate except each row's loyalty split (the first split, which involves all C columns),
  ## i.e. the contrasts between the columns a row's voters do NOT stay with
  loy <- (m_pos + n_neg) == C
  if (E_sd_scale_offdiag != 1) E_sd[!loy] <- E_sd[!loy] * E_sd_scale_offdiag
  E_centre <- match.arg(E_centre)
  if (E_centre == "qi") {
    qi <- rw_qi_centre(mr$V_ilr_model, rm, cm)
    use <- !loy & is.finite(qi)
    E_mean[use] <- qi[use]
    message(sprintf("E_rc prior centred on quasi-independence for %d of %d non-loyalty coordinates", sum(use), sum(!loy)))
  }
  if (bloc_affinity != 0) {
    if (is.null(bloc)) stop("bloc_affinity needs bloc (a named vector of bloc labels)")
    bc <- rw_bloc_coords(mr$V_ilr_model, bloc, dimnames(kc)[[2]], dimnames(kc)[[3]])
    E_mean <- E_mean + bloc_affinity * bc
    message(sprintf("bloc affinity %.2f added to %d coordinates", bloc_affinity, sum(bc != 0)))
  }
  ## for experiments: move the non-loyalty centres away by E_centre_shift prior sds (the sd actually used), alternating in sign, to create
  ## a deliberate conflict between prior and data
  if (E_centre_shift != 0) {
    sgn <- rep(c(1, -1), length.out = length(E_mean))
    E_mean[!loy] <- E_mean[!loy] + E_centre_shift * E_sd[!loy] * sgn[!loy]
    message(sprintf("non-loyalty centres shifted by %.1f prior sds", E_centre_shift))
  }
  if (!is.null(loyalty_mean)) {                              # the loyalty split is the first split of each row: it involves all C columns
    E_mean[loy] <- ifelse(m_pos[loy] == 1, 1, -1) * loyalty_mean
  }
  sig_id <- if (is.null(tiers)) seq_len(Dm1) else rep(as.integer(tiers), R)
  G <- max(sig_id)
  stopifnot(length(sig_id) == Dm1)
  fit <- ssEI::ei_estimate(
    rm, cm,
    E_rc_fixed = rep(0, R * C - 1), sigma_jrc_fixed = rep(.5, R * C - 1),
    fix_E_rc = 0, fix_sigma_jrc = 0,
    E_rc_prior = rep(0, R * C - 1),
    known_cell_values = kc, use_known_cells = 0,
    V_ilr = V_ilr, n_ilr_rows = R,
    fit_type = "soft multinom",
    row_decompose = TRUE,
    rotate_llrep = TRUE,                                    # ignored by the margin model
    margin_reduction = mr,
    ROT_E_rc = diag(Dm1), rotate_E_rc = FALSE,
    link_E_rc = NULL, rotate_lambda = "none",
    sigma_group_id = sig_id, sigma_group_mode = rep("shared", G),             # sigma shared by the areas
    sigma_group_prior_a = rep(log(.3), G), sigma_group_prior_b = rep(0.5, G),
    sigma_ncp = TRUE,
    E_rc_group_id = 1:Dm1, E_rc_group_mode = rep("shared", Dm1),
    E_rc_group_prior_a = unname(E_mean), E_rc_group_prior_b = unname(E_sd),
    E_rc_node_logit = TRUE,
    neutral_logit = "llrep", lambda_centred = FALSE,
    sigma_floor = 0, prior_sigma_c_scale = 2, prior_lambda_raw_scale = 12,
    prior_gamma_shape = 2, prior_gamma_rate = .5,
    raw_seq_cell_weights = TRUE,
    chains = chains, cores = cores, refresh = refresh, iter = iter, warmup = warmup,
    init = rw_init(J, mr$n_pin, (R - 1) * (C - 1), E_mean, G), seed = seed, ...
  )
  ## return_data = TRUE: the Stan data and the starting-value function, for running the model elsewhere (e.g. cmdstanr)
  if (isTRUE(list(...)$return_data)) return(list(data = fit, init = rw_init(J, mr$n_pin, (R - 1) * (C - 1), E_mean, G)))
  if (inherits(fit, "stanfit")) rw_report(fit, kc)
  invisible(fit)
}

## Error index, interval coverage, the sigma and E_rc summaries, and what the sampler did.
rw_report <- function(fit, kc) {
  cv  <- rstan::extract(fit, "cell_values")$cell_values                  # draws x areas x R x C
  est <- apply(cv, 2:4, mean)
  lo  <- apply(cv, 2:4, quantile, 0.05); hi <- apply(cv, 2:4, quantile, 0.95)
  big <- kc >= 50
  cat(sprintf("error index %.2f | 90%% coverage of cells of 50+ %.2f\n",
              50 * sum(abs(est - kc)) / sum(kc), mean((kc >= lo & kc <= hi)[big])))
  s <- rstan::summary(fit, pars = c("E_rc", "sigma_group_mu"))$summary
  if (nrow(s) <= 12) print(round(s[, c("mean", "sd", "n_eff", "Rhat")], 3)) else {
    sg <- s[grepl("^sigma_group_mu", rownames(s)), , drop = FALSE]; print(round(sg[, c("mean", "sd", "n_eff", "Rhat")], 3))
    cat(sprintf("E_rc (%d values): min n_eff %.0f, max Rhat %.3f\n", sum(grepl("^E_rc", rownames(s))),
                min(s[grepl("^E_rc", rownames(s)), "n_eff"]), max(s[grepl("^E_rc", rownames(s)), "Rhat"])))
  }
  cat(sprintf("min n_eff over all of these %.0f | max Rhat %.3f\n", min(s[, "n_eff"]), max(s[, "Rhat"])))
  sp <- do.call(rbind, rstan::get_sampler_params(fit, inc_warmup = FALSE))
  cat(sprintf("step size %.3f | leapfrogs per draw %.0f | divergences %d | draws at maximum tree depth %d\n",
              mean(sp[, "stepsize__"]), mean(sp[, "n_leapfrog__"]), sum(sp[, "divergent__"]),
              sum(sp[, "treedepth__"] >= 10)))
  print(rstan::get_elapsed_time(fit))
}
