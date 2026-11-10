## Sign-restricted contextual effects (monotone-EI style): does stating only the DIRECTION in which a row's split changes with the
## row's local share give honest intervals for the summed (national) table, and what does it do to accuracy?
##
## Usage, from the package directory:  Rscript experiments/structure/sign_context.R <data> <variant> [chains iter warmup]
##   data     senc | nz2020 | scotland7
##   variant  base         no contextual effects
##            free         covariate = standardised log share of the coordinate's own row, on the chosen coordinates; coefficients unrestricted
##            sign         the same with the stated signs
##            sign_fixed   signs, coefficients half-normal(0, 0.5) (no estimated sd)
##            sign_wide    signs, coefficients half-normal(0, 1)
##            smallcell    the small-cell term, no contextual effects; smallcell_free = the term plus `free`
##            all_free     senc: every coordinate, unrestricted (the 13 October run)
##            all_true     every coordinate with a covariate, signs taken from the true tables (a ceiling for sign information)
## Chosen coordinates and stated signs:
##   senc      each row's dem-against-rep split: white <= 0 (more Republican where whites are concentrated), black and natam >= 0
##             (more Democratic where the group is concentrated): contextual reinforcement. The {dem, rep}-against-non split has no covariate.
##   nz2020, scotland7   each loyalty split, <= 0: a candidate's voters are no more loyal on the party vote where the candidate's share
##             is high (personal and tactical votes come from other parties' supporters).
## Set-ups as in claude/ssei_national_totals.md: soft multinomial, rows conditioned, no blocs, QI centre, off-diagonal prior scale 0.5,
## all-party loyalty (NZ) / major-party loyalty (Scotland); senc with two sigma tiers. Results are appended to
## experiments/structure/results/sign_context.csv; draws of the summed table and the coefficients are saved next to it.

suppressMessages({ source("experiments/rowwise_helpers.R"); source("experiments/data_prep.R") })
a <- commandArgs(TRUE); DATA <- a[1]; VARIANT <- a[2]
CHAINS <- if (length(a) >= 3) as.integer(a[3]) else 2L; ITER <- if (length(a) >= 4) as.integer(a[4]) else 1000L; WARMUP <- if (length(a) >= 5) as.integer(a[5]) else 500L
OUT <- "experiments/structure/results"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## true slopes of the within-row model coordinates (0.5 added to cells) on the covariate, weighted by row size
true_slopes <- function(kc, V_ilr, X) {
  R <- dim(kc)[2]; C <- dim(kc)[3]; Vw <- V_ilr[, R:(R * C - 1), drop = FALSE]; rm <- apply(kc, 1:2, sum)
  sapply(seq_len(ncol(Vw)), function(k) { r <- ((which(abs(Vw[, k]) > 1e-10)[1] - 1) %/% C) + 1; ok <- rm[, r] > 0
    if (all(X[ok, k] == 0)) return(0)
    coef(lm(as.vector(log(kc[ok, r, ] + 0.5) %*% Vw[(r - 1) * C + 1:C, k]) ~ X[ok, k], weights = rm[ok, r]))[2] })
}

if (DATA == "senc") {
  data(senc, package = "eiPack"); rown <- c("white", "black", "natam"); coln <- c("dem", "rep", "non"); R <- 3; C <- 3
  kc <- array(0, c(nrow(senc), R, C), dimnames = list(NULL, rown, coln)); cellnm <- outer(c("wh", "bl", "natam"), coln, paste0)
  for (r in 1:R) for (c in 1:C) kc[, r, c] <- senc[[cellnm[r, c]]]
  Vr <- rw_sbp_to_v(cbind(c(1, 1, -1), c(1, -1, 0))); V_within <- matrix(0, R * C, R * (C - 1))
  for (r in 1:R) V_within[(r - 1) * C + 1:C, (r - 1) * (C - 1) + 1:(C - 1)] <- Vr
  V_ilr <- cbind(kronecker(make_helmert_basis(R), matrix(1 / sqrt(C), C, 1)), V_within)
  base_args <- list(kc, V_ilr, alloc = 3L, rem_col = c(1L, 1L, 1L), tiers = c(1, 2))
  coords <- c(2, 4, 6); signs <- c(0, -1, 0, 1, 0, 1); ctx <- "own_row"
} else {
  if (DATA == "nz2020") { d <- prep_nz(2020, 6, "candidate", 1e6); loyal_rows <- NULL }
  if (DATA == "scotland7") { d <- prep_scotland(7, 73)
    loyal_rows <- c("Conservative and Unionist Party [The]", "Labour Party [The]", "Liberal Democrats", "Scottish National Party") }
  R <- length(d$row_names); C <- length(d$col_names); kc <- d$kc; dimnames(kc) <- list(NULL, d$row_names, d$col_names)
  V_ilr <- rw_row_priority_basis(d$row_names, d$col_names, NULL, bisect_order = colSums(d$cm), loyal_rows = loyal_rows)
  base_args <- list(kc, V_ilr, alloc = 3L, rem_col = seq_len(R), tiers = c(1, 1, rep(2, C - 3)), loyalty_mean = qlogis(0.75),
                    loyal_rows = loyal_rows, E_sd_scale_offdiag = 0.5, E_centre = "qi")
  coords <- NULL; signs <- -1; ctx <- "own_row_loyalty"
}
J <- dim(kc)[1]; Dm1 <- R * (C - 1)

## the covariate matrix as rw_fit builds it for "own_row" (needed for the true slopes)
rm <- apply(kc, 1:2, sum)
xr <- sapply(seq_len(R), function(r) { ok <- rm[, r] > 0; x <- numeric(J); l <- log(rm[ok, r] / rowSums(rm)[ok]); if (sd(l) > 0) x[ok] <- (l - mean(l)) / sd(l); x })
X_all <- xr[, rep(seq_len(R), each = C - 1)]
b_true_all <- true_slopes(kc, V_ilr, X_all)

extra <- switch(VARIANT,
  base       = list(),
  free       = list(context = ctx, context_coords = coords),
  sign       = list(context = ctx, context_coords = coords, context_sign = signs),
  sign_fixed = list(context = ctx, context_coords = coords, context_sign = signs, context_tau = FALSE, context_scale = 0.5),
  sign_wide  = list(context = ctx, context_coords = coords, context_sign = signs, context_tau = FALSE, context_scale = 1),
  smallcell  = list(smallcell = 1),
  smallcell_free = list(smallcell = 1, context = ctx, context_coords = coords),
  all_free   = list(context = "own_row"),
  all_true   = list(context = "own_row", context_sign = sign(b_true_all)),
  stop("unknown variant"))
t0 <- Sys.time()
fit <- do.call(rw_fit, c(base_args, extra, list(chains = CHAINS, iter = ITER, warmup = WARMUP, seed = 1234, refresh = 250)))
mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

## scores
cv <- rstan::extract(fit, "cell_values")$cell_values; est <- apply(cv, 2:4, mean)
lo <- apply(cv, 2:4, quantile, 0.05); hi <- apply(cv, 2:4, quantile, 0.95); big <- kc >= 50
tot <- apply(cv, c(1, 3, 4), sum); kt <- apply(kc, 2:3, sum); et <- apply(tot, 2:3, mean)
lt <- apply(tot, 2:3, quantile, 0.05); ht <- apply(tot, 2:3, quantile, 0.95)
rs <- rowSums(kt); width <- 100 * mean(((ht - lt) / rs)[kt > 0])                    # mean 90% interval width of the summed table's row rates, points
sp <- do.call(rbind, rstan::get_sampler_params(fit, inc_warmup = FALSE))
pars <- c("E_rc", "sigma_group_mu", if (!VARIANT %in% c("base", "smallcell")) "E_cov_b")
sm <- rstan::summary(fit, pars = pars)$summary
b_est <- if (!VARIANT %in% c("base", "smallcell")) rstan::summary(fit, pars = "E_cov_b")$summary[, c("mean", "2.5%", "97.5%")] else NULL
res <- data.frame(data = DATA, variant = VARIANT, chains = CHAINS, iter = ITER,
  ei_area = 50 * sum(abs(est - kc)) / sum(kc), ei_total = 50 * sum(abs(et - kt)) / sum(kt),
  cover_area = mean((kc >= lo & kc <= hi)[big]), cover_total = mean((kt >= lt & kt <= ht)[kt > 0]), width_total = width,
  divergences = sum(sp[, "divergent__"]), leapfrogs = mean(sp[, "n_leapfrog__"]), min_neff = min(sm[, "n_eff"], na.rm = TRUE),
  max_rhat = max(sm[, "Rhat"], na.rm = TRUE), minutes = mins)
f <- file.path(OUT, "sign_context.csv")
write.table(res, f, sep = ",", row.names = FALSE, col.names = !file.exists(f), append = file.exists(f))
saveRDS(list(res = res, tot = tot, truth = kt, est = est, b_est = b_est, b_true = b_true_all, summary = sm), file.path(OUT, sprintf("sign_context_%s_%s.rds", DATA, VARIANT)))
print(res, row.names = FALSE)
if (!is.null(b_est)) { cat("coefficients (posterior mean, 95% interval) against the true slopes\n"); print(round(cbind(b_est, true = b_true_all), 3)) }
cat("summed table, row rates: truth / estimate / 90% interval\n")
for (r in seq_len(R)) cat(sprintf("%-10s %s\n", substr(dimnames(kc)[[2]][r], 1, 10),
  paste(sprintf("%.3f/%.3f[%.3f,%.3f]", kt[r, ] / rs[r], et[r, ] / rs[r], lt[r, ] / rs[r], ht[r, ] / rs[r]), collapse = "  ")))
