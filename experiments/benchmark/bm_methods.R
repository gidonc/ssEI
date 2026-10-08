## Competitor methods for the benchmark. Every wrapper takes the area margins and returns the same list:
##   method, est (J x R x C), draws (S x J x R x C or NULL), chain (chain of each draw), conv (bm_conv over every cell, all chains),
##   secs, settings, attempts (one row per run when a run is repeated longer until it converges)
##
## Orientation: rm = row margins (J x R), cm = column margins (J x C); cell [j, r, c] is the number of voters in area j in row r and
## column c. Row r is the "predictor" of the formula `cbind(columns) ~ cbind(rows)` in eiPack and RxCEcolInf, and `votes_election1`
## (origin) in lphom. Categories are renamed r1.. and c1.. so that long party names cannot upset any formula.
##
## Empty rows: nothing is removed. A row with no voters in an area is passed to every method as a zero margin; each method then
## decides what to do with it, and the wrapper checks that its cells in that row come out as 0.
##
## Needs eiPack, RxCEcolInf, lphom, coda, posterior.  Chains run in parallel with parallel::mclapply (not on Windows).

source_dir <- if (exists("BM_DIR")) BM_DIR else "experiments/benchmark"
if (!exists("bm_conv")) source(file.path(source_dir, "bm_score.R"))

## cores for the chains: what a SLURM job was given, else what the machine has, never more than one per chain
bm_cores <- function(chains) {
  slurm <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", NA)))
  min(chains, if (is.na(slurm)) parallel::detectCores() else slurm)
}

bm_codes <- function(rm, cm) {
  rm <- as.matrix(rm); cm <- as.matrix(cm); dimnames(rm) <- list(NULL, paste0("r", seq_len(ncol(rm))))
  dimnames(cm) <- list(NULL, paste0("c", seq_len(ncol(cm))))
  stopifnot(all(rowSums(rm) == rowSums(cm)))
  list(rm = rm, cm = cm, data = data.frame(cm, rm),
       formula = as.formula(sprintf("cbind(%s) ~ cbind(%s)", paste(colnames(cm), collapse = ","), paste(colnames(rm), collapse = ","))),
       fstring = sprintf("%s~%s", paste(colnames(cm), collapse = ","), paste(colnames(rm), collapse = ",")))
}

## cells in an empty row must be 0
bm_check_empty_rows <- function(est, rm) {
  empty <- as.matrix(rm) == 0
  if (!any(empty)) return(0)
  bad <- max(abs(est[rep(empty, times = dim(est)[3])]))     # empty is J x R; recycled over the C columns
  if (bad > 1e-8) warning(sprintf("cells of empty rows are not 0 (largest %.3g)", bad))
  bad
}

## ---------------------------------------------------------------------------
## Baseline: independence within each area, cell = row total x column total / area total
## ---------------------------------------------------------------------------
bm_independence <- function(rm, cm) {
  rm <- as.matrix(rm); cm <- as.matrix(cm); J <- nrow(rm)
  est <- array(0, c(J, ncol(rm), ncol(cm)))
  for (j in seq_len(J)) est[j, , ] <- outer(rm[j, ], cm[j, ]) / sum(rm[j, ])
  list(method = "independence", est = est, draws = NULL, chain = NULL, conv = NULL, secs = 0, settings = list(), attempts = NULL)
}

## ---------------------------------------------------------------------------
## nslphom (Pavia and Romero): linear programming under homogeneity; no draws
## ---------------------------------------------------------------------------
bm_nslphom <- function(rm, cm, iter.max = 10, ...) {
  d <- bm_codes(rm, cm); t0 <- Sys.time()
  fit <- lphom::nslphom(d$rm, d$cm, new_and_exit_voters = "raw", iter.max = iter.max, verbose = FALSE, ...)
  u <- fit$VTM.votes.units                                     # origin x destination x unit, in votes
  stopifnot(identical(dim(u), c(ncol(d$rm), ncol(d$cm), nrow(d$rm))))
  est <- aperm(u, c(3, 1, 2)); dimnames(est) <- NULL
  ## the estimates must reproduce the margins they were given, or the array is the wrong way round
  dev <- max(abs(apply(est, c(1, 2), sum) - d$rm), abs(apply(est, c(1, 3), sum) - d$cm))
  if (dev > 1e-6 * max(d$rm)) warning(sprintf("nslphom units do not reproduce the margins (largest gap %.3g)", dev))
  list(method = "nslphom", est = est, draws = NULL, chain = NULL, conv = NULL,
       secs = as.numeric(Sys.time() - t0, units = "secs"), settings = list(iter.max = iter.max, margin_gap = dev, empty_row_max = bm_check_empty_rows(est, d$rm)),
       attempts = NULL)
}

## ---------------------------------------------------------------------------
## ei.MD.bayes (eiPack): hierarchical multinomial-Dirichlet, Rosen et al. 2001
## ---------------------------------------------------------------------------
## Defaults are the "manual" settings of Pavia and Romero (2023): tuneMD with 10 x 100000 draws, then 1000 draws, thin 100, burn-in
## 100000. Chains are separate runs from random starting values. If convergence over all cells and alphas fails, the run is repeated
## with thin and burn-in doubled, up to max_doublings times (the tuning is kept).
bm_eiMD <- function(rm, cm, chains = 4, sample = 1000, thin = 100, burnin = 1e5, ntunes = 10, tune_draws = 1e5,
                    max_doublings = 0, cores = bm_cores(chains), seed = 1, ...) {
  d <- bm_codes(rm, cm); J <- nrow(d$rm); R <- ncol(d$rm); C <- ncol(d$cm); t_start <- Sys.time()
  tl <- eiPack::tuneMD(d$formula, data = d$data, ntunes = ntunes, totaldraws = tune_draws)
  secs_tune <- as.numeric(Sys.time() - t_start, units = "secs")
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed)
  attempts <- NULL
  for (k in 0:max_doublings) {
    sc <- 2^k; t0 <- Sys.time()
    runs <- parallel::mclapply(seq_len(chains), function(i)
      eiPack::ei.MD.bayes(d$formula, data = d$data, sample = sample, thin = thin * sc, burnin = burnin * sc, tune.list = tl, ...),
      mc.cores = cores, mc.set.seed = TRUE)
    bad <- vapply(runs, function(x) is.null(x) || inherits(x, "try-error"), NA)
    if (any(bad)) stop("ei.MD.bayes failed in chain ", which(bad)[1], ": ", if (is.null(runs[[which(bad)[1]]])) "the worker died (out of memory?)" else runs[[which(bad)[1]]])
    ## cell counts = within-row proportion x row total
    P <- J * R * C; S <- sample * chains; D <- matrix(0, S, P); A <- NULL
    for (i in seq_len(chains)) {
      B <- as.matrix(runs[[i]]$draws$Beta); m <- regmatches(colnames(B), regexec("^beta\\.r([0-9]+)\\.c([0-9]+)\\.([0-9]+)$", colnames(B)))
      stopifnot(all(lengths(m) == 4)); ix <- t(vapply(m, function(x) as.integer(x[2:4]), integer(3)))   # row, col, unit
      cell <- ix[, 3] + J * (ix[, 1] - 1) + J * R * (ix[, 2] - 1)                  # position in a J x R x C array
      D[(i - 1) * sample + seq_len(sample), cell] <- sweep(B, 2, d$rm[cbind(ix[, 3], ix[, 1])], "*")
      a <- as.matrix(runs[[i]]$draws$Alpha); if (is.null(A)) A <- array(NA_real_, c(sample, chains, ncol(a)), dimnames = list(NULL, NULL, colnames(a)))
      A[, i, ] <- a
    }
    draws <- D; dim(draws) <- c(S, J, R, C); chain <- rep(seq_len(chains), each = sample)
    ca <- bm_chain_array(draws, chain)
    conv <- bm_conv(abind::abind(ca, A, along = 3))             # every cell count and every alpha
    attempts <- rbind(attempts, data.frame(scale = sc, thin = thin * sc, burnin = burnin * sc, secs = as.numeric(Sys.time() - t0, units = "secs"),
                                           max_rhat = conv$max_rhat, min_ess_bulk = conv$min_ess_bulk, converged = conv$converged))
    if (conv$converged) break
  }
  est <- apply(draws, 2:4, mean)
  list(method = "ei.MD.bayes", est = est, draws = draws, chain = chain, conv = conv,
       secs = as.numeric(Sys.time() - t_start, units = "secs"),
       settings = list(chains = chains, sample = sample, thin = thin * sc, burnin = burnin * sc, ntunes = ntunes, tune_draws = tune_draws,
                       secs_tune = secs_tune, empty_row_max = bm_check_empty_rows(est, d$rm)), attempts = attempts)
}

## ---------------------------------------------------------------------------
## RxCEcolInf: Greiner and Quinn (2009), hierarchical model with latent cell counts
## ---------------------------------------------------------------------------
## `keep` cell-count draws are saved per chain over num.iters = keep * thin iterations; the first burnin_frac of them are discarded.
## Defaults (keep 1000, thin 1500, burn-in 10%) are 1.5 million iterations, as in the old development scripts. As for ei.MD.bayes, a
## run that has not converged is repeated with thin doubled.
bm_rxc <- function(rm, cm, chains = 4, keep = 1000, thin = 1500, burnin_frac = 0.1, tune_iters = 2e4, tune_runs = 15,
                   max_doublings = 0, cores = bm_cores(chains), seed = 1, ...) {
  d <- bm_codes(rm, cm); J <- nrow(d$rm); R <- ncol(d$rm); C <- ncol(d$cm); t_start <- Sys.time()
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed)
  invisible(utils::capture.output(tn <- RxCEcolInf::Tune(d$fstring, data = d$data, num.iters = tune_iters, num.runs = tune_runs, debug = 0)))
  secs_tune <- as.numeric(Sys.time() - t_start, units = "secs"); attempts <- NULL
  for (k in 0:max_doublings) {
    sc <- 2^k; t0 <- Sys.time(); n_it <- keep * thin * sc
    runs <- parallel::mclapply(seq_len(chains), function(i) {
      invisible(utils::capture.output(ch <- RxCEcolInf::Analyze(d$fstring, rho.vec = tn$rhos, data = d$data, num.iters = n_it,
                                                                burnin = round(burnin_frac * n_it), save.every = max(n_it %/% 1000, 1),
                                                                keepNNinternals = keep, keepTHETAS = 0, debug = 0, ...)))
      attr(ch, "NN.internals")
    }, mc.cores = cores, mc.set.seed = TRUE)
    bad <- vapply(runs, function(x) is.null(x) || inherits(x, "try-error"), NA)
    if (any(bad)) stop("Analyze failed in chain ", which(bad)[1], ": ", if (is.null(runs[[which(bad)[1]]])) "the worker died (out of memory?)" else runs[[which(bad)[1]]])
    n_s <- nrow(runs[[1]]); S <- n_s * chains; D <- matrix(0, S, J * R * C)
    for (i in seq_len(chains)) {
      N <- as.matrix(runs[[i]]); stopifnot(nrow(N) == n_s)
      m <- regmatches(colnames(N), regexec("^NN\\.table ([0-9]+) r([0-9]+)\\.c([0-9]+)$", colnames(N)))
      stopifnot(all(lengths(m) == 4)); ix <- t(vapply(m, function(x) as.integer(x[2:4]), integer(3)))   # unit, row, col
      D[(i - 1) * n_s + seq_len(n_s), ix[, 1] + J * (ix[, 2] - 1) + J * R * (ix[, 3] - 1)] <- N
    }
    draws <- D; dim(draws) <- c(S, J, R, C); chain <- rep(seq_len(chains), each = n_s)
    conv <- bm_conv(bm_chain_array(draws, chain))
    attempts <- rbind(attempts, data.frame(scale = sc, thin = thin * sc, burnin = round(burnin_frac * n_it), secs = as.numeric(Sys.time() - t0, units = "secs"),
                                           max_rhat = conv$max_rhat, min_ess_bulk = conv$min_ess_bulk, converged = conv$converged))
    if (conv$converged) break
  }
  est <- apply(draws, 2:4, mean)
  list(method = "RxCEcolInf", est = est, draws = draws, chain = chain, conv = conv,
       secs = as.numeric(Sys.time() - t_start, units = "secs"),
       settings = list(chains = chains, keep = keep, thin = thin * sc, num_iters = n_it, burnin_frac = burnin_frac, tune_iters = tune_iters,
                       tune_runs = tune_runs, secs_tune = secs_tune, empty_row_max = bm_check_empty_rows(est, d$rm)), attempts = attempts)
}

## ---------------------------------------------------------------------------
## ssEI: the same list from a stanfit returned by rw_fit(), so it goes through the same scoring and convergence summary
## ---------------------------------------------------------------------------
bm_ssei <- function(fit, rm = NULL, label = "ssEI") {
  b <- bm_from_stanfit(fit); est <- apply(b$draws, 2:4, mean)
  list(method = label, est = est, draws = b$draws, chain = b$chain, conv = bm_conv(bm_chain_array(b$draws, b$chain)),
       secs = sum(rstan::get_elapsed_time(fit)), settings = list(empty_row_max = if (is.null(rm)) NA else bm_check_empty_rows(est, rm)),
       attempts = NULL)
}
