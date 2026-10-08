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
## Extending a run that has not converged
## ---------------------------------------------------------------------------
## ei.MD.bayes and RxCEcolInf results carry `state`: what is needed to carry each chain on from where it stopped. Passing a saved result as
## `resume=` runs another segment of draws for every chain from that state and appends it, then convergence is judged on all draws together.
## `max_extensions` does the same automatically: after a run that has not converged, up to that many more segments are added, so a
## non-converged run is continued, never restarted.
##   ei.MD.bayes: exact. start.list (final Alpha and Beta) plus the same tuning, no burn-in.
##   RxCEcolInf:  approximate. Only the latent cell counts and mu can be restored (THETAS.start is broken in the package, Sigma cannot be
##                set), so each extension re-burns in for `resume_burnin_frac` of its iterations and drops those draws.
bm_split_chains <- function(draws, chain) {            # S x J x R x C and chain label -> list of (draws in chain x cells) matrices
  D <- draws; dim(D) <- c(dim(draws)[1], prod(dim(draws)[-1]))
  lapply(sort(unique(chain)), function(i) D[chain == i, , drop = FALSE])
}
bm_join_chains <- function(parts, dims) {              # inverse of bm_split_chains; dims = c(J, R, C)
  D <- do.call(rbind, parts); draws <- D; dim(draws) <- c(nrow(D), dims)
  list(draws = draws, chain = rep(seq_along(parts), times = vapply(parts, nrow, 1L)))
}

## ---------------------------------------------------------------------------
## ei.MD.bayes (eiPack): hierarchical multinomial-Dirichlet, Rosen et al. 2001
## ---------------------------------------------------------------------------
## Defaults are the "manual" settings of Pavia and Romero (2023): tuneMD with 10 x 100000 draws, then 1000 draws, thin 100, burn-in
## 100000. Chains are separate runs from random starting values. Each extension adds `sample` more draws per chain at the same thin.
bm_eiMD <- function(rm, cm, chains = 4, sample = 1000, thin = 100, burnin = 1e5, ntunes = 10, tune_draws = 1e5,
                    max_extensions = 0, resume = NULL, cores = bm_cores(chains), seed = 1, ...) {
  d <- bm_codes(rm, cm); J <- nrow(d$rm); R <- ncol(d$rm); C <- ncol(d$cm); t_start <- Sys.time(); RNGkind("L'Ecuyer-CMRG")
  if (is.null(resume)) {
    set.seed(seed); tl <- eiPack::tuneMD(d$formula, data = d$data, ntunes = ntunes, totaldraws = tune_draws)
    secs_tune <- as.numeric(Sys.time() - t_start, units = "secs"); last <- NULL; parts <- NULL; Aold <- NULL; attempts <- NULL; secs_before <- 0
    n_run <- 1 + max_extensions
  } else {
    stopifnot(identical(resume$method, "ei.MD.bayes"), identical(dim(resume$est), c(J, R, C)), !is.null(resume$state))
    tl <- resume$state$tune; last <- resume$state$last; chains <- length(last); thin <- resume$settings$thin; burnin <- resume$settings$burnin
    ntunes <- resume$settings$ntunes; tune_draws <- resume$settings$tune_draws; secs_tune <- resume$settings$secs_tune
    parts <- bm_split_chains(resume$draws, resume$chain); Aold <- resume$state$alpha; attempts <- resume$attempts
    secs_before <- resume$secs; n_run <- max(max_extensions, 1); set.seed(seed + max(attempts$segment))
  }
  state_of <- function(fit) {                          # final Alpha (R x C) and Beta (R x C x J), as start.list wants them
    B <- as.matrix(fit$draws$Beta); A <- as.matrix(fit$draws$Alpha); i <- nrow(B)
    m <- regmatches(colnames(B), regexec("^beta\\.r([0-9]+)\\.c([0-9]+)\\.([0-9]+)$", colnames(B)))
    ix <- t(vapply(m, function(x) as.integer(x[2:4]), integer(3))); beta <- array(NA_real_, c(R, C, J)); beta[ix] <- B[i, ]
    ma <- regmatches(colnames(A), regexec("^alpha\\.r([0-9]+)\\.c([0-9]+)$", colnames(A))); ia <- t(vapply(ma, function(x) as.integer(x[2:3]), integer(2)))
    alpha <- matrix(NA_real_, R, C); alpha[ia] <- A[i, ]
    list(start.alphas = alpha, start.betas = beta)
  }
  for (k in seq_len(n_run)) {
    t0 <- Sys.time(); first <- is.null(last)
    runs <- parallel::mclapply(seq_len(chains), function(i) {
      a <- list(d$formula, data = d$data, sample = sample, thin = thin, burnin = if (first) burnin else 0, tune.list = tl, ...)
      if (!first) a$start.list <- last[[i]]
      do.call(eiPack::ei.MD.bayes, a)
    }, mc.cores = cores, mc.set.seed = TRUE)
    bad <- vapply(runs, function(x) is.null(x) || inherits(x, "try-error"), NA)
    if (any(bad)) stop("ei.MD.bayes failed in chain ", which(bad)[1], ": ", if (is.null(runs[[which(bad)[1]]])) "the worker died (out of memory?)" else runs[[which(bad)[1]]])
    new_parts <- list(); Anew <- NULL
    for (i in seq_len(chains)) {                       # cell counts = within-row proportion x row total
      B <- as.matrix(runs[[i]]$draws$Beta); m <- regmatches(colnames(B), regexec("^beta\\.r([0-9]+)\\.c([0-9]+)\\.([0-9]+)$", colnames(B)))
      stopifnot(all(lengths(m) == 4)); ix <- t(vapply(m, function(x) as.integer(x[2:4]), integer(3)))   # row, col, unit
      cell <- ix[, 3] + J * (ix[, 1] - 1) + J * R * (ix[, 2] - 1)                  # position in a J x R x C array
      Di <- matrix(0, nrow(B), J * R * C); Di[, cell] <- sweep(B, 2, d$rm[cbind(ix[, 3], ix[, 1])], "*"); new_parts[[i]] <- Di
      a <- as.matrix(runs[[i]]$draws$Alpha); if (is.null(Anew)) Anew <- array(NA_real_, c(nrow(a), chains, ncol(a)), dimnames = list(NULL, NULL, colnames(a)))
      Anew[, i, ] <- a
    }
    last <- lapply(runs, state_of)
    parts <- if (is.null(parts)) new_parts else Map(rbind, parts, new_parts)
    Aold <- if (is.null(Aold)) Anew else abind::abind(Aold, Anew, along = 1)
    j <- bm_join_chains(parts, c(J, R, C)); draws <- j$draws; chain <- j$chain
    conv <- bm_conv(abind::abind(bm_chain_array(draws, chain), Aold, along = 3))      # every cell count and every alpha
    attempts <- rbind(attempts, data.frame(segment = if (is.null(attempts)) 1L else max(attempts$segment) + 1L, draws_per_chain = nrow(parts[[1]]),
                                           secs = as.numeric(Sys.time() - t0, units = "secs"), max_rhat = conv$max_rhat,
                                           min_ess_bulk = conv$min_ess_bulk, converged = conv$converged))
    if (conv$converged) break
  }
  est <- apply(draws, 2:4, mean)
  list(method = "ei.MD.bayes", est = est, draws = draws, chain = chain, conv = conv,
       secs = secs_before + as.numeric(Sys.time() - t_start, units = "secs"),
       settings = list(chains = chains, sample = sample, thin = thin, burnin = burnin, ntunes = ntunes, tune_draws = tune_draws,
                       secs_tune = secs_tune, empty_row_max = bm_check_empty_rows(est, d$rm)),
       attempts = attempts, state = list(tune = tl, last = last, alpha = Aold))
}

## ---------------------------------------------------------------------------
## RxCEcolInf: Greiner and Quinn (2009), hierarchical model with latent cell counts
## ---------------------------------------------------------------------------
## Each segment saves `keep` cell-count draws per chain over num.iters = keep * thin iterations; the first burnin_frac of the first
## segment (resume_burnin_frac of an extension) is discarded. Defaults (keep 1000, thin 1500, burn-in 10%) are 1.5 million iterations,
## as in the old development scripts.
bm_rxc <- function(rm, cm, chains = 4, keep = 1000, thin = 1500, burnin_frac = 0.1, resume_burnin_frac = 0.05, tune_iters = 2e4, tune_runs = 15,
                   max_extensions = 0, resume = NULL, cores = bm_cores(chains), seed = 1, ...) {
  d <- bm_codes(rm, cm); J <- nrow(d$rm); R <- ncol(d$rm); C <- ncol(d$cm); t_start <- Sys.time(); RNGkind("L'Ecuyer-CMRG")
  if (is.null(resume)) {
    set.seed(seed)
    invisible(utils::capture.output(tn <- RxCEcolInf::Tune(d$fstring, data = d$data, num.iters = tune_iters, num.runs = tune_runs, debug = 0)))
    secs_tune <- as.numeric(Sys.time() - t_start, units = "secs"); rhos <- tn$rhos; last <- NULL; parts <- NULL; attempts <- NULL; secs_before <- 0
    n_run <- 1 + max_extensions
  } else {
    stopifnot(identical(resume$method, "RxCEcolInf"), identical(dim(resume$est), c(J, R, C)), !is.null(resume$state))
    rhos <- resume$state$rhos; last <- resume$state$last; chains <- length(last); thin <- resume$settings$thin; keep <- resume$settings$keep
    tune_iters <- resume$settings$tune_iters; tune_runs <- resume$settings$tune_runs; secs_tune <- resume$settings$secs_tune
    burnin_frac <- resume$settings$burnin_frac
    parts <- bm_split_chains(resume$draws, resume$chain); attempts <- resume$attempts; secs_before <- resume$secs
    n_run <- max(max_extensions, 1); set.seed(seed + max(attempts$segment))
  }
  n_it <- keep * thin
  for (k in seq_len(n_run)) {
    t0 <- Sys.time(); first <- is.null(last); bf <- if (first) burnin_frac else resume_burnin_frac
    runs <- parallel::mclapply(seq_len(chains), function(i) {
      a <- list(d$fstring, rho.vec = rhos, data = d$data, num.iters = n_it, burnin = round(bf * n_it), save.every = max(n_it %/% 1000, 1),
                keepNNinternals = keep, keepTHETAS = 0, debug = 0, ...)
      if (!first) { a$NNs.start <- list(last[[i]]$nn); a$mu.vec.cu <- last[[i]]$mu }
      invisible(utils::capture.output(ch <- do.call(RxCEcolInf::Analyze, a)))
      list(nn = attr(ch, "NN.internals"), mu = { M <- as.matrix(ch)[, grep("^mu", colnames(as.matrix(ch))), drop = FALSE]; as.numeric(M[nrow(M), ]) })
    }, mc.cores = cores, mc.set.seed = TRUE)
    bad <- vapply(runs, function(x) is.null(x) || inherits(x, "try-error"), NA)
    if (any(bad)) stop("Analyze failed in chain ", which(bad)[1], ": ", if (is.null(runs[[which(bad)[1]]])) "the worker died (out of memory?)" else runs[[which(bad)[1]]])
    new_parts <- list(); new_last <- list()
    for (i in seq_len(chains)) {
      N <- as.matrix(runs[[i]]$nn)
      m <- regmatches(colnames(N), regexec("^NN\\.table ([0-9]+) r([0-9]+)\\.c([0-9]+)$", colnames(N)))
      stopifnot(all(lengths(m) == 4)); ix <- t(vapply(m, function(x) as.integer(x[2:4]), integer(3)))   # unit, row, col
      Di <- matrix(0, nrow(N), J * R * C); Di[, ix[, 1] + J * (ix[, 2] - 1) + J * R * (ix[, 3] - 1)] <- N; new_parts[[i]] <- Di
      nn <- matrix(0, J, R * C); nn[cbind(ix[, 1], (ix[, 2] - 1) * C + ix[, 3])] <- N[nrow(N), ]     # last latent counts, J x (R*C), column index (r-1)*C + c
      new_last[[i]] <- list(nn = nn, mu = runs[[i]]$mu)
    }
    last <- new_last; parts <- if (is.null(parts)) new_parts else Map(rbind, parts, new_parts)
    j <- bm_join_chains(parts, c(J, R, C)); draws <- j$draws; chain <- j$chain
    conv <- bm_conv(bm_chain_array(draws, chain))
    attempts <- rbind(attempts, data.frame(segment = if (is.null(attempts)) 1L else max(attempts$segment) + 1L, draws_per_chain = nrow(parts[[1]]),
                                           secs = as.numeric(Sys.time() - t0, units = "secs"), max_rhat = conv$max_rhat,
                                           min_ess_bulk = conv$min_ess_bulk, converged = conv$converged))
    if (conv$converged) break
  }
  est <- apply(draws, 2:4, mean)
  list(method = "RxCEcolInf", est = est, draws = draws, chain = chain, conv = conv,
       secs = secs_before + as.numeric(Sys.time() - t_start, units = "secs"),
       settings = list(chains = chains, keep = keep, thin = thin, num_iters = n_it, burnin_frac = burnin_frac, resume_burnin_frac = resume_burnin_frac,
                       tune_iters = tune_iters, tune_runs = tune_runs, secs_tune = secs_tune, empty_row_max = bm_check_empty_rows(est, d$rm)),
       attempts = attempts, state = list(rhos = rhos, last = last))
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
