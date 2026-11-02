## Which space is area variation least tied to composition in? (Tier 1 check on the true tables)
##
## Question (from Kuriwaki & McCartan 2026, "coarsening at random"): every EI method identifies by assuming
## that some description of each area's table varies across areas independently of the area's composition.
## Methods differ in which description:
##   - rate space        King / Goodman / K&M: row rates q_jrc unrelated to the row mix
##   - within-row clr    ssEI's row-conditional hierarchy: log-ratio splits of each row unrelated to the row mix
##   - interaction       raking / IPF: odds ratios common to all areas (area column effects left free)
## This script measures, on the TRUE cells, how much of the real (noise-corrected) between-area variation in
## each space is explained by composition. A space whose variation composition explains least is the space
## where that method's identifying assumption is closest to holding. No model is fitted.
##
## Earlier checks this extends (project notes): ssei_seed_recovery_checks.md (rates and raking column factors
## against the row mix, raw R2 with permutation null), ssei_contextual_effects.md (within-row coordinates
## against own-row share), ssei_sorting_loglinear.md (sorting = the area x column term). New here:
##   (1) one common, noise-corrected measure across the three spaces;
##   (2) the interaction (odds ratios) itself, which the raking assumption is about, tested against both margins.
##
## Definitions, per area j, row r, column c, with a = SMOOTH added to every cell:
##   rate   q_jrc   = (n_jrc + a) / (n_jr + C a)
##   clr    l_jrc   = log(n_jrc + a) - mean_c log(n_jr. + a)            (C coordinates, rank C - 1)
##   column effect   col_jc = unweighted mean over CORE rows of l_jrc    (the area x column / sorting term)
##   interaction     int_jrc = l_jrc - col_jc                            (a function of odds ratios only)
## Core rows are the rows that are non-empty in every area; an unweighted mean over a fixed set of rows keeps
## int a function of odds ratios alone (a size-weighted mean would make it depend on composition mechanically).
## For core rows l_jr = col_j + int_jr exactly.
##
## Blocks: for each row, rate / clr / int, plus the column effect weighted by that row's voters ("col"), so that for
## core rows the clr variation splits into a column-effect part and an interaction part; and one table-level column
## effect block weighted by area totals. Core rows: senc white only (black and natam are empty in some precincts),
## so there the interaction is each row's split relative to white's; Scotland Con, Lab, LD, SNP, spoilt.
##
## Measures, for each block, traced over its C coordinates, weighted by the row's voters n_jr (the table-level
## column effect by the area total):
##   total     T  weighted between-area variance of the coordinates
##   noise     N  mean sampling variance, from a parametric bootstrap that holds the row margins fixed and
##                redraws each row multinomially from its smoothed rates
##   signal    S = T - N, the variation in the underlying propensities
##   explained E  variance of the weighted least-squares fit on composition; E0 its mean under permutation of
##                the composition across areas (what composition "explains" by chance, noise included)
##   comp_share  = (E - E0) / S, the share of the real variation that composition explains. Near 0 means the
##                space's independence assumption holds; the larger it is, the more it fails.
## Composition = ilr of the smoothed row shares ("rows"). For the interaction a second specification adds the
## column shares ("rows+cols"): odds ratios are variation-independent of both margins, so the raking assumption
## predicts no dependence on either. Rates, clr splits and column effects are NOT regressed on the column
## shares, because the column margin is built from them (the relationship would be mechanical).
##
## Usage, from the package directory:  source("experiments/structure/space_dependence.R")
## Needs ei.Datasets (Scotland, NZ) and eiPack's senc data; base R otherwise (dplyr / tidyr / purrr via data_prep.R).
## Writes experiments/structure/results/space_dependence.csv and space_dependence_<dataset>.png.

source("experiments/data_prep.R")

## ---------------------------------------------------------------------------
## Settings
## ---------------------------------------------------------------------------
## Each can be set before sourcing the script.
SMOOTH    <- get0("SMOOTH", ifnotfound = 0.5)       # added to every cell before rates / logs (as in the earlier truth checks)
MIN_ROW   <- get0("MIN_ROW", ifnotfound = 50)       # a row enters its block in area j only with at least this many voters there
MIN_AREAS <- get0("MIN_AREAS", ifnotfound = 15)     # blocks with fewer qualifying areas are skipped
N_BOOT    <- get0("N_BOOT", ifnotfound = 200)       # bootstrap replicates for the noise
N_PERM    <- get0("N_PERM", ifnotfound = 200)       # permutations for the chance-explained null
SEED      <- get0("SEED", ifnotfound = 20261102)
DATASETS  <- get0("DATASETS", ifnotfound = c("senc", "scotland7"))   # also "nz2020" or any "nzYYYY"
OUT_DIR   <- get0("OUT_DIR", ifnotfound = "experiments/structure/results")

## ---------------------------------------------------------------------------
## Data
## ---------------------------------------------------------------------------
load_senc <- function() {
  senc <- NULL
  if (requireNamespace("eiPack", quietly = TRUE)) {
    e <- new.env(); utils::data("senc", package = "eiPack", envir = e); senc <- e$senc
  } else if (file.exists(Sys.getenv("SENC_TAB"))) {                     # fallback: path to eiPack's data/senc.tab.gz
    senc <- utils::read.table(gzfile(Sys.getenv("SENC_TAB")), header = TRUE)
  } else stop("senc not found: install eiPack, or set SENC_TAB to eiPack's data/senc.tab.gz")
  rown <- c("white", "black", "natam"); coln <- c("dem", "rep", "non")
  cellnm <- outer(c("wh", "bl", "natam"), coln, paste0)
  kc <- array(0, c(nrow(senc), 3, 3), dimnames = list(NULL, rown, coln))
  for (r in 1:3) for (c in 1:3) kc[, r, c] <- senc[[cellnm[r, c]]]
  kc
}
short_names <- function(x) {
  x <- sub("Conservative and Unionist Party \\[The\\]", "Con", x); x <- sub("Labour Party \\[The\\]", "Lab", x)
  x <- sub("Liberal Democrats", "LD", x); x <- sub("Scottish National Party", "SNP", x)
  sub("Scottish Green Party", "Green", x)
}
load_dataset <- function(name) {
  if (name == "senc") return(load_senc())
  if (name == "scotland7") {
    d <- prep_scotland(7, 73); kc <- d$kc
    dimnames(kc) <- list(NULL, short_names(d$row_names), short_names(d$col_names)); return(kc)
  }
  if (grepl("^nz[0-9]{4}$", name)) {
    d <- prep_nz(as.integer(sub("nz", "", name)), 6, "candidate", n_areas = 1e6); kc <- d$kc
    dimnames(kc) <- list(NULL, d$row_names, d$col_names); return(kc)
  }
  stop("unknown dataset ", name)
}

## ---------------------------------------------------------------------------
## Coordinates
## ---------------------------------------------------------------------------
helmert <- function(K) {                                       # K x (K-1) orthonormal contrasts
  H <- matrix(0, K, K - 1)
  for (k in 1:(K - 1)) { H[1:k, k] <- 1 / k; H[k + 1, k] <- -1; H[, k] <- H[, k] / sqrt(sum(H[, k]^2)) }
  H
}
## kc: J x R x C counts. core: logical R. Returns rate, clr, int (J x R x C) and col (J x C).
coords <- function(kc, core, a = SMOOTH) {
  C <- dim(kc)[3]
  x <- kc + a
  rate <- x / as.vector(apply(x, c(1, 2), sum))                # recycled over the third index
  lg <- log(x)
  clr <- lg - as.vector(apply(lg, c(1, 2), mean))
  col <- apply(clr[, core, , drop = FALSE], c(1, 3), mean)     # J x C
  int <- clr - aperm(array(col, c(dim(col), dim(kc)[2])), c(1, 3, 2))
  list(rate = rate, clr = clr, int = int, col = col)
}
## Parametric bootstrap holding row margins fixed; returns the per-element sampling variance of each coordinate.
boot_noise <- function(kc, core, B = N_BOOT, a = SMOOTH) {
  J <- dim(kc)[1]; R <- dim(kc)[2]; C <- dim(kc)[3]
  q <- (kc + a) / as.vector(apply(kc + a, c(1, 2), sum))
  nr <- apply(kc, c(1, 2), sum)
  acc <- list(rate = 0, clr = 0, int = 0, col = 0); acc2 <- acc
  for (b in 1:B) {
    kb <- array(0, dim(kc))
    for (j in 1:J) for (r in 1:R) if (nr[j, r] > 0) kb[j, r, ] <- stats::rmultinom(1, nr[j, r], q[j, r, ])
    cb <- coords(kb, core, a)
    for (s in names(acc)) { acc[[s]] <- acc[[s]] + cb[[s]]; acc2[[s]] <- acc2[[s]] + cb[[s]]^2 }
  }
  lapply(setNames(names(acc), names(acc)), function(s) pmax(acc2[[s]] / B - (acc[[s]] / B)^2, 0) * B / (B - 1))
}

## ---------------------------------------------------------------------------
## Variance decomposition for one block
## ---------------------------------------------------------------------------
## Y: n x K coordinates, V: n x K sampling variances, w: n weights, Z: n x p composition covariates.
block_stats <- function(Y, V, w, Z, n_perm = N_PERM) {
  w <- w / sum(w)
  wmean <- colSums(w * Y)
  Tot <- sum(w * sweep(Y, 2, wmean)^2)
  Noi <- sum(w * V)
  fit_ss <- function(Zm) {
    X <- cbind(1, Zm); sw <- sqrt(w)
    f <- X %*% qr.coef(qr(sw * X), sw * Y)
    f[is.na(f)] <- 0
    sum(w * sweep(f, 2, wmean)^2)
  }
  E <- fit_ss(Z)
  E0 <- mean(replicate(n_perm, fit_ss(Z[sample(nrow(Z)), , drop = FALSE])))
  S <- Tot - Noi
  c(total = Tot, noise = Noi, signal = S, signal_share = S / Tot,
    explained = E, explained_null = E0, comp_share = (E - E0) / S, raw_R2 = E / Tot)
}

## ---------------------------------------------------------------------------
## One dataset
## ---------------------------------------------------------------------------
analyse <- function(kc, name) {
  J <- dim(kc)[1]; R <- dim(kc)[2]; C <- dim(kc)[3]
  rn <- dimnames(kc)[[2]]
  nr <- apply(kc, c(1, 2), sum); nc <- apply(kc, c(1, 3), sum); N <- rowSums(nr)
  core <- apply(nr > 0, 2, all)
  message(sprintf("%s: %d areas, %d x %d; core rows for the column effect: %s", name, J, R, C, paste(rn[core], collapse = ", ")))
  cd <- coords(kc, core)
  set.seed(SEED); nv <- boot_noise(kc, core)
  ilr_sh <- function(m) { l <- log(m + SMOOTH); (l - rowMeans(l)) %*% helmert(ncol(m)) }
  Zr <- ilr_sh(nr); Zrc <- cbind(Zr, ilr_sh(nc))
  out <- list()
  add <- function(block, space, spec, ok, Y, V, w, Z) {
    if (sum(ok) < MIN_AREAS) return(invisible())
    if (max(abs(Y[ok, ])) < 1e-12) return(invisible())          # e.g. the interaction of the only core row is 0
    s <- block_stats(Y[ok, , drop = FALSE], V[ok, , drop = FALSE], w[ok], Z[ok, , drop = FALSE])
    out[[length(out) + 1]] <<- data.frame(dataset = name, block = block, space = space, spec = spec, n_areas = sum(ok),
                                          voters = sum(w[ok]), t(s),
                                          sd_per_dim = sqrt(max(s["signal"], 0) / (C - (space != "rate"))),
                                          comp_var = s["explained"] - s["explained_null"],
                                          comp_sd_per_dim = sqrt(max(s["explained"] - s["explained_null"], 0) / (C - (space != "rate"))))
  }
  set.seed(SEED)
  for (r in 1:R) {
    ok <- nr[, r] >= MIN_ROW
    for (sp in c("rate", "clr", "int")) add(rn[r], sp, "rows", ok, cd[[sp]][, r, ], nv[[sp]][, r, ], nr[, r], Zr)
    add(rn[r], "col", "rows", ok, cd$col, nv$col, nr[, r], Zr)        # the column effect, weighted by this row's voters:
                                                                      # for core rows clr = col + int, so this splits clr
    add(rn[r], "int", "rows+cols", ok, cd$int[, r, ], nv$int[, r, ], nr[, r], Zrc)
  }
  add("column effect", "col", "rows", rep(TRUE, J), cd$col, nv$col, N, Zr)
  res <- do.call(rbind, out); rownames(res) <- NULL
  ## voter-weighted summary over rows (each row's block weighted by its voters)
  rw <- res[res$block != "column effect", ]
  summ <- do.call(rbind, lapply(split(rw, list(rw$space, rw$spec), drop = TRUE), function(d) {
    wt <- d$voters / sum(d$voters)
    data.frame(dataset = name, block = "ALL ROWS (voter-weighted)", space = d$space[1], spec = d$spec[1],
               n_areas = NA, voters = sum(d$voters), total = sum(wt * d$total), noise = sum(wt * d$noise),
               signal = sum(wt * d$signal), signal_share = sum(wt * d$signal) / sum(wt * d$total),
               explained = sum(wt * d$explained), explained_null = sum(wt * d$explained_null),
               comp_share = sum(wt * (d$explained - d$explained_null)) / sum(wt * d$signal),
               raw_R2 = sum(wt * d$explained) / sum(wt * d$total), sd_per_dim = NA,
               comp_var = sum(wt * (d$explained - d$explained_null)), comp_sd_per_dim = NA)
  }))
  list(res = rbind(res, summ), coords = cd, nr = nr, core = core)
}

## ---------------------------------------------------------------------------
## Plot: comp_share by row and space
## ---------------------------------------------------------------------------
plot_dataset <- function(res, file) {
  d <- res[res$block != "column effect" & !grepl("^ALL", res$block), ]
  d$key <- ifelse(d$spec == "rows", d$space, "int|both")
  keys <- c("rate", "clr", "col", "int", "int|both")
  rows <- unique(d$block)
  m <- matrix(NA, length(keys), length(rows), dimnames = list(keys, rows))
  for (i in seq_len(nrow(d))) m[d$key[i], d$block[i]] <- d$comp_share[i]
  grDevices::png(file, width = 1100, height = 600, res = 120)
  op <- graphics::par(mar = c(5, 5, 3, 1))
  cols <- c("#4C72B0", "#55A868", "#DD8452", "#C44E52", "#8172B2")
  graphics::barplot(pmax(m, -0.1), beside = TRUE, col = cols, ylim = c(-0.1, 1), las = 1,
                    ylab = "share of real between-area variation\nexplained by composition",
                    main = sprintf("%s: dependence on composition, by space", res$dataset[1]))
  graphics::abline(h = 0, col = "grey40")
  graphics::legend("topright", fill = cols, bty = "n", cex = 0.85,
                   legend = c("rates (King / K&M)", "within-row clr (ssEI hierarchy)", "column effect (sorting term)",
                              "interaction vs row mix", "interaction vs both margins (raking)"))
  graphics::par(op); grDevices::dev.off()
}

## ---------------------------------------------------------------------------
## Run
## ---------------------------------------------------------------------------
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
space_dep <- list()
for (nm in DATASETS) {
  kc <- load_dataset(nm)
  space_dep[[nm]] <- analyse(kc, nm)
  plot_dataset(space_dep[[nm]]$res, file.path(OUT_DIR, sprintf("space_dependence_%s.png", nm)))
}
space_dep_table <- do.call(rbind, lapply(space_dep, `[[`, "res")); rownames(space_dep_table) <- NULL
utils::write.csv(space_dep_table, file.path(OUT_DIR, "space_dependence.csv"), row.names = FALSE)

show <- space_dep_table[, c("dataset", "block", "space", "spec", "n_areas", "signal_share", "comp_share", "raw_R2",
                            "sd_per_dim", "comp_sd_per_dim")]
num <- sapply(show, is.numeric); show[num] <- lapply(show[num], round, 3)
print(show, row.names = FALSE)
