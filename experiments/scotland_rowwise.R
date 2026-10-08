## Scottish Parliament 2007: constituency vote (rows) by regional list vote (columns), soft multinomial, row-by-row margin model.
##
## SIZE 3 (Labour, SNP, Other), 5 (Con, Lab, LD, SNP, Other) or 7 (all seven categories); N_AREAS 30 (the fixed random subset,
## set.seed(1234)) or 73 (every constituency). The default is 7 x 7 on all 73: the case used to compare with other methods. It has
## many empty rows (constituencies with no Green and / or no Other candidate). For a quick run use SIZE <- 5; N_AREAS <- 30.
## Rows and columns are the same categories, so the reference column of each row is its own column. Only the four major parties
## (Con, Lab, LD, SNP) have a loyalty split, with a prior centred on 75%; the spoilt, Other and Green rows have no loyalty split (their
## tree is bisections of all the columns, largest first), and there are no blocs. The sigmas are shared in two tiers: the first two
## splits of each row, and all the rest.
## For the earlier set-up (loyalty in every row, blocs) set LOYAL_ROWS <- NULL; BLOCS <- TRUE before the fit.
##
## Usage, from the package directory, after installing ssEI and ei.Datasets:  source("experiments/scotland_rowwise.R")
## Set ALLOC <- 2L for the odds-ratio logit allocation instead of adjusted row (3).

source("experiments/rowwise_helpers.R")
source("experiments/data_prep.R")

SIZE <- 7; N_AREAS <- 73; ALLOC <- 3L                 # the benchmark case: all seven categories, every constituency
SEED <- 1234; CHAINS <- 4; ITER <- 1000; WARMUP <- 500
E_SD_SCALE_SMALL <- 1; E_SD_SCALE_OFFDIAG <- 1;
## E_SD_SCALE_SMALL / E_SD_SCALE_OFFDIAG < 1 tighten the E_rc prior sd (small-column coordinates / all but the loyalty splits)
E_CENTRE <- "uniform"; BLOC_AFFINITY <- 0; E_CENTRE_SHIFT <- 0;   # E_CENTRE "qi": non-loyalty coordinates centred on quasi-independence (see rw_qi_centre)
BLOCS <- FALSE;   # FALSE: no blocs, the tree is loyalty then bisections of the other columns, largest columns first (the general default)
LOYAL_ROWS <- c("Conservative and Unionist Party [The]", "Labour Party [The]", "Liberal Democrats", "Scottish National Party")
## rows with a loyalty split (NULL: every row). At SIZE 3 and 5 only the ones present are used.

d <- prep_scotland(SIZE, N_AREAS)
R <- length(d$row_names); C <- length(d$col_names)
kc <- d$kc; dimnames(kc) <- list(NULL, d$row_names, d$col_names)
cat(sprintf("Scotland 2007: %d areas, %d x %d, %d voters\n", dim(kc)[1], R, C, sum(kc)))
cat("rows with no votes in some area:", sum(apply(kc, c(1, 2), sum) == 0), "\n")

loyal_rows <- if (is.null(LOYAL_ROWS)) NULL else intersect(LOYAL_ROWS, d$row_names)
V_ilr <- rw_row_priority_basis(d$row_names, d$col_names, if (BLOCS) d$bloc, bisect_order = if (!BLOCS) colSums(d$cm), loyal_rows = loyal_rows)

fit <- rw_fit(kc, V_ilr, alloc = ALLOC,
              rem_col = seq_len(R),                        # each row's own column
              tiers = c(1, 1, rep(2, C - 3)),              # sigma tiers within each row
              loyalty_mean = qlogis(0.75), loyal_rows = loyal_rows,
              E_sd_scale_small = E_SD_SCALE_SMALL, E_sd_scale_offdiag = E_SD_SCALE_OFFDIAG,
              E_centre = E_CENTRE, bloc = if (BLOCS) d$bloc, bloc_affinity = BLOC_AFFINITY, E_centre_shift = E_CENTRE_SHIFT,
              chains = CHAINS, iter = ITER, warmup = WARMUP, seed = SEED)

## Set BM_SAVE <- TRUE before sourcing to save this fit in the form used to compare with other methods (experiments/benchmark)
if (isTRUE(get0("BM_SAVE"))) {
  if (SIZE == 7 && N_AREAS >= 73) { source("experiments/benchmark/bm_lib.R"); bm_ssei_save(fit, kc, "scotland7", label = sprintf("ssEI alloc %d%s%s", ALLOC, if (is.null(LOYAL_ROWS)) "" else " major-loyal", if (BLOCS) "" else " no-blocs")) }
  else message("not saved: the benchmark case is SIZE 7 with all 73 constituencies")
}
