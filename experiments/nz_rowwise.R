## New Zealand: candidate vote by party vote, by electorate, soft multinomial, row-by-row margin model.
##
## YEAR 2002, 2005, ..., 2020. SIZE 6 (Labour, National, Green, NZFirst, ACT, Other) or 4 (Labour, National, Green, Other; only
## electorates with a candidate in every category). ROWS "candidate" puts the candidate vote on the rows, "party" the party vote.
## N_AREAS electorates are drawn at random (set.seed(1234)); use a large number for all of them.
## Rows and columns are the same parties, so the reference column of each row is its own (loyal) column, the loyalty split of each
## row has a prior centred on 75%, and the sigmas are shared in two tiers: loyalty and congruent-vs-cross splits, and all the rest.
##
## Usage, from the package directory, after installing ssEI and ei.Datasets:  source("experiments/nz_rowwise.R")
## Set ALLOC <- 2L for the odds-ratio logit allocation instead of adjusted row (3).

source("experiments/rowwise_helpers.R")
source("experiments/data_prep.R")

YEAR <- 2020; SIZE <- 6; ROWS <- "candidate"; N_AREAS <- 30; ALLOC <- 3L
SEED <- 1234; CHAINS <- 4; ITER <- 1000; WARMUP <- 500
E_SD_SCALE_SMALL <- 1; E_SD_SCALE_OFFDIAG <- 1;
## E_SD_SCALE_SMALL / E_SD_SCALE_OFFDIAG < 1 tighten the E_rc prior sd (small-column coordinates / all but the loyalty splits)
E_CENTRE <- "uniform"; BLOC_AFFINITY <- 0; E_CENTRE_SHIFT <- 0;   # E_CENTRE "qi": non-loyalty coordinates centred on quasi-independence (see rw_qi_centre)
BLOCS <- TRUE;   # FALSE: no blocs, the tree is loyalty then bisections of the other columns, largest columns first (the general default)

d <- prep_nz(YEAR, SIZE, ROWS, N_AREAS)
R <- length(d$row_names); C <- length(d$col_names)
kc <- d$kc; dimnames(kc) <- list(NULL, d$row_names, d$col_names)
cat(sprintf("New Zealand %d: %d electorates, %d x %d, %d voters\n", YEAR, dim(kc)[1], R, C, sum(kc)))
cat("empty rows (areas with no votes for a category):", sum(apply(kc, c(1, 2), sum) == 0), "of", dim(kc)[1] * R, "\n")

V_ilr <- rw_row_priority_basis(d$row_names, d$col_names, if (BLOCS) d$bloc, bisect_order = if (!BLOCS) colSums(d$cm))

fit <- rw_fit(kc, V_ilr, alloc = ALLOC,
              rem_col = seq_len(R),                        # each row's own column
              tiers = c(1, 1, rep(2, C - 3)),              # sigma tiers within each row
              loyalty_mean = qlogis(0.75),
              E_sd_scale_small = E_SD_SCALE_SMALL, E_sd_scale_offdiag = E_SD_SCALE_OFFDIAG,
              E_centre = E_CENTRE, bloc = if (BLOCS) d$bloc, bloc_affinity = BLOC_AFFINITY, E_centre_shift = E_CENTRE_SHIFT,
              chains = CHAINS, iter = ITER, warmup = WARMUP, seed = SEED)

## Set BM_SAVE <- TRUE before sourcing to save this fit in the form used to compare with other methods (experiments/benchmark)
if (isTRUE(get0("BM_SAVE"))) {
  source("experiments/benchmark/bm_lib.R"); nm <- sprintf("nz%d_%d", YEAR, SIZE)
  if (ROWS == "candidate" && dim(kc)[1] == bm_data(nm)$n_areas) bm_ssei_save(fit, kc, nm, label = sprintf("ssEI alloc %d", ALLOC))
  else message("not saved: the benchmark uses every electorate with the candidate vote on the rows")
}
