## Scottish Parliament 2007: constituency vote (rows) by regional list vote (columns), soft multinomial, row-by-row margin model.
##
## SIZE 3 (Labour, SNP, Other), 5 (Con, Lab, LD, SNP, Other) or 7 (all seven categories); N_AREAS 30 (the fixed random subset,
## set.seed(1234)) or 73 (every constituency). Rows and columns are the same parties, so the reference column of each row is its own
## (loyal) column, the loyalty split of each row has a prior centred on 75%, and the sigmas are shared in two tiers:
## loyalty and congruent-vs-cross splits, and all the rest.
##
## Usage, from the package directory, after installing ssEI and ei.Datasets:  source("experiments/scotland_rowwise.R")
## Set ALLOC <- 2L for the odds-ratio logit allocation instead of adjusted row (3).

source("experiments/rowwise_helpers.R")
source("experiments/data_prep.R")

SIZE <- 5; N_AREAS <- 30; ALLOC <- 3L
SEED <- 1234; CHAINS <- 4; ITER <- 1000; WARMUP <- 500
E_SD_SCALE_SMALL <- 1; E_SD_SCALE_OFFDIAG <- 1;
## E_SD_SCALE_SMALL < 1 tightens the E_rc prior sd on coordinates that involve only small columns

d <- prep_scotland(SIZE, N_AREAS)
R <- length(d$row_names); C <- length(d$col_names)
kc <- d$kc; dimnames(kc) <- list(NULL, d$row_names, d$col_names)
cat(sprintf("Scotland 2007: %d areas, %d x %d, %d voters\n", dim(kc)[1], R, C, sum(kc)))
cat("rows with no votes in some area:", sum(apply(kc, c(1, 2), sum) == 0), "\n")

V_ilr <- rw_row_priority_basis(d$row_names, d$col_names, d$bloc)

fit <- rw_fit(kc, V_ilr, alloc = ALLOC,
              rem_col = seq_len(R),                        # each row's own column
              tiers = c(1, 1, rep(2, C - 3)),              # sigma tiers within each row
              loyalty_mean = qlogis(0.75),
              E_sd_scale_small = E_SD_SCALE_SMALL, E_sd_scale_offdiag = E_SD_SCALE_OFFDIAG,
              chains = CHAINS, iter = ITER, warmup = WARMUP, seed = SEED)
