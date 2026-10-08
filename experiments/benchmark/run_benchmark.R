## Run the competitor methods on one election, score them and save the results.
##
## From the package directory:
##   ELECTION <- "senc"; METHODS <- c("independence", "nslphom", "ei.MD.bayes", "RxCEcolInf"); QUICK <- TRUE
##   source("experiments/benchmark/run_benchmark.R")
## QUICK = TRUE uses short chains to check that everything runs and is scored sensibly. QUICK = FALSE uses the published settings and
## repeats ei.MD.bayes / RxCEcolInf with longer chains (up to MAX_DOUBLINGS times) until every cell count and alpha has converged.
## ssEI: run the experiment script (experiments/senc_rowwise.R, scotland_rowwise.R, nz_rowwise.R) with BM_SAVE <- TRUE and the fit is saved in
## the same folder in the same form; ELECTION <- ...; source("experiments/benchmark/bm_compare.R") then tabulates everything saved for that election.
##
## Results go to experiments/benchmark/results/<election>_<method>.rds (not tracked by git), or to $BM_RESULTS if that is set.

BM_DIR <- "experiments/benchmark"

## Command line, for batch jobs:  Rscript experiments/benchmark/run_benchmark.R <election> <methods> <quick|full> [chains] [max_doublings]
##   Rscript experiments/benchmark/run_benchmark.R scotland7 ei.MD.bayes,RxCEcolInf full 4 4
## Methods are comma separated (independence, nslphom, ei.MD.bayes, RxCEcolInf). The variables below can also be set before source().
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) >= 1) ELECTION <- .args[1]
if (length(.args) >= 2) METHODS <- strsplit(.args[2], ",")[[1]]
if (length(.args) >= 3) QUICK <- identical(.args[3], "quick")
if (length(.args) >= 4) CHAINS <- as.integer(.args[4])
if (length(.args) >= 5) MAX_DOUBLINGS <- as.integer(.args[5])
source(file.path(BM_DIR, "bm_lib.R"))
if (!exists("ELECTION")) ELECTION <- "senc"
if (!exists("METHODS")) METHODS <- c("independence", "nslphom", "ei.MD.bayes", "RxCEcolInf")
if (!exists("QUICK")) QUICK <- TRUE
if (!exists("MAX_DOUBLINGS")) MAX_DOUBLINGS <- 3
if (!exists("CHAINS")) CHAINS <- 4
if (!exists("CORES")) CORES <- bm_cores(CHAINS)

d <- bm_data(ELECTION)
cat(sprintf("%s: %d areas, %d x %d, %d voters, %d empty row margins, %d empty column margins\n", ELECTION, d$n_areas, nrow(d$kc[1, , ]),
            ncol(d$kc[1, , ]), sum(d$kc), d$n_empty_rows, d$n_empty_cols))

## settings: QUICK for a check that everything runs; otherwise the published settings with repeats until convergence.
## Override any argument of bm_eiMD / bm_rxc with EIMD_ARGS / RXC_ARGS, e.g. EIMD_ARGS <- list(thin = 20, burnin = 20000, tune_draws = 20000)
eimd_args <- if (QUICK) list(sample = 500, thin = 10, burnin = 2000, ntunes = 5, tune_draws = 5000, max_doublings = 0) else list(max_doublings = MAX_DOUBLINGS)
rxc_args  <- if (QUICK) list(keep = 500, thin = 40, tune_iters = 2000, tune_runs = 5, max_doublings = 0) else list(max_doublings = MAX_DOUBLINGS)
if (exists("EIMD_ARGS")) eimd_args <- modifyList(eimd_args, EIMD_ARGS)
if (exists("RXC_ARGS"))  rxc_args  <- modifyList(rxc_args, RXC_ARGS)

results <- list()
for (m in METHODS) {
  cat("running", m, "...\n")
  res <- switch(m,
    independence = bm_independence(d$rm, d$cm),
    nslphom = bm_nslphom(d$rm, d$cm),
    ei.MD.bayes = do.call(bm_eiMD, c(list(d$rm, d$cm, chains = CHAINS, cores = CORES), eimd_args)),
    RxCEcolInf = do.call(bm_rxc, c(list(d$rm, d$cm, chains = CHAINS, cores = CORES), rxc_args)),
    stop("unknown method ", m))
  bm_save(res, ELECTION); results[[m]] <- res
  if (!is.null(res$attempts)) print(res$attempts)
}
print(bm_report(results, d$kc))
