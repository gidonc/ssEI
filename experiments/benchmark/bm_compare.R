## Tabulate every saved result for one election, the true cells being rebuilt from the data.
##   ELECTION <- "scotland7"; source("experiments/benchmark/bm_compare.R")
BM_DIR <- "experiments/benchmark"
source(file.path(BM_DIR, "bm_lib.R"))
if (!exists("ELECTION")) stop("set ELECTION first")
d <- bm_data(ELECTION)
tab <- bm_report(bm_load(ELECTION), d$kc)
print(tab[order(tab$ei_area), ])
