# Simulation sweep (Hamilton)

Simulated 4 x 4 tables, 60 areas, with known truth. The grid is set in `make_grid.R`:
- association: 0.3 (census-like) or 1 (election-like);
- sorting: 0, 0.15 or 0.3 (column shifts correlated with composition);
- assimilation: 0 or 1 (small groups vote like their neighbours; half assimilated at a 5% share);
- segregation: 0.4, 0.8 or 1.6;
- 3 seeds.

That gives 108 datasets x 5 methods = 540 tasks.

## Setup (once)

```r
## R >= 4.2 with rstan; from the repository root (branch `contextual`):
install.packages(c("remotes", "quadprog", "dplyr"))
remotes::install_github("cran/eiPack")            # King's ei.MD.bayes
## install the package from this checkout (compiles the Stan model: about 5-10 minutes, 4+ GB of memory)
remotes::install_local(".", force = TRUE)
```

## Run

```bash
git checkout contextual && git pull
Rscript experiments/sweep/make_grid.R
sbatch --array=1-540%60 experiments/sweep/sweep.slurm     # %60 = at most 60 tasks at once
## when done (or part-way):
Rscript experiments/sweep/collect.R sweep_out
```

Each task writes one `.rds` to `sweep_out/` and skips work whose file already exists, so resubmitting the array resumes. `collect.R` writes `sweep_out/results_long.csv` (one row per method, dataset and seed: area error, national error, coverage, Rhat, minutes) and prints the seed-averaged table.

Quick local test of one task: `Rscript experiments/sweep/sweep_task.R 1 sweep_test`.

## Rank-1 against full-rank sorting (grid_rank.csv)

`make_grid_rank.R` writes `grid_rank.csv` (198 tasks):
- tasks 1-90: reference methods (seconds each). These include Goodman with a contextual sorting term common to all rows, fitted three ways:
  - all composition directions (`goodman_ctx`);
  - one estimated direction (`goodman_rr1`);
  - the true direction (`goodman_orc`, an oracle).
- tasks 91-198: ssei_est, ssei_kingqi and King on the rank-1 datasets.

The full-rank King and ssEI fits are already in `results/results_long.csv`: same seeds, identical data. Rank-1 datasets share compositions, area sizes and the sorting draw with their full-rank twins.

```bash
sbatch --array=1-198%60 <your working Slurm script, with the grid and output directory below>
#   Rscript --vanilla experiments/sweep/sweep_task.R $SLURM_ARRAY_TASK_ID /nobackup/<user>/ssei_rank_out experiments/sweep/grid_rank.csv
Rscript --vanilla experiments/sweep/collect.R /nobackup/<user>/ssei_rank_out
```
