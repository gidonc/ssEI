# Running sweeps on Hamilton

No edits on Hamilton are needed. Each step is one command.

## Each time

1. **RStudio Server (Open OnDemand), in the ssEI project.** Pull: Git pane > Pull, or `git pull` in the RStudio terminal.
2. **Only if the package itself changed** (R/, inst/stan/, DESCRIPTION). From the shell:
   ```bash
   bash experiments/hamilton/sweep.sh install
   ```
   This runs as a Slurm job. It reinstalls ssEI into `/nobackup/<user>/R-libs` when the package source differs from what is installed, and says "up to date" otherwise. Compiling takes about 5-10 minutes, and the compiler output goes to `<library>/ssEI_install.log`.

   In RStudio, `source("experiments/hamilton/install.R")` after Session > Restart R does the same, with the output kept out of the console.
3. **Hamilton shell** (OnDemand shell, or `ssh` and `cd ~/R/projects/ssEI`):
   ```bash
   bash experiments/hamilton/sweep.sh submit base       # grid experiments/sweep/grid_base.csv
   bash experiments/hamilton/sweep.sh status base       # progress, failed logs, missing tasks
   bash experiments/hamilton/sweep.sh rerun base        # only the tasks with no result (e.g. after time-outs)
   bash experiments/hamilton/sweep.sh collect base --push   # results_long.csv -> experiments/sweep/results/base_long.csv, committed and pushed
   ```
   Without `--push`, collect copies the file into the checkout and you can commit it from the RStudio Git pane.

## What the script does

- **R version.** It loads the R module the installed ssEI was built with. It refuses to submit if ssEI does not load, or was installed from different package source than the checkout. It then tells you to run `install.R`.
- **Code snapshot.** Tasks run from a copy of `experiments/` taken at submission, so pulling new code during a sweep does not affect it.
- **Output locations.**
  - results: `/nobackup/<user>/ssei_NAME_out`
  - logs: `/nobackup/<user>/slurmout/ssei_NAME/`
  - completed tasks are skipped on resubmission.
- **Slurm defaults.** Partition `shared`, 2 cores, 4G, 45 minutes, at most 60 tasks at once. Override them on the command line:
  ```bash
  bash experiments/hamilton/sweep.sh submit base --time=01:30:00 --mem=8G
  ```
- **Overrides by environment variable.**
  - `SSEI_LIB`: library path.
  - `SSEI_R_MODULE`: R module name, if the automatic choice fails.

## Grids

| NAME | File | Contents |
|---|---|---|
| main | grid.csv | 540 tasks: the first sweep |
| rank | grid_rank.csv | 198 tasks: rank-1 sorting |
| form | grid_form.csv | 72 reference tasks: sorting linear in shares |
| base | grid_base.csv | 216 tasks: low-dimensional base (`ssei_fixed`) and kappa-scaled (`ssei_scaled`) on the main datasets |

To make a new grid, add `make_grid_NAME.R` writing `grid_NAME.csv` (columns as in the others, plus `task`). Commit both, pull on Hamilton, then submit NAME.
