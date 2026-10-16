#!/bin/bash
## Slurm array job for the ssEI simulation sweep (Hamilton). Submit from the repository root:
##   Rscript experiments/sweep/make_grid.R          # writes experiments/sweep/grid.csv and prints the number of tasks
##   sbatch --array=1-540%60 experiments/sweep/sweep.slurm
## Each task fits one method to one simulated dataset (2 chains on 2 cores; most tasks take 2-10 minutes).
#SBATCH --job-name=ssei_sweep
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=sweep_logs/%A_%a.out
## Add the partition / account lines your Hamilton project needs, e.g.
## #SBATCH -p shared
## #SBATCH -A <your account>

module load r 2>/dev/null || true          # adjust to the R module you use on Hamilton
export OMP_NUM_THREADS=1
mkdir -p sweep_logs sweep_out
Rscript experiments/sweep/sweep_task.R "$SLURM_ARRAY_TASK_ID" sweep_out experiments/sweep/grid.csv
