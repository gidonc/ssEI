#!/bin/bash -l
## One sweep task; submitted by experiments/hamilton/sweep.sh (which sets RMOD, R_LIBS, CODE, GRID, OUT). Not for direct use.
#SBATCH -p shared
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH --time=00:45:00
module load "$RMOD"
export OMP_NUM_THREADS=1 R_LIBS
cd "$CODE"
echo "task $SLURM_ARRAY_TASK_ID on $(hostname), $(date)"
Rscript --vanilla experiments/sweep/sweep_task.R "$SLURM_ARRAY_TASK_ID" "$OUT" "$GRID"
echo "ended $(date)"
