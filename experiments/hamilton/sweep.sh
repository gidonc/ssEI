#!/bin/bash
## Simulation sweeps on Hamilton, one command per step. Run from a Hamilton shell (any directory):
##   bash experiments/hamilton/sweep.sh submit  NAME [extra sbatch options, e.g. --time=01:30:00]
##   bash experiments/hamilton/sweep.sh status  NAME
##   bash experiments/hamilton/sweep.sh rerun   NAME [extra sbatch options]   (submits only the tasks with no result yet)
##   bash experiments/hamilton/sweep.sh collect NAME [--push]
## NAME picks the grid experiments/sweep/grid_NAME.csv ("main" = grid.csv). Results: /nobackup/$USER/ssei_NAME_out;
## logs: /nobackup/$USER/slurmout/ssei_NAME. Tasks run from a snapshot of experiments/ taken at submission, so pulling new code
## while a sweep runs does not affect it. The R module is the one the installed ssEI was built with (override: SSEI_R_MODULE).
set -eo pipefail
REPO=$(cd "$(dirname "$0")/../.." && pwd)
CMD=${1:?command: submit, status, rerun or collect}; NAME=${2:?name of the grid, e.g. base}; shift 2
LIB=${SSEI_LIB:-/nobackup/$USER/R-libs}
GRID=$REPO/experiments/sweep/grid_$NAME.csv; [ "$NAME" = main ] && GRID=$REPO/experiments/sweep/grid.csv
ROOT=${SSEI_OUT_ROOT:-/nobackup/$USER}
OUT=$ROOT/ssei_${NAME}_out; LOGS=$ROOT/slurmout/ssei_$NAME; CODE=$OUT/code
[ -f "$GRID" ] || { echo "no grid file $GRID"; exit 1; }
N=$(( $(wc -l < "$GRID") - 1 ))
type module > /dev/null 2>&1 || source /etc/profile > /dev/null 2>&1 || true
load_r() {
  BUILT=$(sed -n 's/^Built: R \([0-9.]*\);.*/\1/p' "$LIB/ssEI/DESCRIPTION" 2> /dev/null)
  RMOD=${SSEI_R_MODULE:-r/$BUILT}
  [ -n "$BUILT" ] || [ -n "$SSEI_R_MODULE" ] || { echo "ssEI is not installed in $LIB: in RStudio run source(\"experiments/hamilton/install.R\")"; exit 1; }
  module load "$RMOD" > /dev/null 2>&1 || { echo "could not load module $RMOD (set SSEI_R_MODULE=<module name>)"; exit 1; }
  export R_LIBS=$LIB
}
missing() {   # tasks with no result, as Slurm ranges: 3-7,9,
  local s="" e="" i
  for i in $(seq 1 "$N"); do
    if compgen -G "$OUT/$(printf 't%04d' "$i")_*.rds" > /dev/null; then [ -z "$s" ] || { [ "$s" = "$e" ] && echo -n "$s," || echo -n "$s-$e,"; s=""; }
    else [ -z "$s" ] && s=$i; e=$i; fi
  done
  [ -z "$s" ] || { [ "$s" = "$e" ] && echo -n "$s," || echo -n "$s-$e,"; }; }
submit() {   # $1 = array spec; rest = sbatch options
  local arr=$1; shift
  load_r; Rscript --vanilla "$REPO/experiments/hamilton/check.R" "$REPO" || exit 1
  mkdir -p "$OUT" "$LOGS"; rm -rf "$CODE"; mkdir -p "$CODE"; cp -r "$REPO/experiments" "$CODE/"
  sbatch --array="$arr" -J "ssei_$NAME" -o "$LOGS/%A_%a.out" \
    --export=ALL,RMOD="$RMOD",R_LIBS="$LIB",CODE="$CODE",OUT="$OUT",GRID="$CODE/experiments/sweep/$(basename "$GRID")" "$@" \
    "$REPO/experiments/hamilton/task.slurm"
  echo "results -> $OUT ; logs -> $LOGS ; check progress: bash experiments/hamilton/sweep.sh status $NAME"
}
case $CMD in
  submit) submit "1-$N%60" "$@" ;;
  rerun)  M=$(missing); [ -n "$M" ] || { echo "all $N tasks have results"; exit 0; }; echo "rerunning tasks ${M%,}"; submit "${M%,}%60" "$@" ;;
  status)
    DONE=$(compgen -G "$OUT/t*.rds" | wc -l || true); Q=$(squeue -u "$USER" -n "ssei_$NAME" -h 2> /dev/null | wc -l)
    echo "$NAME: $DONE of $N tasks have results; $Q array tasks still queued or running"
    BAD=$(grep -l -E "Error|Execution halted|DUE TO TIME LIMIT|CANCELLED|oom-kill" "$LOGS"/*.out 2> /dev/null | head -20 || true)
    [ -z "$BAD" ] || { echo "logs reporting problems (first 20):"; echo "$BAD"; echo "first problem:"; grep -h -m 3 -E "Error|halted|LIMIT|CANCELLED|oom" $(echo "$BAD" | head -1); }
    if [ "$Q" -eq 0 ] && [ "$DONE" -lt "$N" ]; then M=$(missing); echo "nothing queued and tasks missing: ${M%,}  -> bash experiments/hamilton/sweep.sh rerun $NAME"; fi ;;
  collect)
    load_r; Rscript --vanilla "$REPO/experiments/sweep/collect.R" "$OUT" | tail -3
    mkdir -p "$REPO/experiments/sweep/results"; cp "$OUT/results_long.csv" "$REPO/experiments/sweep/results/${NAME}_long.csv"
    echo "copied to experiments/sweep/results/${NAME}_long.csv"
    if [ "${1:-}" = "--push" ]; then cd "$REPO"; git add "experiments/sweep/results/${NAME}_long.csv"
      git commit -q -m "Sweep results: $NAME ($(date +%F))" && git push -q && echo "pushed to $(git rev-parse --abbrev-ref HEAD)"; fi ;;
  *) echo "unknown command $CMD"; exit 1 ;;
esac
