#!/bin/bash
# bench.sh <binary> <tag> [runs=3] [threads=8]: run the pilot quantification `runs` times with the best
# configuration, report each run's wall / CPU / peak RSS and the minimum, and compare the last run's
# output with the baseline.
set -uo pipefail
B=$1; TAG=$2; RUNS=${3:-3}; T=${4:-8}
P=/scratch5/zahra/tmp/benchmark_sim/pilot/ONT_12; MAP=/scratch5/zahra/tmp/benchmark_sim/pilot/t2t.tsv
OUT=/scratch5/zahra/tmp/optim/bench_$TAG; mkdir -p $OUT; cd $OUT
best_wall=999999; best_user=999999; best_rss=999999999
for i in $(seq 1 $RUNS); do
  /usr/bin/time -v -o time_$i $B quant -i $P/af_0.95 -m $MAP -t $T -o . -r parsimony-em --lambda-size=0.0 --tau-delta=0.0 --true-umi-mode off --ed-distance-model general --ed-upper-bound 2 --dedup-edges bayes --ed-far-window 10 --model-truncation --truncation-interval-mode all --truncation-end-kernels absolute --truncation-cell-prior 0.25 --truncation-positional-power 0.25 > quant_$i.log 2>&1
  rc=$?
  wall=$(grep Elapsed time_$i | awk '{print $NF}'); user=$(grep "User time" time_$i | awk '{print $NF}'); rss=$(grep "Maximum resident" time_$i | awk '{print $NF}')
  ws=$(echo $wall | awk -F: '{ if (NF==3) print $1*3600+$2*60+$3; else print $1*60+$2 }')
  echo "$TAG run $i: exit $rc wall ${wall} (${ws}s) user ${user}s rss ${rss} kB"
  best_wall=$(echo "$ws $best_wall" | awk '{print ($1<$2)?$1:$2}'); best_user=$(echo "$user $best_user" | awk '{print ($1<$2)?$1:$2}'); best_rss=$(echo "$rss $best_rss" | awk '{print ($1<$2)?$1:$2}')
done
echo "$TAG best: wall ${best_wall}s user ${best_user}s rss ${best_rss} kB"
/nfshomes/zzare/miniconda3/envs/single-cell/bin/python -I /scratch5/zahra/tmp/optim/cmp_quants.py /scratch5/zahra/tmp/optim/baseline_pilot $OUT | tail -n 2
