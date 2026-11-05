#!/bin/bash
# verify_all.sh: run the original and the optimized binary on every data set of the per-cell report, in
# both deduplication configurations, and compare the quantification matrices.
#   config bayes : --dedup-edges bayes --ed-far-window 10 (the requested best configuration)
#   config rule  : the report's PC13 flags (general ed2, rule edges)
# usage: verify_all.sh <out root>
set -uo pipefail
OUT=$1; mkdir -p $OUT
BASE=/scratch5/zahra/tmp/optim/bin/alevin-fry-base
NEW=/scratch5/zahra/tmp/optim/bin/alevin-fry-step10
SIMMAP=/fs/cbcb-lab/rob/students/zahra/nanosim_data/reference/new_files/LRGASP.t2t.from_transcriptome
EXPMAP=/fs/cbcb-lab/rob/students/zahra/sg-nex-data/Reference_sequin/LRGASP.t2t.from_transcriptome
TR="-r parsimony-em --lambda-size=0.0 --tau-delta=0.0 --true-umi-mode off --ed-distance-model general --ed-upper-bound 2 --model-truncation --truncation-interval-mode all --truncation-end-kernels absolute --truncation-cell-prior 0.25 --truncation-positional-power 0.25"
declare -A IN MAP
IN[sim_new_ONT]=/scratch1/zahra/single-cell-quant/new_results/sim_new_data_exp_end_all_models/0.95/ONT/forward/alevin-fry-likelihood;      MAP[sim_new_ONT]=$SIMMAP
IN[sim_new_PacBio]=/scratch1/zahra/single-cell-quant/new_results/sim_new_data_exp_end_all_models/0.95/PacBio/forward/alevin-fry-likelihood; MAP[sim_new_PacBio]=$SIMMAP
IN[sim_old_ONT]=/scratch1/zahra/single-cell-quant/new_results/sim_old_data_exp_end_all_models/0.95/ONT/forward/alevin-fry-likelihood;      MAP[sim_old_ONT]=$SIMMAP
IN[sim_old_PacBio]=/scratch1/zahra/single-cell-quant/new_results/sim_old_data_exp_end_all_models/0.95/PacBio/forward/alevin-fry-likelihood; MAP[sim_old_PacBio]=$SIMMAP
IN[exp_ONT]=/scratch5/zahra/single-cell-quant/new_results/exp_all_end_all_models/0.95/ONT/forward/alevin-fry-likelihood;                   MAP[exp_ONT]=$EXPMAP
IN[exp_PacBio]=/scratch5/zahra/single-cell-quant/new_results/exp_all_end_all_models/0.95/PacBio/forward/alevin-fry-likelihood;             MAP[exp_PacBio]=$EXPMAP

run_dataset() {
  ds=$1; in=${IN[$ds]}; map=${MAP[$ds]}
  for cfg in bayes rule; do
    extra=""; [ $cfg = bayes ] && extra="--dedup-edges bayes --ed-far-window 10"
    for which in base new; do
      b=$BASE; [ $which = new ] && b=$NEW
      d=$OUT/$ds/${cfg}_$which; mkdir -p $d; cd $d
      /usr/bin/time -v -o quant_time $b quant -i $in -m $map -t 8 -o . $TR $extra > quant.log 2>&1
      echo "$ds $cfg $which exit $? $(grep Elapsed quant_time | awk '{print $NF}') rss $(grep 'Maximum resident' quant_time | awk '{print $NF}') kB" >> $OUT/runs.log
    done
    res=$(/nfshomes/zzare/miniconda3/envs/single-cell/bin/python -I /scratch5/zahra/tmp/optim/cmp_quants.py $OUT/$ds/${cfg}_base $OUT/$ds/${cfg}_new 2>&1 | tail -n 3 | tr '\n' ' ')
    echo "$ds $cfg: $res" >> $OUT/compare.log
  done
}
export -f run_dataset; export OUT BASE NEW TR; export -A IN MAP 2>/dev/null
# bash cannot export associative arrays to xargs children: run the chains as background jobs in this shell
for ds in sim_new_ONT sim_new_PacBio sim_old_ONT sim_old_PacBio exp_ONT exp_PacBio; do run_dataset $ds & done
wait
echo "===== runs ====="; cat $OUT/runs.log; echo "===== comparisons ====="; cat $OUT/compare.log
