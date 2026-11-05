#!/bin/bash
# End-to-end test of a binary on the pilot BAM: convert -> generate-permit-list -> collate -> quant,
# then compare the quantification with the baseline (names-in-RAD files, scratch binary).
# usage: pipeline_test.sh <binary> <out dir>
set -euo pipefail
B=$1; OUT=$2; mkdir -p $OUT
P=/scratch5/zahra/tmp/benchmark_sim/pilot/ONT_12
WL=/scratch5/zahra/tmp/benchmark_sim/pilot/3M-february-2018.rc.txt
MAP=/scratch5/zahra/tmp/benchmark_sim/pilot/t2t.tsv
d=$OUT/af_0.95; mkdir -p $d
/usr/bin/time -v -o $d/convert_time $B convert --threads 8 --filter_best --bam $P/T.tagged.bam --output $d/map.rad --score-threshold 0.95 > $d/convert.log 2>&1
/usr/bin/time -v -o $d/permit_time $B generate-permit-list --input $d --expected-ori fw --output-dir $d --unfiltered-pl $WL > $d/permit.log 2>&1
touch $d/unmapped_bc_count.bin
/usr/bin/time -v -o $d/collate_time $B collate -i $d -r $d -t 8 > $d/collate.log 2>&1
q=$OUT/quant; mkdir -p $q; cd $q
/usr/bin/time -v -o quant_time $B quant -i $d -m $MAP -t 8 -o . -r parsimony-em --lambda-size=0.0 --tau-delta=0.0 --true-umi-mode off --ed-distance-model general --ed-upper-bound 2 --dedup-edges bayes --ed-far-window 10 --model-truncation --truncation-interval-mode all --truncation-end-kernels absolute --truncation-cell-prior 0.25 --truncation-positional-power 0.25 > quant.log 2>&1
for f in convert permit collate; do printf "%-8s %s  RSS %s kB\n" $f "$(grep Elapsed $d/${f}_time | awk '{print $NF}')" "$(grep 'Maximum resident' $d/${f}_time | awk '{print $NF}')"; done
printf "%-8s %s  RSS %s kB\n" quant "$(grep Elapsed quant_time | awk '{print $NF}')" "$(grep 'Maximum resident' quant_time | awk '{print $NF}')"
ls -la $d/map.rad $d/map.collated.rad $P/af_0.95/map.rad $P/af_0.95/map.collated.rad | awk '{print $5, $9}'
/nfshomes/zzare/miniconda3/envs/single-cell/bin/python -I /scratch5/zahra/tmp/optim/cmp_quants.py /scratch5/zahra/tmp/optim/baseline_pilot $q
