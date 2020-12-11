#!/bin/bash
# for call peak for EdU/BrdU

echo -e "Help: script case_bam ctrl_bam name "
case=$1
#ctrl=/home/workdir/chen/proj/work/LY_related/DRC-seq/HU022/data/LY072a/LY072a_ass_trimmed_dedup.bam
#name=$2
ctrl=$2
name=$3

source activate py2
SICER.py -t $case  -c $ctrl --mapq 30 -rt 0 -w 5000 -g 3 -fs 150 -gs 0.793 > ${name}_peak.bed
conda deactivate py2
