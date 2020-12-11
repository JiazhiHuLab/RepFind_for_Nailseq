#!/bin/bash
# 20190320 :update: use >10k HUEdU peak
# sysinfo_page - A script to call peak from HUEdU-seq according to Andre
case=$1
ctrl=$2
name=${3:-huseq}
bwfile=${4:-bwfile} # We compute mean value, so the bwfile must be evenly windowed 

fc_cutoff=400
maxgap=10000 # K562
#maxgap=20000 # GM12878
minsize=10000

#--- input and output ---#
# input bams of case and ctrl
# output peaks, initiation zones(merged in 20kb)
# compute the reads HU-seq's RPKM for each regions 
# Optional output: can output the mean value of a bigwig files for each 

# call peak and remove chrY and chrM peaks
# -slocal 5000 -llocal 50000
source activate py26
macs14 -t $case -c $ctrl -f BAM -g hs -p 1e-5   \
    --nolambda --nomodel --keep-dup=all -n $name 

#column 8 is the fold_enrichment

grep -v "#" ${name}_peaks.xls | awk -v x=$fc_cutoff 'NR>2 {if($8>=x){print $1,$2,$3,$8}}' OFS="\t" | grep -v "chrY" | grep -v "chrM" > ${name}_peaks.bed


# remove peaks overlaped with blacklist
blacklist=/home/workdir/chen/databases/genomes/hg19/blacklists/DAC_BlackListed_regions_ENCFF000KJP.bed
bedtools intersect -v -a ${name}_peaks.bed -b $blacklist > ${name}_peaks_rmBlack.bed


# identify initiation zone 
bedtools merge -i ${name}_peaks_rmBlack.bed -d $maxgap  |  awk -v x=$minsize '($3-$2)>=x'  > ${name}_initzone.bed 

## ------- compute HUEdU RPKM in initiation zones ----------#
x=`samtools view -c $case`
cut -f 1,2,3 ${name}_initzone.bed | bedtools intersect -a stdin -b $case -c -wa | awk -vx=$x '{print $1,$2,$3, $4*1e+09/(x*($3-$2)) }' OFS='\t' >  ${name}_initzone_rpkm.bdg


## --------- compute E-B signal in HUEdU initiation zones ---------------------#
if ["$4"!=""]; then
    bwtool summary ${name}_initzone.bed $bwfile stdout | awk '{print $1,$2,$3,$8}' OFS='\t' > ${name}_${bwfile%.bw}_mean.bdg
fi

conda deactivate py26
