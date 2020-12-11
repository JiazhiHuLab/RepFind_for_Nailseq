#!/bin/bash
# Name: Chen Ai
# Date: 20190402 update for extending EdUrich
# Used for: process EdU and BrdU labeled sequencing reads, Nail-seq
# Input: 
#       1. table with edu and brdu bams path
#       2. window size and slide size
#       3. name prefix
#       4. default save files in current directory       
# Output: 
#     1. E-B signals: normalized batch and smoothed
#     2. EdU, BrdU normalized library size, RPM
#     3. E-B local maxima peak

set -eu -o pipefail

SCRIPT=`basename ${BASH_SOURCE[0]}`

# Help function
function HELP {
    echo -e "Program: RepFind Nail-seq EdU-BrdU peak"
    echo -e "Version: 1.0"
    echo -e "Author: Chen Ai"
    echo -e "\nUsage: $SCRIPT [options] "
    echo -e "\nOptions:  "
    echo -e "-1 --Input EdU bamfile, REQUIRED"
    echo -e "-2 --Input BrdU bamfile, REQUIRED"
    echo -e "-g  --Genome [hg19/mm9/mm10], REQUIRED"
    echo -e "-n  --Output prefix of EdU-BrdU, name_EdU-BrdU.bdg . REQUIRED"
    echo -e "-w --Window [INT]bp default 5000 bp."
    echo -e "-b --Broadpeak of output Nailseq_call_broadPeak.sh "
    echo -e "-h  --Displays this help message.\n"
    echo -e "Example: $SCRIPT -1 EdU.bam -2 BrdU.bam -g hg19 -o EdU-BrdU.peak "
    exit 1
}

NUMARGS=$#
if [ $NUMARGS -eq 0 ]; then
    echo -e \\n"Error: Number of arguments: $NUMARGS\n"
    HELP
fi

# Check the options
while getopts "1:2:g:n:o:p:r:l:h"  arg
do
    case $arg in 
        1) 
            edubam=$OPTARG 
            echo -e "EdU bam: $bam1"
            ;;
        2) 
            brdubam=$OPTARG
            echo -e "BrdU bam: $bam2"
            ;;
        g) 
            genome=$OPTARG
            echo -e "Genome is $genome"
            ;;
        n)  name=$OPTARG
            ;;
        w)  win=$OPTARG
            ;;
        b)  broadpeak=$OPTARG
            ;;
        h)  HELP 
            ;;
        \?) 
            echo "Unknown Arugument"
            HELP
            exit 2
            ;;
    esac
done

case $genome in 
    hg19) 
        bw2_idx=/home/ubuntu/genomes/bowtie2_indexes/hg19
        bwa_ref=/home/mengzhu/database/bwa_indexes/hg19/hg19
        chrom_ref=/home/ubuntu/genomes/hg19/annotation/ChromInfo.txt
        blacklist=/home/workdir/chen/databases/genomes/blacklists/DAC_BlackListed_regions_ENCFF000KJP.bed
        ;;
    mm10)
        bw2_idx=/home/ubuntu/genomes/bowtie2_indexes/mm10
        bwa_ref=/home/mengzhu/database/bwa_indexes/mm10/mm10
        chrom_ref=/home/ubuntu/genomes/mm10/annotation/ChromInfo.txt
        blacklist=/home/workdir/chen/databases/genomes/blacklists/mm10.blacklist.bed.gz
        ;;
    mm9)
        # mouse genome : mm9
        bw2_idx=/home/ubuntu/genomes/bowtie2_indexes/mm9
        bwa_ref=/home/workdir/chen/databases/bwa_indexes/mm9.fa
        chrom_ref=/home/ubuntu/genomes/mm9/annotation/ChromInfo.txt
        blacklist=/home/workdir/chen/databases/genomes/blacklists/mm9.blacklist.bed.gz
        ;;
    *)
        echo -e "Error: Unknown genome assembly"
        HELP
        exit 2
        ;;
esac   

mywindow=./chr_win${win}_rmblack.bed
script_dir=`dirname ${BASH_SOURCE[0]}`
dir=$script_dir/src 


#---------------------------- make windows  ----------------------------------#

#if [[ ! -f $mywindow ]]; then
bedtools makewindows -w $win -g $chrom_ref | bedtools intersect -a stdin -b $blacklist -v -wa | sort -k1,1 -k2,2n > $mywindow
#fi


#--------------------------- get RPM bdg for edu and brdu --------------------#
# prepare input bams path
# input_table=gm78_eb_071.txt 
#bams=`cut -f 3 $input_table`
#edubam=$(echo $bams| cut -f1 -d " ")
#brdubam=$(echo $bams| cut -f2 -d " ")

if [[ ! -f $edubam ]]; then
    echo " Error: $edubam : is not found"
fi

if [[ ! -f $brdubam ]]; then
    echo " Error: $brdubam : is not found"
fi


# convert to RPM bedgraph 
function getbdg {
    filename=${1##*/}
    # if [[ ! -f  ${1%.bam}.bed ]]; then
    #     echo "Convert bam to bed ! \n"
    #     bamToBed -i $1 | cut -f 1,2,3,4,5,6 | sort -T . -k1,1 -k2,2n -S 5G > ${1%.bam}.bed
    # fi
    #
    x=`samtools view -c ${1}`
    bedtools intersect -sorted -c -b ${1} -a $mywindow | awk -vx=$x '{print $1,$2,$3,$4*1e+06/x}' OFS="\t" > $2.bdg
}

getbdg $edubam ${name}_edu
getbdg $brdubam ${name}_brdu


#------------------------- Operations between edu.bdg and brdu.bdg ---------------------#
paste ${name}_edu.bdg ${name}_brdu.bdg | awk '{if ($8>0 && $4>0){print $1,$2,$3,$4-$8} }' OFS="\t" > ${name}_E-B.bdg 
echo -e "chr\tstart\tstop\t"`ls *E-B.bdg`  | sed 's/\ /\t/g' > ${name}_E-B.txt
cat ${name}_E-B.bdg >> ${name}_E-B.txt


#------- convert to non-slide bw files ---------#
# awk 'NR%5==1'  ${name}_E-B.bdg   > ${name}_E-B_nonslide.bdg
# bedGraphToBigWig ${name}_E-B_nonslide.bdg $hg19 ${name}_E-B_nonslide.bw

#--------- call EdUrich seed -------
signal=${name}_E-B.bdg 
awk '$4>=0.3 {print $1,$2,$3,$4,"s1"}' OFS='\t' $signal > ${name}_edurich_ind.bed
awk '$4<=-0.3 {print $1,$2,$3,$4,"s2"}' OFS='\t' $signal > ${name}_brdurich_ind.bed
awk '{ if($4>-0.3 && $4<0.3) {print $1,$2,$3,$4,"s3"}}' OFS='\t' $signal > ${name}_cold_ind.bed

bedtools cluster -i ${name}_edurich_ind.bed -d $win > ${name}_edurich_cluster_dist${win}.bed
bedtools cluster -i ${name}_brdurich_ind.bed -d $win > ${name}_brdurich_cluster_dist${win}.bed

name2=${name}_edurich
Rscript $dir/Nailseq_callEdUrich_Rsub.R ${name}_edurich_cluster_dist${win}.bed $name2
bedtools merge -i ${name2}_flt_bylen.bed -d 5000 > ${name2}_seed.bed

#name2=${name}_brdurich
#Rscript $dir/Nailseq_callEdUrich_Rsub.R ${name}_brdurich_cluster_dist${win}.bed $name2
#bedtools merge -i ${name2}_flt_bylen.bed -d 5000 > ${name2}_seed.bed

#--------- extending EdUrich ------
$dir/Nailseq_extendEdUrich.py -s $signal -b ${name}_edurich_seed.bed -o ${name}_edurich_v1.bed

bedtools intersect -a ${name}_edurich_v1.bed -b $broadpeak -wa | awk '($3-$2)>20000' > ${name}_edurich_final.bed 

#-------------- clean files ------------#
#${name}_edurich_final.bed 
rm  ${name}_edurich_ind.bed ${name}_brdurich_ind.bed ${name}_cold_ind.bed ${name}_edurich_cluster_dist${win}.bed ${name}_brdurich_cluster_dist${win}.bed ${name2}_seed.bed ${name}_edurich_v1.bed 




