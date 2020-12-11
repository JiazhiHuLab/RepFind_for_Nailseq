#!/bin/bash
# Date:2019-04-30
# Author: Chen Ai

set -euo pipefail

SCRIPT=`basename ${BASH_SOURCE[0]}`

# Initiate the variables to default values.
name=NA 
max_threads=4
barcode_len=17
blue_primer=TGTAGAGCACGCGTGG

# prepare the scripts
script_dir=`dirname ${BASH_SOURCE[0]}`
py_trim_se=${script_dir}/src/01.trim_reads_get_bar_inR1.py
py_dedup_se=${script_dir}/src/02.dedup_bam.py
#py_filter_len=/home/workdir/chen/proj/scripts/Nailseq_Fetseq_common/03.fileter_bam_by_length.py


# Help function
function HELP {
    echo -e "Program: Nailseq pipeline"
    echo -e "Version: 3.0"
    echo -e "Author: Chen Ai"
    echo -e "\nUsage: $SCRIPT [options] "
    echo -e "\nOptions:  "
    echo -e "-1 --Input fq1, REQUIRED"
    echo -e "-2 --Input fq2, REQUIRED"
    echo -e "-g  --Genome [hg19/hg38/mm9/mm10], REQUIRED"
    echo -e "-n  --Name, prefix name of the output files. REQUIRED"
    echo -e "-o  --OutputDIR, DIR of output. REQUIRED"
    echo -e "-p  --Maximum_threads, default: 4 "
    echo -e "-r  --Length of random molecular barcode: default 17"
    echo -e "-l  --Ligased blue primer sequence: defualt: TGTAGAGCACGCGTGG"
    echo -e "-h  --Displays this help message.\n"
    echo -e "Example: $SCRIPT -f1 R1.fq -f2 R2.fq -g hg19 -n test -o test_dir "
    exit 1
}

#Check the number of arguments. If none are passed, print help and exit.
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
            fq1=$OPTARG 
            echo -e "Input Fq1: $fq1"
            ;;
        2) 
            fq2=$OPTARG
            echo -e "Input Fq2: $fq2"
            ;;
        g) 
            genome=$OPTARG
            echo -e "Genome is $genome"
            ;;
        n)  name=$OPTARG
            ;;
        o)  outpath=$OPTARG
            ;;
        p)  max_threads=$OPTARG
            ;;
        r)  barcode_len=$OPTARG
            ;;
        l)  blue_primer=$OPTARG  
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

if [ ! -f $fq1 ] ; then
    echo -e  "$fq1 does not exist ! "
    exit 2
fi
if [ ! -f $fq2 ] ; then
    echo -e  "$fq2 does not exist ! "
    exit 2
fi
if [ ! -d $outpath ]; then
    mkdir -p $outpath
fi


#---------- set genome assembly ---------#
case $genome in 
    hg19) 
        bw2_idx=/home/ubuntu/genomes/bowtie2_indexes/hg19
        bwa_ref=/home/mengzhu/database/bwa_indexes/hg19/hg19
        chrom_ref=/home/ubuntu/genomes/hg19/annotation/ChromInfo.txt
        ;;
    hg38)
        bwa_ref=/home/mengzhu/database/bwa_indexes/hg38/hg38
        chrom_ref=/home/ubuntu/genomes/hg38/annotation/ChromInfo.txt
        ;;
    mm10)
        bw2_idx=/home/ubuntu/genomes/bowtie2_indexes/mm10
        bwa_ref=/home/mengzhu/database/bwa_indexes/mm10/mm10
        chrom_ref=/home/ubuntu/genomes/mm10/annotation/ChromInfo.txt
        ;;
    mm9)
        # mouse genome : mm9
        bw2_idx=/home/ubuntu/genomes/bowtie2_indexes/mm9
        bwa_ref=/home/workdir/chen/databases/bwa_indexes/mm9.fa
        chrom_ref=/home/ubuntu/genomes/mm9/annotation/ChromInfo.txt
        ;;
    *)
        echo -e "Error: Unknown genome assembly"
        HELP
        exit 2
        ;;
esac   



fq1=$(readlink -e $fq1)
fq2=$(readlink -e $fq2)

outpath_abs=$(readlink -e $outpath)
cd $outpath_abs

tmppath=$outpath_abs/tmp 
if [ ! -d $tmppath ]; then
    mkdir -p $tmppath
fi

qcpath=$outpath_abs/qc
if [ ! -d $qcpath ]; then
    mkdir -p $qcpath
fi


#------ !!! IMPORTANT step: Stitch the Raw Reads -----#
pear -v 10 -j $max_threads -f $fq1 -r $fq2 -o $tmppath/$name
### output the following files:
ass_se=$tmppath/${name}.assembled.fastq
unass_r1=$tmppath/${name}.unassembled.forward.fastq
unass_r2=$tmppath/${name}.unassembled.reverse.fastq


#------ MAP blue primer to reads --------#
# prepare the next output file names
tmp_ass_bluesam=$tmppath/${name}_ass_blue.sam
tmp_unass_blue1=$tmppath/${name}_unass_blue1.sam 
tmp_unass_blue2=$tmppath/${name}_unass_blue2.sam 

bowtie2 --local --no-1mm-upfront -D 15 -R 2 -N 1 -L 10 -i C,6 --ma 2 --mp 6,2 --np 1 --rdg 5,3 --rfg 5,3 --score-min C,20 -p 4 --reorder -t \
        -x /home/workdir/chen/proj/recent/adapter_index/bluep -U $ass_se -S  $tmp_ass_bluesam
bowtie2 --local --no-1mm-upfront -D 15 -R 2 -N 1 -L 10 -i C,6 --ma 2 --mp 6,2 --np 1 --rdg 5,3 --rfg 5,3 --score-min C,20 -p 4 --reorder -t \
        -x /home/workdir/chen/proj/recent/adapter_index/bluep -U $unass_r1 -S $tmp_unass_blue1
bowtie2 --local --no-1mm-upfront -D 15 -R 2 -N 1 -L 10 -i C,6 --ma 2 --mp 6,2 --np 1 --rdg 5,3 --rfg 5,3 --score-min C,20 -p 4 --reorder -t \
        -x /home/workdir/chen/proj/recent/adapter_index/bluep -U $unass_r2 -S $tmp_unass_blue2


#----- TRIM the blue primer and adapter sequence; and get the bacode -----#fter Stitched reads, if Read2 also has blue primer, remove this reads
samtools view -h -F 16 $tmp_ass_bluesam > $tmppath/${name}_ass_blue_filtered.sam
rm $tmp_ass_bluesam
#----- TRIM the blue primer and adapter sequence; and get the bacode -----#
ass_base=$tmppath/${name}_ass
unass_base=$tmppath/${name}_unass

python $py_trim_se -A $tmppath/${name}_ass_blue_filtered.sam  -B $ass_base  -l $barcode_len -p  $blue_primer


#-------- MAPPING ASSEMBLE using se mode --------# 
# the results all include "trimmed"
ass_trim_fq=${ass_base}"_trimmed_R1.fq"
ass_tb_bar=${ass_base}"_read_barcode.pkl"
ass_trim_sam=${ass_base}"_trimmed.sam"
ass_trim_bam=${ass_base}"_trimmed.bam"


bwa mem -k 20 -t $max_threads $bwa_ref $ass_trim_fq > $ass_trim_sam 2> $tmppath/bwa.ass.log
samtools view -bS -q 30 $ass_trim_sam | samtools  sort -o $ass_trim_bam
rm $ass_trim_sam
samtools index $ass_trim_bam
python $py_dedup_se -b $ass_tb_bar -i $ass_trim_bam -L ${tmppath}/${name}_ass_dedup.log
# result: bed[:-4]"_trimmed_dedup.bed"  bed[:-4]"_trimmed_dup.bed" 

ass_dedup_bam=${ass_base}"_trimmed_dedup.bam" #final 
ass_dedup_bdg=${ass_base}"_trimmed_dedup_all.bdg"
ass_dedup_plus_bdg=${ass_base}"_trimmed_dedup_plus.bdg"
ass_dedup_minus_bdg=${ass_base}"_trimmed_dedup_minus.bdg"

samtools index $ass_dedup_bam
# -------- by strand ---------#
bedtools genomecov -bg -ibam $ass_dedup_bam -g $chrom_ref > $ass_dedup_bdg
bedtools genomecov -bg -strand + -ibam $ass_dedup_bam -g $chrom_ref > $ass_dedup_plus_bdg
bedtools genomecov -bg -strand - -ibam $ass_dedup_bam -g $chrom_ref > $ass_dedup_minus_bdg




#-------- quality control ---------#
fastqc -o $qcpath -t $max_threads $fq1 $fq2 $ass_trim_fq 
samstat  $ass_trim_bam $ass_dedup_bam 
bamtools stats -in $ass_dedup_bam > ass_bam_stat.txt

zcat $fq1 | awk 'NR%4==2{c++; }
          END{
                print c;
              }' > ${name}"_1_fastq_R1_stat.txt"

zcat $fq2 | awk 'NR%4==2{c++; }
        END{
              print c;
            }' > ${name}"_1_fastq_R2_stat.txt"

samtools view -c $ass_trim_bam > ${name}"_2_ass_before_dedup_stat.txt"
samtools view -c $ass_dedup_bam > ${name}"_3_ass_dedup_stat.txt"



#--------Clean files -----------#
#samtools index $ass_dedup_bam
mv ${ass_dedup_bam}* $ass_dedup_bdg $ass_dedup_plus_bdg $ass_dedup_minus_bdg $outpath_abs
mv ass_bam_stat.txt unass_bam_stat.txt $outpath_abs
#rm -rf tmp/*.fastq tmp/*fail*.fq tmp/*.sam 
rm -rf tmp 
