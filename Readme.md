# Transcription shapes DNA replication initiation to preserve genome integrity.

We develop a high-throughput nucleoside analog incorporation sequencing assay, NAIL-seq, and identify thousands of early replication initiation zones (ERIZs) in both mouse and human cells. The identified ERIZs fall in open chromatin compartments while are mutually exclusive with transcription elongation and occupy mainly non-transcribed regions adjacent to transcribed regions. Furthermore, we reveal that RNA polymerase II actively redistributes the chromatin-encircled mini-chromosome maintenance (MCM) complex but not the origin-recognition complex (ORC) to actively restrict early DNA replication initiation outside of transcribed regions. The coupling of RNA polymerase II and MCM is further validated by detected MCM accumulation and DNA replication initiation after RNA polymerase II stalling via anchoring nuclease-dead Cas9 at the transcribed genes. Importantly, we also find that the orchestration of DNA replication initiation by transcription can efficiently prevent gross DNA damage.

## Datasets
Liu Y, Ai C, Gan T, et al., Hu J: Transcription shapes DNA replication initiation to preserve genome integrity. [GSE174680](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174680)

## RepFind for Nail-seq analysis

Nail-seq is a sequencing method to identify the early replication initiation zones with two strategies.
The pipeline is under GNU General Public License v3.0.

- EdU/BrdU dual labeling sequencing.
- EdU-seq under HU treatment.


## EdU/BrdU dual labeling sequencing

### Step 1: mapping to the genome

Example: 

- For EdU sample: 

```
Nailseq_preprocess_stitch.sh -1 EdU_R1.fq -2 EdU_R2.fq -g hg19 -n EdU -o OutDir -p 4  
```

- For BrdU sample:  

```
Nailseq_preprocess_stitch.sh -1 BrdU_R1.fq -2 BrdU_R2.fq -g hg19 -n BrdU -o OutDir -p 4 
```
(Please set the genome assembly configuration in line 100-119 of Nailseq_preprocess_stitch.sh) 


### Step 2: detect broad peak by SICER

Example:

```
Nailseq_call_broadpeak.sh EdU.bam Ctrl.bam EdU_broad
```
The control sample is the input DNA without enrichment.

(Please install a SICER.py in python 2.7 environment)

### Step 3: detect the EdU-BrdU rich peak and signal 
Use the bam file got from Step 1 and broadPeak got from Step2 for calliing EdU-BrdU rich peak.

Example:
```
Nailseq_call_edu-brdu_signal_peak.sh -1 EdU.bam -2 BrdU.bam -g hg19 -n EdU-BrdU -w 5000 -b EdU_broad.peak
```

### Step 4: detect the local maxima of peak
Because inside the peaks, there are some the multiple summits.
This program is to identify the local maximum peak within a broad peak. 

Algorithm: 
First, the signal in the peak is smoothed using loess smooth function.
Then the difference between the neighboring smoothed signal is calculated from the signal (using R function diff).
Finally, the location of the local maxima is the transition sites of the difference (from positive to negative).



## EdU-seq under HU treatment
### Step 1: mapping to genome
- For EdU_HU sample: 
```
Nailseq_preprocess_stitch.sh -1 EdUHU_R1.fq -2 EdUHU_R2.fq -g hg19 -n EdU -o OutDir -p 4  
```

### Step 2: detect the Peaks
First, the peaks are detected by macs14 with narrowPeak mode.
Next, the peaks are merged < maxgap are merged together
```
Nailseq_call_EdUHU_peak.sh EdUHU.bam Ctrl.bam HUSeq EdUHU.bw
```


## Installation

Dependencies:

#### Python
packages: argparse, subprocess, gzip, scipy, pandas, re, pickle, pysam >= 0.15.4 
SCICER.py(https://github.com/dariober/SICERpy)

#### R
tidyverse


## Usage

### Nailseq\_preprocess_stitch.sh
```
Nailseq_preprocess_stitch.sh

Program: Nailseq pipeline
Version: 1.0
Author: Chen Ai

Usage: Nailseq_preprocess_stitch.sh [options]

Options:
-1 --Input fq1, REQUIRED
-2 --Input fq2, REQUIRED
-g  --Genome [hg19/hg38/mm9/mm10], REQUIRED
-n  --Name, prefix name of the output files. REQUIRED
-o  --OutputDIR, DIR of output. REQUIRED
-p  --Maximum_threads, default: 4
-r  --Length of random molecular barcode: default 17
-l  --Ligased blue primer sequence: defualt: TGTAGAGCACGCGTGG
-h  --Displays this help message.

```

### Nailseq_call\_edu-brdu\_signal\_peak.sh
```
Program: RepFind Nail-seq EdU-BrdU peak
Version: 1.0
Author: Chen Ai

Usage: Nailseq_call_edu-brdu_signal_peak.sh [options]

Options:
-1 --Input EdU bamfile, REQUIRED
-2 --Input BrdU bamfile, REQUIRED
-g  --Genome [hg19/mm9/mm10], REQUIRED
-n  --Output prefix of EdU-BrdU, name_EdU-BrdU.bdg . REQUIRED
-w --Window [INT]bp default 5000 bp.
-b --Broadpeak of output Nailseq_call_broadPeak.sh
-h  --Displays this help message.

Example: Nailseq_call_edu-brdu_signal_peak.sh -1 EdU.bam -2 BrdU.bam -g hg19 -o EdU-BrdU.peak

```

