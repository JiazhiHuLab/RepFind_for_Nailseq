#!/usr/bin/env python
"""
Update: strict condition for blueprimer: > (n-1)bp match
Date: 20190611
Author: Chen Ai
input: adapter sam file

output:
1. trimmed fastq file
2. barcode_tb


"""
import argparse
import subprocess
import os
import sys
import gzip
import pandas as pd
import re
import pickle
import time
import pysam

parser = argparse.ArgumentParser(description='A program to get barcode, trim fastq.')
parser.add_argument('-A','--Aread', help='Input adapter sam file',dest='adapter',required=True)
parser.add_argument('-B','--base_name', help='Output files base_name',dest='base_name',required=True)
parser.add_argument('-l','--barcode_len', help='Length of barcode', dest='blen',default=14, type=int)
parser.add_argument('-p','--blue_primer', help='Blueprimer sequence', dest='primer', default='TGTAGAGCACGCGTGG')
parser.add_argument('-Z','--gzip',help='Gzip indicator',dest='gzip',action="store_true")
args = parser.parse_args()

t0 = time.time()
def submitter(commander):
        """Submits commands directly to the command line and waits for the process to finish."""
        submiting = subprocess.Popen(commander,shell=True)
        submiting.wait()


complements = {'A':'T', 'T':'A', 'G':'C', 'C':'G','N':'N'}
def reverse_complement(x):
    xrev = x[::-1]
    xrevcomp = ''.join([complements[z] for z in xrev])
    return xrevcomp

# prepare blue primer features
barcode_len = args.blen
bluep = args.primer

# min length of insert sequence
min_len = 40
min_bar_len = barcode_len - 2
min_bluep_len = len(bluep) - 1

# prepare base name
base_name = args.base_name
savefq1 = base_name+"_trimmed_R1.fq"
barcode_dict = base_name+"_read_barcode.pkl"
infile=args.adapter

# open file handles
adp_in = pysam.AlignmentFile(infile, 'r')
save1 = open(savefq1, 'w')

# prepare the reads for output
matchreads = 0
readcount = 0
reads1 = [] # to save trimmed fastq lines
# barcode dict
read_bar = {} # key: rname; value: barcode
#bar_read = {} # key: barcode; value: readname

for read in adp_in.fetch():
    if not read.cigarstring:
        continue
    readcount += 1
    rname = read.qname
    r1 = read.seq
    rmatch1 = read.cigarstring   # example: 17S10M123S

    pattern = re.compile(r"\d+\w")
    m = pattern.findall(rmatch1)   # example: ['17S', '10M', '123S']

    if 'M' in rmatch1:

        # add condition for blueprimer length:
        matched_blue_len = int(re.search(r"\d+", m[1]).group(0) )
        if matched_blue_len < min_bluep_len:
            continue

        # sequence structure: Barcode + Blueprimer + Insert-seq
        # get the barcode position in the R1 read
        end = int(re.search(r"\d+", m[0]).group(0))
        start = max(0, end - barcode_len)
        # get the barcodes
        barcode = r1[start:end]
        read_bar[rname] = barcode

        # add condition for barcode length:
        if (end-start) < min_bar_len :
            continue

        ### trim the random barcode and blue primer of R1
        q = read.qual
        read.seq = read.seq[end + matched_blue_len:]
        read.qual = q[end + matched_blue_len:]
        if read.query_length < min_len:
            continue
        matchreads += 1
        read_for_fq = "@%s\n%s\n+\n%s\n" % (read.qname, read.seq, read.qual)
        reads1.append(read_for_fq)
    else:
        pass

    if readcount == 25000:
        # write fastq file
        save1.writelines(reads1)
        readcount = 0
        reads1 = []


if readcount>0 :
    save1.writelines(reads1)

print ("!!! Blue primer reads: %d \n" % matchreads)

adp_in.close()
save1.close()

with open(barcode_dict, 'wb') as f:
    pickle.dump(read_bar, f, pickle.HIGHEST_PROTOCOL)

# df = pd.DataFrame({'Qname':read_bar.keys(),'Barcode':read_bar.values()})
# df = df.sort_values(by=["Qname"])
# df.to_csv(barcode_tb, sep="\t", index=None)

t1 = time.time()
print ("Time used %s seconds " % (t1-t0))

