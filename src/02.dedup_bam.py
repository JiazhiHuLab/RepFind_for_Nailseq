#!/usr/bin/env python2
##################################
#  dedup the trimmed alignment
#  save the barcode and read relationship
#  Author: Chen Ai
#  Update: 2019-04-30
#  Date: May 18th, 2018
##################################
import argparse
import subprocess
import os
import sys
import gzip
import pickle
import re
import time
import pysam
import pandas as pd

parser = argparse.ArgumentParser(description='A program to remove pcr duplication in the trimmed sam.')
parser.add_argument('-b', '--barcode', help="input the read name and barcode pkl.", dest="barcode", required=True)
parser.add_argument('-i', '--inputBam', help="input the trimmed bam file to dedup ", dest="bamfile", required=True)
parser.add_argument('-L', '--logfile', help="log file", dest="logfile")
args = parser.parse_args()

t0 = time.time()

# import read_bar dict: {readname:barcode_seq}
x = []
with open(args.barcode, 'rb') as f:
    try:
        x.append(pickle.load(f))
    except EOFError:
        sys.exit("Fail to open %s \n" % args.barcode)
read_bar = x[0]

# prepare output file
basename = args.bamfile[:-4]
outf1 = basename + "_dedup.bam"
outf2 = basename + "_dup.bam"

in1 = pysam.AlignmentFile(args.bamfile, 'rb')
out_unique = pysam.AlignmentFile(outf1, 'wb', template= in1 )
out_dup = pysam.AlignmentFile(outf2, 'wb', template= in1 )


bcd_reads = {}  # barcode : read_name
bcd_counts = {}
totreads = 0
rmreads = 0


for read in in1.fetch():
    if not read.cigarstring: # whether the read is aligned , if not aligned cigarstring is none
        continue
    rname = read.qname
    r1_chr = read.reference_id
    r1_start = read.reference_start
    r1_end = read.reference_end

    bcd = read_bar[rname]
    pos = tuple([r1_chr, r1_start, r1_end])
    try:
        bcd_counts[bcd] += 1
        for bpos in bcd_reads[bcd]:
            if bpos == pos:
                out_dup.write(read)
                rmreads += 1
                break
        else:
            out_unique.write(read)
            bcd_reads[bcd].append((r1_chr, r1_start, r1_end))

    except KeyError:
        bcd_counts[bcd] = 1
        bcd_reads[bcd] = [(r1_chr, r1_start, r1_end)]
        out_unique.write(read)
    totreads += 1


in1.close()
out_unique.close()
out_dup.close()

logout = open(args.logfile,'w')
print('In the remove duplication ' +  \
    'total=' + str(totreads) + \
    '\nremove=' + str(rmreads), file=logout)

logout.close()

t1 = time.time()
print ("Time used %d seconds " % (t1-t0))
