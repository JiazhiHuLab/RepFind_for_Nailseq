#!/usr/bin/env python
from optparse import OptionParser
import pandas as pd
import numpy as np
import pybedtools
from pybedtools import BedTool


parser = OptionParser()
parser.add_option('-s','--sig', dest='signalfile', help='Input the E-B.bdg for signal file.')
parser.add_option('-b','--bed', dest='bedfile', help='Input the bed file of E-B>0.3 . which is output of ')
parser.add_option('-o','--out', dest='outbed', help='Output of the extending EdUrich peak.')

options, args = parser.parse_args()


hg19="/home/workdir/chen/databases/genomes/hg19.chrom.sizes"
cutoff_intvl = 100000 # take the interval < 100k into consideration

sig = pybedtools.BedTool(options.signalfile)
seed = pybedtools.BedTool(options.bedfile)
out = options.outbed

#----- get interval between two EdUrich seeds -----#
seed = seed.sort(g=hg19)
bw_seed = seed.complement(g=hg19)
bw_ov = bw_seed.intersect(sig, wao=True)
group_bw_ov = bw_ov.groupby(g=[1,2,3] , c=7, o='sum', prec=2)
group_bw_ov2 = bw_ov.groupby(g=[1,2,3] , c=7, o='collapse', prec=2)

df_v1 = group_bw_ov.to_dataframe()
df_v2 = group_bw_ov2.to_dataframe()
df = pd.concat([df_v1.iloc[:,:3], df_v1.name, df_v2.name], axis=1)
df.columns = ['chrom', 'start','end','sum','score']

# limit length: filter peaks the too large:
df  = df.loc[(df.end-df.start)<= cutoff_intvl, :]

# --- First: selected by the interval sum:
df['sum'] = df['sum'].replace(".", np.nan)
df['sum'] = df['sum'].astype(float).fillna(0.0)

# filter by sum of E-B in the interval
df2 = df.loc[df['sum'] > 0,:]

# ---- Second:  filter by number < 2
tmp = df['score'].str.split(',').tolist()
c=[0]*len(tmp)
for i,ilist in zip(range(len(tmp)),tmp):
    try:
        for j in ilist:
            if float(j)<= 0 :
                c[i]+=1
    except (TypeError, ValueError):
        c[i] = 10
df['lt0_count'] = c
df3 = df.loc[df['lt0_count'] <=2, :]

# ---- Third: union of two selections:
x=list(set(df2.index) | set(df3.index))
df_final = df.loc[x,:]

bw_final = BedTool.from_dataframe(df_final.iloc[:,:3])

# merge the original seeds and the intervals with <2 E-B
peaks = BedTool.cat(seed,bw_final).sort(g=hg19).merge(d=0)
peaks.saveas(out)



