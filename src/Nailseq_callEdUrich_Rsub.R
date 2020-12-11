#!/usr/bin/Rscript
# for statistic of edu rich region
#setwd("/home/workdir/chen/proj/work/LY_related/DRC-seq/EBdouble/gm78_notreat_calldomain")
args = commandArgs(trailingOnly = T)
file=args[1]
outname=args[2]
.libPaths()
#library(tidyverse)
suppressPackageStartupMessages(require(tidyverse))

# outstat
summ_out = paste0(outname,"_summary_of_cluster_wincout.txt")
outbed = paste0(outname,"_flt_bylen.bed")

df= read_tsv(file, col_names = F)
cut_off=1  # 5x5000 = 25000

win_count = df %>% group_by(X6) %>%
  summarise(count=n())
#summary(win_count$count)
filter_win_count = win_count %>% filter(count> cut_off)
filter_win_count
sm=summary(filter_win_count$count)
cat(sm,file = summ_out)

df2 = df %>% filter(X6 %in% filter_win_count$X6)
write_tsv(df2, path = outbed, col_names = F)

# df3 = df %>% filter( !(X6 %in% filter_win_count$X6))
# head(df3)
# df3 = df3 %>% filter(X4>1.20)  # 1.20 is from quantile(0.999)
# write_tsv(df3, "gm78_071_edurich_cluster_dist5000_flt_bysignal.bed", col_names = F)
