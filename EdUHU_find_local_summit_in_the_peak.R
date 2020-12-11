#!/usr/bin/env Rscript
# Date: 20190607
# Author: Chen Ai
# FOR: find local maxima in each ERIZ of HUEdU
# input: 1. HUEdU bdg 2. ERIZ
# output: local peak maxima of ERIZ
# Lessons: local maxima: 1 deviate , sign, You need to learn the essence of the job
suppressPackageStartupMessages(library(rtracklayer))

get_local_maxima <- function(x, cutoff=0.01){
  # diff(x) 一阶导
  delta = diff(x)
  # 一阶导符号的变化 +1 到-1 是 -2 求出最大值
  si = sign(delta)
  idx = which(diff(si)== -2) +1
  
  # 但是有一些小突起，避免影响结果:设置一个cutoff，步长为5k内变化未超过的就不算
  final_idx = c()
  for (i in idx){
    up = x[i] - x[max(1,i-5)]
    down = x[i] - x[min(length(x), i+5)]
    m = mean(up , down)
    if (m > cutoff){
      final_idx = append(final_idx, i)
    }
  }
  #idx = which(diff(sign(diff(x)))==-2) + 1
  return(final_idx)
}

get_smooth <- function(bdg, eriz){
  e1 = bdg[findOverlaps(bdg, eriz)@from]
  if (width(eriz) > 400000 ){
    lspan = 0.2
  }else{
    lspan = 0.4
  }
  df = data.frame(x=start(e1), y = e1$score)
  res = loess(y~x, df, span=lspan)
  
  # get local minima 
  idx = get_local_maxima(res$fitted)
  peaks = e1[idx]
  peaks = resize(peaks, fix="center", width = 5000)
  eriz_peaks = peaks
  
  # get smooth bdg
  e2 = e1
  e2$score = res$fitted
  eriz_loes_bdg = e2
  return(list(peak=peaks, bdg=e1, smoothbdg=e2))
}


#setwd("/home/workdir/chen/proj/recent/analysis/initiationSignals")
message("#-------------\nProgram: finding local maximum in a peak. \n
  Usage: Rscript thisscript.R win1k.bdg peaks.bed outname\n")
args = commandArgs(trailingOnly = T)
if (length(args)!= 3){
  message("Program: finding local maximum in a peak. \n
  Usage: Rscript thisscript.R win1k.bdg peaks.bed outname\n")
}


bdg = import.bedGraph(args[1])
eriz = import.bed(args[2])
name = args[3]
eriz$name = paste("peak", seq_len(length(eriz)), sep="_")

eriz_bdg = GRangesList()
eriz_peaks = GRangesList()
eriz_loes_bdg = GRangesList()
for (i in seq_len(length(eriz))){
  res = get_smooth(bdg=bdg, eriz=eriz[i])
  eriz_bdg[[i]] = res$bdg
  eriz_loes_bdg[[i]] = res$smoothbdg
  eriz_peaks[[i]] = res$peak
}
names(eriz_bdg) = paste("ERIZ", seq_len(length(eriz)), sep="_")
names(eriz_peaks) = paste("ERIZ", seq_len(length(eriz)), sep="_")
names(eriz_loes_bdg) = paste("ERIZ", seq_len(length(eriz)), sep="_")

b1 = makeGRangesFromDataFrame(as.data.frame(eriz_bdg), keep.extra.columns = T) 
b2 = makeGRangesFromDataFrame(as.data.frame(eriz_peaks), keep.extra.columns = T) 
b3 = makeGRangesFromDataFrame(as.data.frame(eriz_loes_bdg), keep.extra.columns = T) 

export.bedGraph(b1,  paste0(name,"_inPeak.bdg"))
export.bed(b2, paste0(name, "local_max_inPeak.bed"))
export.bedGraph(b3, paste0(name, "smoothed_in_loes.bdg"))


