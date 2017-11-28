#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Two arguments must be supplied (input file, and name)", call.=FALSE)
} 
df<-read.delim(args[1])
cat(paste0(args[2],",", median(df$End-df$Start),"\n"))
