#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("One arguments must be supplied (output file)", call.=FALSE)
} 
df1<-read.delim("temp/CDSmedian.csv",header = T,sep = ",")
df2<-read.delim("temp/checkm_binstats_parsed.csv",header = T,sep = ",")
df3<-read.delim("temp/found16S.csv",header = T,sep = ",")
df4<-read.delim("temp/checkm.tsv",header = T,sep = "\t",strip.white = T,check.names = T)
dc1<-merge.data.frame(x = df1,y = df2,by = "bin")
dc2<-merge.data.frame(x = df3,y = df4,by = "bin")
dc3<-merge.data.frame(x = dc1,y = dc2,by = "bin")
write.csv(x = dc3,file = args[1],row.names = F)
