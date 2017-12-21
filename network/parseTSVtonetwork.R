args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one arguments must be supplied (input file)", call.=FALSE)
} 
# Load library
library(dplyr)

# Load data
d<-read.delim(file = args[1],header = T,sep = "\t")

# Subset to only the reads that have more than 1 mapping otherwise they have no connections
dsub<-d[d$READ %in% d$READ[duplicated(d$READ)],c(1,3)]

# Create a dataframe with all read and contig-contig combinations
combinations<-merge(x = dsub,y = dsub,by = "READ")
# Get rid of duplicate combinations due to order e.g. 2-3, and 3-2.
combinations[,c(2,3)]<-t(apply(combinations[,c(2,3)],1,sort))
# Count the number of each combination
sum_combinations<-combinations[!duplicated(combinations),] %>% count(CONTIG.x,CONTIG.y)
# Remove connections within the same contig
sum_combinations<-sum_combinations[-which(sum_combinations$CONTIG.x==sum_combinations$CONTIG.y),]
# Remove combinations that have a count of 0 (Not entirely sure that they do occur at all)
sum_combinations<-sum_combinations[which(sum_combinations$n>0),]
# Format output for the network function
output<-paste0(sum_combinations$CONTIG.x,"\t",sum_combinations$CONTIG.y,"\t",sum_combinations$n)
output<-c("scaffold1\tscaffold2\tconnections",output)

if (is.null(args[2])) {
  filename<-"nanopore_mapping_connections.tsv"
} else { filename<- args[2]}

# Save the output to a file
write(x = output,file = filename)

# Test data: 2 reads mapping to four contigs 
# Read 1 maps to contig 1 twice and contig 2
# Read 2 maps to contig 3 and contig 4
# --> 2 connections in total with a count of 1 each
# Desired output from the test data
# scaffold1 scaffold2 connections
# 1 2 1
# 3 4 1
# dat<-data.frame(READ=c(1,1,1,2,2),CONTIG=c(1,1:4))           
# combinations<-merge(x = dat,y = dat,by = "READ")
# combinations[,c(2,3)]<-t(apply(combinations[,c(2,3)],1,sort))
# 
# counts<-combinations[!duplicated(combinations),] %>% count(CONTIG.x,CONTIG.y)
# counts[!(counts$CONTIG.x==counts$CONTIG.y),]