##########################################################################
##########################################################################
# Project: Philipp's small RNA  
# Script purpose: prepare files for GEO submission
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Jun 25 12:56:35 2020
##########################################################################
##########################################################################
## samples to upload for Jorge's project
library(openxlsx)
setwd('/Volumes/groups/cochella/jiwang/Projects/Paula/GEO')

copyDir = 'geo_submission_25june'
fromDir = '../R8024_R7846_R7708/ngs_raw'

## inputs
samples = read.xlsx('sample_to_sumbit.xlsx', colNames = FALSE)
colnames(samples) = 'sampleID'

##########################################
# collect raw data  
##########################################
system(paste0("mkdir -p ", copyDir))
file.format = '*.bam'

files.all = list.files(path = fromDir, pattern = file.format, full.names = TRUE, recursive = TRUE)
files.used = list.files(path = copyDir, pattern = file.format, full.names = TRUE)
sId = samples$sampleID

for(n in 1:length(sId))
{
  kk = grep(sId[n], files.used)
  jj = grep(sId[n], files.all)
  
  cat(sId[n], " -- ")
  if(length(kk) != 1){

    if(length(kk) == 0){
      cat('Not Found -- ')
      if(length(jj)==1) {
        cat('one file to copy')
        system(paste0('cp ', files.all[jj], ' ', copyDir))
      }
      if(length(jj)==0) cat(" >> no file to move")
      if(length(jj) > 1) cat(" >>too many files to copy")
    }
    if(length(kk) > 1){
      cat("Too Many !!", files.used[kk])
    }
  }else{
      cat('copied already')
  }
  
  cat("\n")
}

system(paste0('cp ', files.all[jj], ' ', copyDir))
########################################################
########################################################
# Section : compile the metadata file for raw data
# 
########################################################
########################################################
raw.md5 = read.table("md5sum_raw.txt", header = FALSE, sep = " ")

# here pair-end raw data and 2 fastq for one samples
samples = matrix(NA, ncol = 2, nrow = nrow(raw.md5))
samples = data.frame(samples)
colnames(samples) = c('file', 'md5sum')


for(n in 1:nrow(raw.md5)){
  cat(n, '--')
  cat(raw.md5$V3[n], '\n')
  samples$file[n] = basename(as.character(raw.md5$V3[n]))
  samples$md5sum[n] = as.character(raw.md5$V1[n])
 
}



#samples$file = NA
#samples$md5sum = NA
for(n in 1:nrow(samples)){
  cat(n, '--')
  kk = grep(samples$sampleID[n], raw.md5$V3)
  if(length(kk)==1){
    cat(raw.md5$V3[kk], '\n')
    samples$file[n] = basename(as.character(raw.md5$V3[kk]))
    samples$md5sum[n] = as.character(raw.md5$V1[kk])
  }else{
    cat("ERROR\n")
  }
}


write.xlsx(samples, file = "raw_files_md5sum.xlsx", row.names = FALSE, quote = FALSE)

########################################################
########################################################
# Section : save count table for processed data
# 
########################################################
########################################################
load(file='/Volumes/groups/cochella/jiwang/Projects/Philipp/smallRNA_analysis_philipp/results/R8024_R7846_R7708_all/Design_Raw_readCounts_UMIfr_miRNAs_R8024_R7846_R7708_20190723.Rdata')

read.count = all[, -1];
rownames(read.count) = all$gene
sel.samples.with.spikeIns = match(samples$sampleID, design.matrix$SampleID)

read.count = read.count[, sel.samples.with.spikeIns]
read.count = read.count[-c(9:11), ]

colnames(read.count) = paste0(samples$sampleID, '.umifr.count')
write.table(read.count, file = 'processed/all_samples_umifr_count.txt', sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)

