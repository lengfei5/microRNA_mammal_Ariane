##################################################
## Project: Paula's mir-1 project
## Script purpose: Quant-seq data analysis
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Wed Jan 24 12:50:51 2018
##################################################
#library("openxlsx")
require('DESeq2') 

RNAfunctions = "/Volumes/groups/cochella/jiwang/scripts/functions/RNAseq_functions.R"
RNA_QCfunctions =  "/Volumes/groups/cochella/jiwang/scripts/functions/RNAseq_QCs.R"

### data verision and analysis version
version.Data = 'Paula_Quantseq_R10603_quantseq_ce'
version.analysis = paste0("_", version.Data, "_20201127")

# Counts.to.Use = "UMIfr"
Save.Tables = TRUE
check.quality.by.sample.comparisons = FALSE


### Directories to save results
#design.file = "../exp_design/R10331_quantseq_hg_KO_lines.xlsx"
dataDir = "../../../Paula/R10603_quantseq_redo_pairend/"

resDir = paste0("../results/", version.Data, "/")
tabDir =  paste0(resDir, "tables/")
RdataDir = paste0(resDir, "/Rdata/")

if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}
if(!dir.exists(RdataDir)){dir.create(RdataDir)}

##########################################
# Import Sample information and table of read counts
# mainly manully 
##########################################
design = data.frame(SampleID = c('135694', '135695', '135696', '135697'), 
                    condition = c('wt', 'wt', 'mir1.mutant', 'mir1.mutant'))


##########################################
# processing count table of umi and read
##########################################
# table for read counts and UMI
Dir_umi = paste0(dataDir, "htseq_counts_BAMs_umi")
Dir_read = paste0(dataDir, "htseq_counts_BAMs")

source(RNAfunctions)

aa1 <- list.files(path = Dir_umi, pattern = "*umiDedup.txt", full.names = TRUE)
aa1 = merge.countTables.htseq(aa1)
colnames(aa1)[-1] = paste0(colnames(aa1)[-1], ".UMI")

aa2 <- list.files(path = Dir_read, pattern = "*.txt", full.names = TRUE)
aa2 = merge.countTables.htseq(aa2)
colnames(aa2)[-1] = paste0(colnames(aa2)[-1], ".readCount")

aa <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "gene", all = TRUE), list(aa1, aa2))

aa = aa[grep('^__', aa$gene, invert = TRUE), ]

## compare read counts vs. umi counts
source(RNAfunctions)
Compare.UMI.vs.readCounts = TRUE
if(Compare.UMI.vs.readCounts){
  pdfname = paste0(resDir, "readCounts_vs_UMI_normalized", version.analysis, ".pdf")
  pdf(pdfname, width = 10, height = 8)
  
  compare.readCount.UMI(design, aa, normalized = FALSE)
  
  dev.off()
}

#save(design, aa, file = paste0(RdataDir, 'Design_Raw_readCounts_UMI_All_incl_Paula_Emilio', version.analysis, '.Rdata'))
save(design, aa, file=paste0(RdataDir, 'Design_Raw_readCounts_UMI', version.analysis, '.Rdata'))

######################################
######################################
## Section: spike-in and piRNA normalization
# optionally double check the data quality with sample comparisons
## save the normalized tables
######################################
######################################
Counts.to.Use = "UMI"
QC.for.cpm = TRUE
EDA.with.normalized.table = FALSE

load(file=paste0(RdataDir, 'Design_Raw_readCounts_UMI', version.analysis, '.Rdata'))
source(RNAfunctions)
source(RNA_QCfunctions)

if(Counts.to.Use == 'readCounts'){
  all = process.countTable(all=aa, design = design, special.column = ".readCount", ensToGeneSymbol = TRUE)
}else{
  if(Counts.to.Use == "UMI"){
    all = process.countTable(all=aa, design = design, special.column = "UMI", ensToGeneSymbol = TRUE)
  }else{
    cat("Error : no counts found for ", Counts.to.Use, "for miRNAs \n")
  }
}

all = all[which(!is.na(all$gene)), ]
raw = ceiling(as.matrix(all[, -1]))
raw[which(is.na(raw))] = 0
rownames(raw) = all$gene

###
### specify parameters for DESEeq2 and pairwise comparisons
###
require(DESeq2)
source(RNA_QCfunctions)
lowlyExpressed.readCount.threshold = 5

##########################################
# quality control  
##########################################
if(QC.for.cpm){
  #treat = length(unique(design$treatment[kk]));
  #index.qc = c(3, 5)[which(c(length(unique(design.matrix$genotype)), length(unique(design.matrix$promoter)))>1)]
  # samples.sels = setdiff(c(1:nrow(design)), which(design$condition == "none"))
  
  samples.sels = c(1:nrow(design))
  index.qc = c(1, 2)
  
  source(RNA_QCfunctions)
  
  pdfname = paste0(resDir, "/Data_qulity_assessment_AllSamples", version.analysis, "_", Counts.to.Use, ".pdf")
  pdf(pdfname, width = 12, height = 10)
  
  Check.RNAseq.Quality(read.count=raw[, samples.sels], design.matrix = design[samples.sels, index.qc])
  
  dev.off()
  
}

##########################################
# calculate scaling factor and normalization
##########################################
if(EDA.with.normalized.table){
  require(DESeq2)
  samples.sels = setdiff(c(1:nrow(design)), which(design$condition == "none"))
  
  raw = ceiling(as.matrix(all[, (samples.sels+1)]))
  raw[which(is.na(raw))] = 0
  rownames(raw) = all$gene
  
  #source(RNAfunctions)
  dds <- DESeqDataSetFromMatrix(raw, 
                                DataFrame(design[samples.sels, ]), 
                                design = ~ condition + stage)
  lowlyExpressed.readCount.threshold = 20
  dds <- dds[ rowSums(counts(dds)) >= lowlyExpressed.readCount.threshold, ]
  dds <- estimateSizeFactors(dds)
  
  fpm = fpm(dds, robust = TRUE)
  
  if(Save.Tables){
    xx = data.frame(fpm, stringsAsFactors = FALSE)
    write.csv(xx, file = paste0(tabDir, "Table_normalized_for_", Counts.to.Use,  version.analysis, ".csv"), 
              row.names = TRUE)
  }
  
}

########################################################
########################################################
# Section: pairwise comparisons, each condition vs UN 
# 
########################################################
########################################################
samples.sels = c(1:nrow(design))
design.sels = design[samples.sels, ]

pdfname = paste0(resDir, "/Data_KO_vs_WT_", Counts.to.Use, ".pdf")
pdf(pdfname, width = 12, height = 10)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)

##  start DE analysis
dds <- DESeqDataSetFromMatrix(raw[, samples.sels], DataFrame(design[samples.sels, ]), design = ~ condition)
dds$condition = relevel(dds$condition, "wt")

dds <- dds[ rowSums(counts(dds)) >= lowlyExpressed.readCount.threshold, ]
dds <- estimateSizeFactors(dds)

cpm = fpm(dds, robust = TRUE)

#colnames(xx) = colnames(xx) = c('wt', 'mutant', 'rescue', 'mir35ko.20degree')
#pairs(log2(xx), upper.panel = panel.fitting, lower.panel=NULL, cex = 0.4, main = 'Gastrulation')

colnames(cpm) = paste0(colnames(cpm), ".normDESeq2")

dds = estimateDispersions(dds, fitType = 'local')

plotDispEsts(dds, ylim=c(0.001, 10), cex=1.0)
abline(h=c(0.1, 0.01), col = 'red', lwd=1.2)

dds = nbinomWaldTest(dds, betaPrior = TRUE)
resultsNames(dds)

res = results(dds, contrast=c("condition", 'mir1.mutant', 'wt'))
colnames(res) = paste0(colnames(res), '_', 'mir1.mutant', '.vs.wt')

res = data.frame(res[, c(2, 5, 6)])

plot(res$log2FoldChange_mir1.mutant.vs.wt, -log10(res$pvalue_mir1.mutant.vs.wt), cex = 0.7)
kk = grep('vha|dct-1|tbc-7', rownames(res))
points(res$log2FoldChange_mir1.mutant.vs.wt[kk], -log10(res$pvalue_mir1.mutant.vs.wt)[kk], col = 'red', cex = 1.2, pch = 16)
text(res$log2FoldChange_mir1.mutant.vs.wt[kk], -log10(res$pvalue_mir1.mutant.vs.wt)[kk], labels = rownames(res)[kk])
#xx = res
#xx = xx[order(-res$log2FoldChange_mir1.mutant.vs.wt), ]

xx = data.frame(cpm, res, stringsAsFactors = FALSE)
xx = xx[order(-xx$log2FoldChange_mir1.mutant.vs.wt), ]

write.csv(xx,
          file = paste0(resDir, "DESeq2.norm_all.mir1KO.vs.wt_", Counts.to.Use,  version.analysis, ".csv"), 
          row.names = TRUE)


dev.off()


########################################################
########################################################
# Section : Gene set enrichment analysis 
# here using the package B "enrichplot" recommended by Maria Fisher 
########################################################
########################################################
library(enrichplot)
library(clusterProfiler)
library(DOSE)
data(geneList)



edo2 <- gseNCG(geneList, nPerm=1000)

gseaplot2(edo2, geneSetID = 1, title = edo2$Description[1], pvalue_table = TRUE)




