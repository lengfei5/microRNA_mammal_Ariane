##################################################
## Project: Paula's mir-1 project
## Script purpose: Quant-seq data analysis
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Wed Jan 24 12:50:51 2018
##################################################
require('DESeq2')

RNAfunctions = "/Volumes/groups/cochella/jiwang/scripts/functions/RNAseq_functions.R"
RNA_QCfunctions =  "/Volumes/groups/cochella/jiwang/scripts/functions/RNAseq_QCs.R"

### data verision and analysis version
version.Data = 'Paula_Quantseq_R10603_quantseq_ce_redoDemultiplexing'
version.analysis = paste0("_", version.Data, "_20201201")

# Counts.to.Use = "UMIfr"
Save.Tables = TRUE
check.quality.by.sample.comparisons = FALSE


### Directories to save results
#design.file = "../exp_design/R10331_quantseq_hg_KO_lines.xlsx"
dataDir = "../../../Paula/R10603_quantseq_redemultiplexed_pairend/"

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

write.table(raw, file = paste0(tabDir, "QuantSeq_all_samples_umi_count.txt"), sep = '\t', col.names = TRUE, quote = FALSE, 
          row.names = TRUE)

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

plotMA(dds, ylim = c(-1.5, 1.5))

res = results(dds, contrast=c("condition", 'mir1.mutant', 'wt'))
colnames(res) = paste0(colnames(res), '_', 'mir1.mutant', '.vs.wt')

res = data.frame(res[, c(1, 2, 5, 6)])

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

saveRDS(res, file = paste0(resDir, 'DESeq2_result_saved.rds'))

########################################################
########################################################
# Section : Valcano plot and MA plot with labels
# 
########################################################
########################################################
library(ggplot2)
library(ggrepel)
set.seed(42)

res = readRDS(file = paste0(resDir, 'DESeq2_result_saved.rds'))
res = data.frame(res, stringsAsFactors = FALSE)

res$pvalue.log10 = -log10(res$pvalue_mir1.mutant.vs.wt)
index.hight = unique(c(grep('vha', rownames(res)), which(rownames(res) == 'dct-1'), which(rownames(res) == 'tbc-7')))
  
res$gene = NA
res$gene[index.hight] = rownames(res)[index.hight] 

p1 <- ggplot(res, aes(log2FoldChange_mir1.mutant.vs.wt, pvalue.log10, label = gene)) +
  geom_point(color = ifelse(!is.na(res$gene), "red", "gray70"), size = ifelse(!is.na(res$gene), 3, 1)) +
  geom_text_repel(force = 5,
                  #arrow = arrow(length = unit(0.03, "npc"), type = "closed", ends = "first"),
                  direction     = "both", 
                  point.padding = 0.1,
                  segment.color = "grey40") 

#p1 + ggsave(paste0(resDir, "volcano_plot_vha_complex_subunits.pdf"), width = 10, height = 5, dpi = "print")

res$mean.log2 = log2(res$baseMean_mir1.mutant.vs.wt)
ggplot(res, aes(mean.log2, log2FoldChange_mir1.mutant.vs.wt, label = gene)) +
  geom_point(color = ifelse(!is.na(res$gene), "red", "gray70"), size = ifelse(!is.na(res$gene), 3, 1)) +
  geom_text_repel(force = 1,
                  #arrow = arrow(length = unit(0.03, "npc"), type = "closed", ends = "first"),
                  direction     = "both", 
                  point.padding = NA,
                  segment.color = "grey40") +
  geom_hline(yintercept=0, linetype="solid", color = "black")

########################################################
########################################################
# Section : Gene set enrichment analysis 
# here using the package B "enrichplot" recommended by Maria Fisher 
########################################################
########################################################
library(enrichplot)
library(clusterProfiler)
library(openxlsx)
library(ggplot2)

#res = read.csv(file = paste0(resDir, "DESeq2.norm_all.mir1KO.vs.wt_", Counts.to.Use,  version.analysis, ".csv"), row.names = 1)
res = readRDS(file = paste0(resDir, 'DESeq2_result_saved.rds'))
geneList = res$log2FoldChange_mir1.mutant.vs.wt
names(geneList) = rownames(res)
geneList = geneList[which(!is.na(geneList))]
geneList = geneList[order(-geneList)]


##########################################
# vha complex
##########################################
vha = read.xlsx(paste0(resDir, 'vha_genes.xlsx'), sheet = 1)

vha.complex = data.frame(term = 'vha.complex', gene = vha$vha.genes, stringsAsFactors = FALSE)

gsea.vha <- GSEA(geneList, TERM2GENE = vha.complex, TERM2NAME = NA, nPerm = 10000)
#edo2 <- gseNCG(geneList, nPerm=1000)
gseaplot1 = gseaplot2(gsea.vha, geneSetID = 1, title = 'vha complex', pvalue_table = TRUE, subplots = c(1:2),
          ES_geom = "line")

gseaplot1 +
  ggsave(paste0(resDir, "GSEA_vha.complex.pdf"), width = 10, height = 5, dpi = "print")

##########################################
# mir1 targets
##########################################
targets = read.xlsx(paste0(resDir, 'mir1_predicted_targets_TargetScanworm.xlsx'), sheet = 1)
targets = targets[!is.na(match(targets$Gene.symbol, names(geneList))), ]

#kk = which(as.numeric(targets$Aggregate.PCT) >0.7)

mir1.targets = data.frame(term = 'mir1.targets', gene = targets$Gene.symbol,
                          stringsAsFactors = FALSE)
mir1.targets = mir1.targets[which(!is.na(mir1.targets$gene)), ]

xx = res
xx = xx[match(mir1.targets$gene, rownames(xx)), ]
length(which(xx$pvalue_mir1.mutant.vs.wt<0.05 & xx$log2FoldChange_mir1.mutant.vs.wt>0))
length(which(xx$pvalue_mir1.mutant.vs.wt<0.05 & xx$log2FoldChange_mir1.mutant.vs.wt<0))

write.table(xx, file = paste0(resDir, 'mir1.targets_up_downregulated.txt'), sep = '\t', col.names = TRUE, row.names = TRUE)


gsea.mir1 <- GSEA(geneList, TERM2GENE = mir1.targets, TERM2NAME = NA, nPerm = 1000)

#edo2 <- gseNCG(geneList, nPerm=1000)
gseaplot2 = gseaplot2(gsea.mir1, geneSetID = 1, title = 'mir1.targets', pvalue_table = TRUE, subplots = c(1:2),
                      ES_geom = "line")

gseaplot2 +
  ggsave(paste0(resDir, "GSEA_mir1.targets.pdf"), width = 10, height = 5, dpi = "print")



##########################################
# go term enrichment analysis
##########################################
library(enrichplot)
library(clusterProfiler)
library(openxlsx)
library(ggplot2)
library(stringr)
library(org.Ce.eg.db)

annot = read.csv(file = "/Volumes/groups/cochella/jiwang//annotations/BioMart_WBcel235_noFilters.csv", 
                 header = TRUE, stringsAsFactors = FALSE)

annot$uniprot = annot$Gene.description
xx = annot$uniprot
annot$uniprot = sapply(xx, function(x) unlist(strsplit(str_match(x, "(?<=\\[).+?(?=\\])")[1, 1], ':'))[3])

ensname = annot$Gene.stable.ID[match(rownames(res), annot$Gene.name)]
uniprotnames = annot$uniprot[match(rownames(res), annot$Gene.name)]

pval.cutoff = 0.05

pdfname = paste0(resDir, "/GO.Terms_enrichment_mir1.mutant_upregulated_downregulated_genes.pdf")
pdf(pdfname, width = 12, height = 10)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)

kk.up = which(res$pvalue_mir1.mutant.vs.wt < pval.cutoff & res$log2FoldChange_mir1.mutant.vs.wt > 0)
kk.down = which(res$pvalue_mir1.mutant.vs.wt < pval.cutoff & res$log2FoldChange_mir1.mutant.vs.wt < 0)
cat('nb of upregulated genes : ', length(kk.up), '\n')
cat('nb of downregulated genes : ', length(kk.down), '\n')

ego.up <-  enrichGO(gene         = ensname[kk.up],
                 universe     = ensname,
                 OrgDb         = org.Ce.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
head(ego.up)
barplot(ego.up, showCategory=20) + ggtitle("go term for upregulated genes")

ego.down <-  enrichGO(gene         = ensname[kk.down],
                    universe     = ensname,
                    OrgDb         = org.Ce.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
head(ego.down)
barplot(ego.down, showCategory=20) + ggtitle("go term for downregulated genes")

write.csv(ego.up, file = paste0(resDir, "GO_term_enrichmenet_for_upregulated_genes_pval_0.05.csv"), 
          row.names = TRUE)

write.csv(ego.down, file = paste0(resDir, "GO_term_enrichmenet_for_downregulated_genes_pval_0.05.csv"), 
          row.names = TRUE)

# kegg.up <- enrichKEGG(gene         = uniprotnames[kk.up],
#                  organism     = 'cel',
#                  universe     = uniprotnames,
#                  keyType = 'uniprot',
#                  pvalueCutoff = 0.05)
# head(kegg.up)
# barplot(kegg.up, showCategory=20) + ggtitle("kegg for upregulated genes")
# 
# 
# kegg.down <- enrichKEGG(gene         = uniprotnames[kk.down],
#                       organism     = 'cel',
#                       universe     = uniprotnames,
#                       keyType = 'uniprot',
#                       pvalueCutoff = 0.2)
# head(kegg.down)
# 
# barplot(kegg.down, showCategory=20) + ggtitle("kegg for upregulated genes")

dev.off()


