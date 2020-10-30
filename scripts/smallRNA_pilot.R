##########################################################################
##########################################################################
# Project: Ariane's small RNA-seq in human cell line
# Script purpose:
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Feb 25 15:32:02 2020
##########################################################################
##########################################################################
library("openxlsx")
require('DESeq2')
miRNAfunctions = "/Volumes/groups/cochella/jiwang/scripts/functions/miRNAseq_functions.R"
RNA_QCfunctions =  "/Volumes/groups/cochella/jiwang/scripts/functions/RNAseq_QCs.R"
source(miRNAfunctions)

### data verision and analysis version
version.Data = 'miRNAs_mouse_R10227'
version.analysis = paste0("_", version.Data, "_20200928")

Counts.to.Use = "readCounts"
Save.Tables = FALSE
check.quality.by.sample.comparisons = FALSE
Save.Tables.correctedBackground = TRUE

### Directories to save results
design.file = "../../R10227_miRNA_mm/NGS_design_R10227.xlsx"
dataDir = "../../R10227_miRNA_mm/nf_processing_srbc/result"
# srbcDir = paste0(dataDir, "R8024_R7846_R7708_srbc_v2")
#oldDir = paste0(dataDir, "R8024_R7846_R7708_old_v2")

resDir = "../results/miRNA_mouse_R10227/"
tabDir =  paste0(resDir, "tables/")

if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}

##########################################
# Import Sample information and table of read counts
# mainly manully 
##########################################
if(file.exists(design.file)){
  design = read.xlsx(design.file, sheet = 1, colNames = TRUE)
  design = data.frame(design, stringsAsFactors = FALSE)
  design = design[c(1:9),]
  
  jj = which(colnames(design) == 'Sample.ID')
  design = design[, c(jj, setdiff(c(1:ncol(design)), jj))]
  design = data.frame(design$Sample.ID, design$Genetic.Background, design$Brief.Sample.Description)
  
  design$cells = NA
  design$treatment = NA
  
  design$cells[grep('whole spleen', design$design.Brief.Sample.Description)] = 'spleen'
  design$cells[grep('CD43', design$design.Brief.Sample.Description)] = 'CD43.spleen'
  design$cells[grep('GFP', design$design.Brief.Sample.Description)] = 'GFP.spleen'
  design$treatment[grep(', oxidized', design$design.Brief.Sample.Description)] = 'oxidized'
  design$treatment[grep(', not oxidized', design$design.Brief.Sample.Description)] = 'not.oxidized'
  
  
  # design$batch = sapply(design$design.Brief.Sample.Description, function(x) gsub('#', '',unlist(strsplit(as.character(x), ' '))[1]))
  # design$genotype = sapply(design$design.Brief.Sample.Description, function(x) unlist(strsplit(as.character(x), ' '))[2])
  # design$genotype[grep('24', design$genotype)] = '24'
  # kk = which(design$genotype != 'wt')
  # design$genotype[kk] = paste0("henn1.exp_", design$genotype[kk], 'h')
  # design$treatment = 'treatment'
  # design$treatment[grep('-', design$design.Brief.Sample.Description)] = 'notreatment'
  # 
  # design = design[, c(1, 3:5)]
  design = design[, c(1, 4:5)]
  colnames(design) = c('SampleID', 'cells', 'treatment')
  
}

##########################################
# processing count table, 
# piRNAs total nb of reads and other stat numbers
# spike-in 
##########################################
# table for read counts and UMI
aa = read.delim(paste0(dataDir,  "/countTable.txt"), sep = "\t", header = TRUE)

# spike-ins 
spikes = read.delim(paste0(dataDir,  "/spikeInTable.txt"), sep="\t", header = TRUE, row.names = 1)

source(miRNAfunctions)
pdfname = paste0(resDir, "readCounts_vs_UMI", version.analysis, ".pdf")
pdf(pdfname, width = 10, height = 6)
compare.readCounts.umiFr.umiNum(design, aa, spikes)

dev.off()

# start to process table and merge them into one
# kk = grep("piRNA_|pash", aa$Name)
# if(length(kk)>0){ piRNAs = aa[kk, ]; aa = aa[-kk,];}

source(miRNAfunctions)

if(Counts.to.Use == 'readCounts'){
  all = process.countTable(all=aa, design = design, select.counts = "Total.count")
}else{
  if(Counts.to.Use == "UMIfr"){
    all = process.countTable(all=aa, design = design, select.counts = "Total.UMIfr.count")
  }else{
    cat("Error : no counts found for ", Counts.to.Use, "for miRNAs \n")
  }
}

xx = as.matrix(all[, -1])
xx[which(is.na(xx))] = 0
stat.miRNAs = (apply(xx, 2, sum))

spikes = data.frame(gene=rownames(spikes), spikes, stringsAsFactors = FALSE)
if(Counts.to.Use == "readCounts"){
  spikes = process.countTable(all=spikes, design = design, select.counts = "Total.spikeIn")
}else{
  if(Counts.to.Use == "UMIfr"){
    spikes = process.countTable(all=spikes, design = design, select.counts = "Total.UMI.spikeIn")
  }else{
    cat("Error : no counts found for ", Counts.to.Use, "for spike-ins \n")
  }
}

total.spikes = floor(apply(as.matrix(spikes[, -1]), 2, sum))

#xx = rbind(pash, all)
#xx = rbind(mirtrons, xx)
xx = rbind(spikes, all)

all = xx 

design.matrix = data.frame(design, stat.miRNAs, total.spikes)

save(all, design.matrix, file=paste0(resDir, 'Design_Raw_readCounts_', Counts.to.Use,  version.analysis, '.Rdata'))


########################################################
########################################################
# Section : Quality control
# 1) spik-in check
# 2) read counts normalized by DESeq2
# 3) read counts normalized by spike-ins 
########################################################
########################################################
load(file=paste0(resDir, 'Design_Raw_readCounts_',Counts.to.Use,  version.analysis, '.Rdata'))

read.count = all[, -1];
sel.samples.with.spikeIns = c(1:nrow(design))

raw = floor(as.matrix(read.count[, sel.samples.with.spikeIns]))
raw[which(is.na(raw))] = 0
rownames(raw) = all$gene
#dds <- DESeqDataSetFromMatrix(raw, DataFrame(design.matrix), design = ~ treatment + stage)

##########################################
# spike-in normalization
##########################################
index.spikeIn = grep("spikeIn", rownames(raw))
spike.concentrations = c(0.05, 0.25, 0.5, 1.5, 2.5, 3.5, 5, 25) ## the concentration is amol/mug-of-total-RNA

o1 = order(design.matrix$cells, design$treatment)

## calculate scaling factor using spike-ins
source(miRNAfunctions)

pdfname = paste0(resDir, "/Spike_in_signals_normalized_DESeq_", Counts.to.Use,  version.analysis, ".pdf")
pdf(pdfname, width = 16, height = 10)
par(mfrow=c(2,2))

res.spike.in = calculate.scaling.factors.using.spikeIns(raw[, o1], concentrations = spike.concentrations, index.spikeIn = index.spikeIn, read.threshold = 5)

dev.off()

norms = res.spike.in$norms4DESeq2
cpm = res.spike.in$cpm;
res = res.spike.in$normalization.spikeIn
colnames(cpm) = paste0(colnames(cpm), ".cpm")
colnames(res) = paste0(colnames(res), ".amol.per.mugRNA.normBySpikeIns")

#res = data.frame(res)
write.csv(res, file = paste0(tabDir, paste0("normalized_signals_scalling.factors_using_spikeIns", version.analysis, ".csv")), 
          row.names = TRUE)
write.csv(cpm, file = paste0(tabDir, paste0("cpm", version.analysis, ".csv")), 
          row.names = TRUE)

##########################################
# QC for read counts normalized with DESeq2
##########################################
source(RNA_QCfunctions)

index.qc = c(1, 3, 4)
kk = order(design.matrix$genotype, design$treatment, design.matrix$batch)

pdfname = paste0(resDir, "Data_qulity_assessment_DESeq2", version.analysis, ".pdf")
pdf(pdfname, width = 12, height = 10)
Check.RNAseq.Quality(read.count=raw[, kk], design.matrix = design.matrix[kk, index.qc])
dev.off()

kk = order(design.matrix$genotype, design$treatment, design.matrix$batch)
ss = apply(raw[,kk], 2, sum)
kk = kk[which(ss>500000)]

pdfname = paste0(resDir, "Data_qulity_assessment_filteredSamples_DESeq2", version.analysis, ".pdf")
pdf(pdfname, width = 12, height = 10)
Check.RNAseq.Quality(read.count=raw[, kk], design.matrix = design.matrix[kk, index.qc])
dev.off()

pdfname = paste0(resDir, "Data_qulity_assessment_filteredSamples_SpikeIns", version.analysis, ".pdf")
pdf(pdfname, width = 12, height = 10)
Check.RNAseq.Quality(read.count=raw[, kk], design.matrix = design.matrix[kk, index.qc], norms = norms[kk])
dev.off()

########################################################
########################################################
# Section : Enrichment analysis for different time points
# 
########################################################
########################################################
load(file=paste0(resDir, 'Design_Raw_readCounts_',Counts.to.Use,  version.analysis, '.Rdata'))
design = design.matrix
read.count = all[, -1];
read.count = floor(as.matrix(read.count))
read.count[which(is.na(read.count))] = 0
rownames(read.count) = all$gene


ss = apply(read.count, 2, sum)
kk = which(ss>500000)
design = design[kk, ]
read.count =  read.count[, kk]

tcs = unique(design$genotype)
length(tcs)

require('DESeq2')
Check.data.quality = TRUE
Filter.lowly.expressed.using.predefined.miRNA.list = FALSE
Save.Comparison = TRUE
index.N2 = which(design$genotype=="wt") 
read.count = read.count[grep('spikeIn', rownames(read.count), invert = TRUE), ] # not consider spike-in here

for(n in c(1:length(tcs)))
{
  # n = 6
  specifity = tcs[n];
  
  specDir = paste0(tabDir, specifity, "/")
  if(!dir.exists(specDir)){dir.create(specDir)}
  
  kk = which(design$genotype==specifity)
  # kk = kk[order(design$promoter[kk])]
  
  # find genome types and promoters
  genos = unique(design$genotype[kk]);
  #promoters = unique(design$promoter[kk]); 
  
  if(specifity != 'wt') { kk = unique(c(kk, index.N2)) }
  
  design.matrix = data.frame(sample=colnames(read.count)[kk], design[kk, ])
  raw = as.matrix(read.count[,kk])
  
  ## save plots for enrichment analysis
  pdfname = paste0(specDir, "QCs_assessment_Enrichment_analysis_", specifity, ".pdf") #save all plots during data processing
  pdf(pdfname, width = 12, height = 8)
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
  
  #treat = length(unique(design$treatment[kk]));
  #index.qc = c(3, 5)[which(c(length(unique(design.matrix$genotype)), length(unique(design.matrix$promoter)))>1)]
  if(length(unique(design.matrix$genotype))>1){
    index.qc = c(1, 4, 5)
  }else{
    index.qc = c(1, 5)
  }
  
  design.matrix.QC = design.matrix[, index.qc]
  
  source(RNA_QCfunctions)
  if(Check.data.quality){ Check.RNAseq.Quality(read.count=read.count[, kk], design.matrix = design.matrix.QC);}
  
  #dev.off()
  
  countData = raw
  if(specifity != 'wt'){
    dds <- DESeqDataSetFromMatrix(countData, DataFrame(design.matrix), design = ~ genotype + treatment + genotype:treatment)
    dds$genotype <- relevel(dds$genotype, ref=specifity);
    dds$treatment = relevel(dds$treatment, ref="notreatment");
  }else{
    dds <- DESeqDataSetFromMatrix(countData, DataFrame(design.matrix), design = ~ treatment )
  }
  
  ## normalization read counts and cpm
  dds= estimateSizeFactors(dds)
  cpm = fpm(dds, robust = TRUE)
  
  dds = estimateDispersions(dds, fitType = "parametric")
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
  plotDispEsts(dds, ylim=c(0.001, 10), cex=1.0)
  
  colnames(cpm) = paste0(colnames(cpm), "_",  dds$treatment,  '.normalized.read.counts.DESeq2')
  #colnames(cpm0) = paste0(colnames(cpm0), "_",  dds$treatment,  '.cpm')
  dds = nbinomWaldTest(dds, betaPrior = FALSE)
  resultsNames(dds)

  if(specifity != 'wt')
  {
    ## the ones showing difference treated vs untreated in WT background
    res = results(dds, contrast=c("treatment","treatment","notreatment"))
    summary(res)
    res0 = res;
    colnames(res0) = paste0(colnames(res0), ".without.wt")
    kk.mature = c(1:nrow(res))
    plot(res$log2FoldChange[kk.mature], -log10(res$pvalue[kk.mature]), xlab='log2(FoldChange)', ylab='-log10(pvalue)', cex=0.8, 
         main=paste0(specifity, "--",  " (NOT using WT)"))
    abline(v=0, lwd=2.0, col='black')
    abline(h=c(5, 10), lwd=1.0, lty=2, col='blue')
    text(res$log2FoldChange[kk.mature], -log10(res$pvalue[kk.mature]), rownames(res)[kk.mature], cex=0.7, offset = 0.3, pos = 3)
    
    ## the one showing N2-specific treatment effect is smaller than WT-specific treatment effect 
    ## (refer to the DEseq2 mannul or tutorial https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
    res2 = results(dds, name = "genotypewt.treatmenttreatment", lfcThreshold = 0, altHypothesis = "less")
    
    ## merge both res1 and res2
    res1 = as.data.frame(res); 
    res2 = as.data.frame(res2);
    
    ## change the sign (log2(FC.WT.treated) - log2(FC.WT.untreated)) - (log2(FC.N2.treated) - log2(FC.N2.untreated))
    res2$log2FoldChange = -res2$log2FoldChange; 
    ii.enrich = which(res1$log2FoldChange>0) 
    res1[ii.enrich, ] = res2[ii.enrich, ] ## replaced by res2 if enriched; keep if depleted
    
    ## merge the table with comparison with N2 and the one without
    res = data.frame(res1, res0, stringsAsFactors = FALSE)
    
  }else{
    res <- results(dds, contrast = c("treatment", "treatment", "notreatment"));
    summary(res)
    res = as.data.frame(res);
  }
  
  kk.mature = c(1:nrow(res))
  plot(res$log2FoldChange[kk.mature], -log10(res$pvalue[kk.mature]), xlab='log2(FoldChange)', ylab='-log10(pvalue)', 
       cex=0.8, main=paste0(specifity))
  abline(v=seq(-1, 1, by=0.5), lwd=1.0, col='gray')
  abline(h=c(3, 5, 10), lwd=1.0, lty=2, col='blue')
  text(res$log2FoldChange[kk.mature], -log10(res$pvalue[kk.mature]), rownames(res)[kk.mature], cex=0.7, offset = 0.3, pos = 3)
  #ii = match(c('lsy-6', 'mir-791', 'mir-790', 'mir-793'), rownames(res));
  #points(res$log2FoldChange[ii], -log10(res$pvalue)[ii], cex=1.5, col='darkgreen', pch=16)
  
  ## save the comparion in excel form
  if(Save.Comparison){
    #res.sig <- subset(res, padj < FDR.cutoff)
    res.sig = data.frame(cpm, res, stringsAsFactors = FALSE);
    #res.sig = res;
    #res.sig = as.data.frame(res.sig);
    res.sig = res.sig[order(-res$log2FoldChange), ]
    #res = res[o1, ]
    #cpm = cpm[o1, ]
    #kk.mature = grep("cel-",rownames(res.sig), invert = TRUE)
    #kk.star = grep("cel-",rownames(res.sig), invert = FALSE)
    write.csv(res.sig, 
              file=paste0(specDir, 'miRNA_Enrichment_Analysis_', specifity, version.analysis, '.csv'),     
              row.names = TRUE, quote = FALSE)
    #write.csv(res.sig[kk.star, ], 
    #          file=paste0(specDir, 'miRNA_Enrichment_Analysis_', specifity, '_', cc, '_', prot, '_Star_', version.analysis,'.csv'),     
    #          row.names = TRUE, quote = FALSE)
  }
  
  dev.off();
  
}


########################################################
########################################################
# Section : legacy code
# 
########################################################
########################################################

## enrichment analysis is done with samples with the same genotype and promoter
for(cc in genos)
{
  # cc = genos[1]; 
  cat(specifity, "--", cc, "--\n")
  
  jj = which(design.matrix$genotype==cc)
  
  if(length(jj)>=2)
  {
    ## add N2 genotype if WT bakcground for interaction term; 
    if(cc=="WT" & specifity != "whole.body") jj = unique(c(jj, which(design.matrix$genotype=="N2"))) 
    
    
    
    ## filter lowly expressed miRNA first before estimating scaling factor and dispersion parameters
    jj.expressed = NULL
    jj.expressed = match(rownames(dds), expressed.miRNAs$miRNA)
    sels = !is.na(jj.expressed)
    cat("nb of expressed miRNAs --", sum(sels), "\n")
    #cat("nb of expressed genes -- ", sum(expressed.miRNAs$mature[sels]))
    dds = dds[sels, ]
    index.sel = jj.expressed[sels]
    jj.mature = expressed.miRNAs$mature[index.sel]
    rownames(dds)[jj.mature] = as.character(expressed.miRNAs$gene[index.sel[jj.mature]])
    cpm0 = cpm0[sels, ]
    
    #cat("size factor is -- ", sizeFactors(dds), "\n")
    
    ## estimate scaling factor and dispersion parameter
    kk.mature = grep("cel-",rownames(dds), invert = TRUE)
    dds_sf0 = dds[kk.mature, ]
    dds_sf0 <- estimateSizeFactors(dds_sf0)
    sizeFactors(dds) = sizeFactors(dds_sf0)
    cat("size factor is -- ", sizeFactors(dds), "\n")
    
    
    
    #ii.test = which(rownames(cpm)=="lsy-6"| rownames(cpm)=="cel-lsy-6-3p");
    #log2(mean(cpm[ii.test, c(3:4)])/mean(cpm[ii.test, c(1:2)]))
    
  }
  
}