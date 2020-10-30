##################################################
## Project: Ariane's project
## Script purpose: Quant-seq data analysis
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Wed Jan 24 12:50:51 2018
##################################################
library("openxlsx")
require('DESeq2') 

RNAfunctions = "/Volumes/groups/cochella/jiwang/scripts/functions/RNAseq_functions.R"
RNA_QCfunctions =  "/Volumes/groups/cochella/jiwang/scripts/functions/RNAseq_QCs.R"

### data verision and analysis version
version.Data = 'Quantseq_R10331_quantseq_hg'
version.analysis = paste0("_", version.Data, "_20201030")

# Counts.to.Use = "UMIfr"
Save.Tables = TRUE
check.quality.by.sample.comparisons = FALSE


### Directories to save results
design.file = "../exp_design/R10331_quantseq_hg_KO_lines.xlsx"
dataDir = "../../R10331_quantseq_hg/"

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
if(file.exists(design.file)){
  design = read.xlsx(design.file, sheet = 1, colNames = TRUE, skipEmptyRows = TRUE)
  
  design = data.frame(design, stringsAsFactors = FALSE)
  design = design[which(!is.na(design$Sample.ID)), ]
  
  jj = which(colnames(design) == 'Sample.ID')
  design = design[, c(jj, setdiff(c(1:ncol(design)), jj))]
  
  # select columns to keep
  col2keep = c("Sample.ID", "Brief.Sample.Description", "Brief.Sample.Description")
  design = design[, match(col2keep, colnames(design))]
  colnames(design) = c('SampleID', 'strain', 'stage', 'treatment')
  design$condition = NA
  
  Extract.Emilio.Paula.data = FALSE
  if(Extract.Emilio.Paula.data){
    ## Emilio and Paula's data
    jj = which(design$SampleID >= 100003 | design$strain == "MLC1384" | design$strain == "MT17810")
    design$condition[jj] = design$treatment[jj] 
  }else{
    ## just select Philipp's data
    kk = which(design$SampleID < 100003 & design$strain != "MLC1384" & design$strain != "MT17810")
    design = design[kk, ]
  }
 
  ####
  ## manually prepare the design infos
  ####
  design$strain[grep('Arabidopsis', design$strain)] = "Ath"
  design$strain[grep('H2O', design$strain)] = "H2O.control"
  
  design = design[which(design$strain != "Ath" & design$strain != "H2O.control"), ]
  
  design$stage[grep("cells", design$stage)] = "2cells"
  design$stage[grep('fold', design$stage)] = "2.3.fold"
  #design$stage[grep("-", design$stage)] = "none"
  
  design$condition[which(design$strain=="N2")] = "wt"
  design$condition[which(design$strain=="MLC860")] = "pash1.ts"
  design$condition[which(design$strain=="MLC1795")] = "pash1.ts.mirtron"
  design$condition[which(design$strain=="MLC1726")] = "drosha.pash1.aid.RNAi"
  design$condition[which(design$strain=="MLC1729")] = "drosha.pash1.aid.pash1.RNAi.mirtron"
  #design$condition[which(design$strain=="MLC1729" & design$treatment == "pash-1 RNAi")] = "drosha.pash1.aid.pash1.RNAi.mirtron"
  
  design$condition[which(design$strain == "MT14533" & design$treatment == "20 degree")] = 'mir35.ko.20degree'
  design$condition[which(design$strain == "MT14533" & design$treatment == "25 degree")] = 'mir35.ko.25degree'
  #design$condition[which(is.na(design$condition))] = 'none'
  
  design = design[, -which(colnames(design)=="treatment")]
  
}

##########################################
# processing count table, 
# piRNAs total nb of reads and other stat numbers
# spike-in 
##########################################
# table for read counts and UMI
Dir_umi = paste0(dataDir, "R8043_R8521_htseq_counts_BAMs_umi")
Dir_read = paste0(dataDir, "R8043_R8521_htseq_counts_BAMs")

source(RNAfunctions)

aa1 <- list.files(path = Dir_umi, pattern = "*out_umiDedup", full.names = TRUE)
aa1 = merge.countTables.htseq(aa1)
colnames(aa1)[-1] = paste0(colnames(aa1)[-1], ".UMI")

aa2 <- list.files(path = Dir_read, pattern = "*.out", full.names = TRUE)
aa2 = merge.countTables.htseq(aa2)
colnames(aa2)[-1] = paste0(colnames(aa2)[-1], ".readCount")

aa <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "gene", all = TRUE), list(aa1, aa2))

## compare read counts vs. umi counts
source(RNAfunctions)
Compare.UMI.vs.readCounts = FALSE
if(Compare.UMI.vs.readCounts){
  pdfname = paste0(resDir, "readCounts_vs_UMI_normalized", version.analysis, ".pdf")
  pdf(pdfname, width = 10, height = 8)
  
  compare.readCount.UMI(design, aa, normalized = TRUE)
  
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
QC.for.cpm = FALSE
EDA.with.normalized.table = FALSE
Add.miRNA.Targets = TRUE

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
lowlyExpressed.readCount.threshold = 10
require(DESeq2)
source(RNA_QCfunctions)

##########################################
# process miRNA targets 
##########################################
if(Add.miRNA.Targets){
  targets = read.xlsx('../data/examples_miRNA_targets.xlsx', sheet = 2)
  targets = targets[-c(1), -c(1)]
  colnames(targets) = as.character(targets[1, ])
  targets = targets[-c(1:4), ]
  targets = targets[, -ncol(targets)]
  colnames(targets) = c('miR35', 'miR51', 'let7', 'lin4', 'miR58', 'miR1', 'miR35.miRanda')
  
  length(intersect(targets[,1], targets[, 2]))
  length(intersect(targets[,3], targets[, 2]))
}

##########################################
# quality control  
##########################################
if(QC.for.cpm){
  #treat = length(unique(design$treatment[kk]));
  #index.qc = c(3, 5)[which(c(length(unique(design.matrix$genotype)), length(unique(design.matrix$promoter)))>1)]
  # samples.sels = setdiff(c(1:nrow(design)), which(design$condition == "none"))
  
  samples.sels = which(design$stage == "2cells"| design$stage == 'Gastrulation')
  index.qc = c(1, 4, 3)
  
  source(RNA_QCfunctions)
  
  pdfname = paste0(resDir, "/Data_qulity_assessment_2cells_gastrulaiton", version.analysis, "_", Counts.to.Use, ".pdf")
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
# Section for Q1 and Q2
# Pairwise comparisons of mutant, rescue vs wt for different time points, mainly addressing Q1 and Q2
########################################################
########################################################
compares = list(list("2.3.fold", c("pash1.ts.mirtron", "drosha.pash1.aid.pash1.RNAi.mirtron", "drosha.pash1.aid.mirtron")),
                list("L1", c("pash1.ts.mirtron", "drosha.pash1.aid.pash1.RNAi.mirtron", "drosha.pash1.aid.mirtron"))
                #list("Gastrulation", c("pash1.ts", "pash1.ts.mirtron")), 
                #list("Gastrulation", c("drosha.pash1.aid.RNAi", "drosha.pash1.aid.pash1.RNAi.mirtron"))
                )

for(n in 1:length(compares)){
  
  stage.sel = compares[[n]][[1]]
  cond.sel = compares[[n]][[2]]
  compName = paste0(c(stage.sel, cond.sel, "N2"), collapse = "_")
  outDir = paste0(resDir, compName)
  
  cat("time to compare -- ", stage.sel, "\n")
  cat("conditions to compare -- ", cond.sel, "\n")
  cat("output Directory -- ", outDir, "\n")
  if(!dir.exists(outDir)) dir.create(outDir)
  
  # select samples for camparisons
  samples.sels = which(design$condition == "wt")
  for(m in 1:length(cond.sel)) samples.sels = c(samples.sels, which(design$condition == cond.sel[m]))
  samples.sels = intersect(which(design$stage == stage.sel), samples.sels)
  
  # check QC
  pdfname = paste0(outDir, "/Data_QC", version.analysis, "_", Counts.to.Use, "_", compName, ".pdf")
  pdf(pdfname, width = 12, height = 10)
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
  
  Check.RNAseq.Quality(read.count=raw[, samples.sels], design.matrix = design[samples.sels, c(1, 4)])
  # start DE analysis
  dds <- DESeqDataSetFromMatrix(raw[, samples.sels], DataFrame(design[samples.sels, ]), design = ~ condition)
  dds <- dds[ rowSums(counts(dds)) >= lowlyExpressed.readCount.threshold, ]
  dds <- estimateSizeFactors(dds)
  
  cpm = fpm(dds, robust = TRUE)
  colnames(cpm) = paste0(colnames(cpm), ".normDESeq2")
  
  dds = estimateDispersions(dds, fitType = "local")
  plotDispEsts(dds, ylim=c(0.001, 10), cex=1.0)
  abline(h=c(0.1, 0.01), col = 'red', lwd=1.2)
  
  dds = nbinomWaldTest(dds, betaPrior = TRUE)
  resultsNames(dds)
  
  res = cpm
  index.upper = c()
  index.lower = c()
  padj.cutoff = 0.1
  
  for(ii in 1:length(cond.sel)){
    res.ii = results(dds, contrast=c("condition",cond.sel[ii], "wt"))
    summary(res.ii)
    plotMA(res.ii, ylim = c(-2, 2), main = paste0("MA plot -- ", cond.sel[ii], " vs wt"))
    
    ## special cases:  
    # list("Gastrulation", c("pash1.ts", "pash1.ts.mirtron")), 
    # list("Gastrulation", c("drosha.pash1.aid.RNAi", "drosha.pash1.aid.pash1.RNAi.mirtron"))
    if(stage.sel == "Gastrulation" & (cond.sel[ii] == "pash1.ts.mirtron" | cond.sel[ii] == "drosha.pash1.aid.pash1.RNAi.mirtron")){
      upper.ii = which(res.ii$padj > padj.cutoff)
      lower.ii = upper.ii
    }else{
      upper.ii = which(res.ii$log2FoldChange > 0 &  res.ii$padj < padj.cutoff)
      lower.ii = which(res.ii$log2FoldChange < 0 &  res.ii$padj < padj.cutoff)
    }
    
    if(ii == 1) {
      index.upper = c(index.upper, upper.ii)
      index.lower = c(index.lower, lower.ii)
    }else{
      index.upper = intersect(index.upper, upper.ii)
      index.lower = intersect(index.lower, lower.ii)
    }
    
    colnames(res.ii) = paste0(colnames(res.ii), paste0("_", cond.sel[ii], ".vs.N2"))
    res =  data.frame(res, res.ii[, c(1, 2, 5,6)])
  }
  
  dev.off()
  
  write.csv(res, 
            file = paste0(outDir, "/NormalizedTable_DEanalysis_for_", Counts.to.Use, "_", compName, version.analysis, ".csv"), 
            row.names = TRUE)
  
  #write.csv(res[index.upper, ],
  #          file = paste0(outDir, "/NormalizedTable_DEanalysis_Upregulated_for_", Counts.to.Use, "_", compName, version.analysis, ".csv"), 
  #          row.names = TRUE)
  #write.csv(res[index.lower, ], 
  #          file = paste0(outDir, "/NormalizedTable_DEanalysis_Lowregulated_for_", Counts.to.Use, "_", compName, version.analysis, ".csv"), 
  #          row.names = TRUE)
  
  head(res[index.upper, grep("log2Fold|padj", colnames(res))])
  head(res[index.lower, grep("log2Fold|padj", colnames(res))])
  
  if(miRNA.Targets){
    #res.up = res[index.upper, ]
    #res.down = res[index.lower, ]
    
    for(nn in 1:ncol(targets)){
      # n = 1
      tags = targets[,nn]
      tags =  tags[which(!is.na(tags)==TRUE)]
      
      index.all = match(tags, rownames(res)); index.all = index.all[which(!is.na(index.all))]
      #index.up = match(tags, rownames(res.up)); index.up = index.up[which(!is.na(index.up))]
      #index.down = match(tags, rownames(res.down)); index.down = index.down[which(!is.na(index.down))]
      write.csv(res[index.all, ], 
                file = paste0(outDir, "/NormalizedTable_DEanalysis_for_", Counts.to.Use, "_mirTargets_for_", colnames(targets)[nn],  ".csv"), 
                row.names = TRUE)
    }
    
  }
  
}

########################################################
########################################################
# Section Q3: Compare two time points, gastrulation vs. 2cell stage, separately for
# N2, two mutants and two rescues 
# Note:
# 1) after checking the scaterplot, all samples for 2cell stage look very similar, we here just focus on the gastrulation
# 2) so here the code is for either 2 cell or gastrulation stage
########################################################
########################################################
stage.sel = '2cells' # 'Gastrulaiton' or '2cell'
stage.sel = 'Gastrulation' 
cond.sel = c('wt',
             'drosha.pash1.aid.pash1.RNAi.mirtron',
             'drosha.pash1.aid.RNAi', 
             'mir35.ko.20degree' 
             #'mir35.ko.25degree', 
             #'pash1.ts', 'pash1.ts.mirtron'
             )
Post.Candidates.Selection = FALSE

outDir = paste0(resDir, "tables_", stage.sel)

cat("time to compare -- ", stage.sel, "\n")
cat("conditions to compare -- ", cond.sel, "\n")
cat("output Directory -- ", outDir, "\n")

if(!dir.exists(outDir)) dir.create(outDir)

# select samples for camparisons
samples.sels = c()
for(cc in cond.sel) {
  for(tt in stage.sel)
    samples.sels = c(samples.sels, which(design$condition == cc & design$stage == tt))
}

samples.sels = unique(samples.sels)
design.sels = design[samples.sels, ]

# check QC
pdfname = paste0(outDir, "/Data_QC_wt_mutant_rescue_mir35KO_", stage.sel, version.analysis, "_", Counts.to.Use, ".pdf")
pdf(pdfname, width = 12, height = 10)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)

#Check.RNAseq.Quality(read.count=raw[, samples.sels], design.matrix = design[samples.sels, c(1, 4, 3)])

#dev.off()

##  start DE analysis
dds <- DESeqDataSetFromMatrix(raw[, samples.sels], DataFrame(design[samples.sels, ]), design = ~ condition)
dds$condition = relevel(dds$condition, "wt")

dds <- dds[ rowSums(counts(dds)) >= lowlyExpressed.readCount.threshold, ]
dds <- estimateSizeFactors(dds)

cpm = fpm(dds, robust = TRUE)

xx = data.frame(apply(cpm[, which(design.sels$condition=='wt')], 1, mean), 
                apply(cpm[, which(design.sels$condition=="drosha.pash1.aid.RNAi")], 1, mean),
                apply(cpm[, which(design.sels$condition=="drosha.pash1.aid.pash1.RNAi.mirtron")], 1, mean),
                apply(cpm[, which(design.sels$condition=="mir35.ko.20degree")], 1, mean)
                #apply(cpm[, which(design.sels$condition=="mir35.ko.25degree")], 1, mean)
)

colnames(xx) = colnames(xx) = c('wt', 'mutant', 'rescue', 'mir35ko.20degree')
#pairs(log2(xx), upper.panel = panel.fitting, lower.panel=NULL, cex = 0.4, main = 'Gastrulation')

colnames(cpm) = paste0(colnames(cpm), ".normDESeq2")

dds = estimateDispersions(dds)

plotDispEsts(dds, ylim=c(0.001, 10), cex=1.0)
abline(h=c(0.1, 0.01), col = 'red', lwd=1.2)

dds = nbinomWaldTest(dds, betaPrior = TRUE)
resultsNames(dds)

#res = cpm

res.ii = results(dds, contrast=c("condition", 'drosha.pash1.aid.RNAi', 'wt'))
colnames(res.ii) = paste0(colnames(res.ii), "_mutant.vs.wt")
res = data.frame(res.ii[, c(2, 5, 6)])

res.ii = results(dds, contrast=c("condition", 'drosha.pash1.aid.pash1.RNAi.mirtron', 'drosha.pash1.aid.RNAi'))
colnames(res.ii) = paste0(colnames(res.ii), "_rescue.vs.mutant")
res = data.frame(res, res.ii[, c(2, 5, 6)])

res.ii = results(dds, contrast=c("condition", 'drosha.pash1.aid.pash1.RNAi.mirtron', 'wt'))
colnames(res.ii) = paste0(colnames(res.ii), "_rescue.vs.wt")
res = data.frame(res, res.ii[, c(2, 5, 6)])

res.ii = results(dds, contrast=c("condition", 'mir35.ko.20degree', 'wt'))
colnames(res.ii) = paste0(colnames(res.ii), "_mir35KO.20degree.vs.wt")
res = data.frame(res, res.ii[, c(2, 5, 6)])

# res.ii = results(dds, contrast=c("condition", 'mir35.ko.25degree', 'wt'))
# colnames(res.ii) = paste0(colnames(res.ii), "_mir35KO.25degree.vs.wt")
# res = data.frame(res, res.ii[, c(2, 5, 6)])

res.ii = results(dds, contrast=c("condition", 'mir35.ko.20degree', 'drosha.pash1.aid.RNAi'))
colnames(res.ii) = paste0(colnames(res.ii), "_mir35KO.20degree.vs.mutant")
res = data.frame(res, res.ii[, c(2, 5, 6)])

# res.ii = results(dds, contrast=c("condition", 'mir35.ko.25degree', 'drosha.pash1.aid.RNAi'))
# colnames(res.ii) = paste0(colnames(res.ii), "_mir35KO.25degree.vs.mutant")
# res = data.frame(res, res.ii[, c(2, 5, 6)])
# 
# res.ii = results(dds, contrast=c("condition", 'pash1.ts', 'wt'))
# colnames(res.ii) = paste0(colnames(res.ii), "_pash1.ts.vs.wt")
# res = data.frame(res, res.ii[, c(2, 5, 6)])
# 
# res.ii = results(dds, contrast=c("condition", 'pash1.ts.mirtron', 'wt'))
# colnames(res.ii) = paste0(colnames(res.ii), "_pash1.ts.mirtron.vs.wt")
# res = data.frame(res, res.ii[, c(2, 5, 6)])


if(Add.miRNA.Targets){
  tgs.mapping = data.frame(matrix(NA, nrow = nrow(res), ncol = ncol(targets)))
  colnames(tgs.mapping) = paste0('targets_', colnames(targets))
  rownames(tgs.mapping) = rownames(res)
  
  for(n in 1:ncol(targets)){
    # n = 1
    tags = targets[,n]
    tags =  tags[which(!is.na(tags)==TRUE)]
    
    index.all = match(tags, rownames(res)); 
    index.all = index.all[which(!is.na(index.all))]
    #index.up = match(tags, rownames(res.up)); index.up = index.up[which(!is.na(index.up))]
    #index.down = match(tags, rownames(res.down)); index.down = index.down[which(!is.na(index.down))]
    tgs.mapping[index.all, n] = TRUE
  }
}


xx = data.frame(cpm, res, tgs.mapping)
xx = xx[grep('^__', rownames(xx), invert = TRUE), ]

write.csv(xx,
          file = paste0(outDir, "/GeneAll_wt_mutant_rescue_mir35KO_", stage.sel,"_", Counts.to.Use,  version.analysis, ".csv"), 
          row.names = TRUE)


if(Post.Candidates.Selection){
  jj1 = which(res$pvalue_mutant.vs.wt<0.01 & res$log2FoldChange_mutant.vs.wt >0  & 
                res$pvalue_rescue.vs.mutant<0.01 & res$log2FoldChange_rescue.vs.mutant < 0 )
  jj2 = which(res$pvalue_mutant.vs.wt<0.01 & res$log2FoldChange_mutant.vs.wt < 0  & 
                res$pvalue_rescue.vs.mutant<0.01 & res$log2FoldChange_rescue.vs.mutant > 0 )
  
  plot(res$log2FoldChange_mutant.vs.wt[jj1], res$log2FoldChange_rescue.vs.wt[jj1]);
  abline(0, 1, col='red')
  
  plot(res$log2FoldChange_mutant.vs.wt[jj1], res$log2FoldChange_mir35KO.20degree.vs.wt[jj1]);
  abline(0, 1, col='red');abline(h=0, col='red')
  
  plot(res$log2FoldChange_mutant.vs.wt[jj1], res$log2FoldChange_mir35KO.25degree.vs.wt[jj1]);
  abline(0, 1, col='red'); abline(h=0, col='red')
  
  write.csv(xx[jj1, ], 
            file = paste0(outDir, "/GeneList_UpInMutant_DownInRescue_inGastrulation_", Counts.to.Use, "_", compName, version.analysis, ".csv"), 
            row.names = TRUE)
  write.csv(xx[jj2, ],
            file = paste0(outDir, "/GeneList_DownInMutant_UpInRescue_inGastrulation_", Counts.to.Use, "_", compName, version.analysis, ".csv"), 
            row.names = TRUE)
}

dev.off()

##########################################
# mir35/51 targets enrichment test 
##########################################
Enrichment.Analysis.mir35.51.targets = TRUE
if(Enrichment.Analaysis.mir35.51.targets){
  
  ## Hypergeometric test for up-genelist
  jj.ups = jj1
  pvals = rep(NA, ncol(tgs.mapping))
  names(pvals) = colnames(tgs.mapping)
  
  for(n in 1:ncol(tgs.mapping)){
    #n = 1
    nb.target.up = length(which(tgs.mapping[jj.ups, n] == TRUE))
    nb.up = length(jj.ups)
    total.targets = length(which(tgs.mapping[,n] == TRUE))
    total.notargets = nrow(tgs.mapping) - total.targets
    pvals[n] = phyper(q = nb.target.up, m = total.targets, n = total.notargets, k = nb.up, lower.tail = FALSE, log.p = FALSE)
  }
   
  write.csv(pvals,
            file = paste0(outDir, "/Enrichment_pvalues_mirTargets_in_Gene_UpInMutant_DownInRescue_inGastrulation_", Counts.to.Use, "_", compName, version.analysis, ".csv"), 
            row.names = TRUE)
  
  
  ## Subsampling the same nb of genes as mir35 targets 
  yy = data.frame(apply(cpm[, grep("N2", colnames(cpm))], 1, mean), 
                  apply(cpm[, grep("MLC1726", colnames(cpm))], 1, mean),
                  apply(cpm[, grep("MLC1729", colnames(cpm))], 1, mean)
                  )
  colnames(yy) = c('wt', 'mutant', 'rescue')
  yy  = log2(yy+2^-6)
  
  nn = 2
  tags = targets[,nn]
  tags =  tags[which(!is.na(tags)==TRUE)]
  
  index.all = match(tags, rownames(yy)); 
  index.all = index.all[which(!is.na(index.all))]
  
  boxplot(yy[index.all, ])
  
  apply(yy[index.all, ], 2, mean)
  
  nb.sampling = 1000
  
  sims = matrix(NA, ncol = 3, nrow = (nb.sampling+1))
  colnames(sims) = colnames(yy)
  rownames(sims) = c(colnames(targets)[nn], paste0('random_', c(1:nb.sampling)))
  sims[1, ] = apply(yy[index.all, ], 2, mean)
  
  nb.occur = 0
  offset = 0.5
  set.seed(2019)
  for(n in 1:nb.sampling)
  {
    ss = yy[sample(c(1:nrow(yy)), size = length(index.all), replace = FALSE), ]
    avgs = apply(ss, 2, mean)
    sims[(n+1), ] = avgs
    if(avgs[2] > (avgs[1] + 1) & avgs[2] > (avgs[3] + 0.5)) {
      nb.occur = nb.occur +1;
      #print(avgs)
    }
  }
  
  write.csv(sims,
            file = paste0(outDir, "/RandomSubsampling_GeneSet_toTest_", colnames(targets)[nn], ".csv"), 
            row.names = TRUE)
  
}


##########################################
# Q4 for maternal RNA decay
# Step0 define maternal mRNAs using wt
# Step1 identify mRNAs regulated by mir35/51 using 1) wt, mutant and rescue 2) predicted targets from targetscan
##########################################
Test.Candidates.from.Stoeckium = FALSE
Add.mirs.Targets = TRUE
if(Test.Candidates.from.Stoeckium){
  maternal = read.xlsx("../data/Stoeckium-embj-2014.xlsx", sheet = 1, colNames = TRUE)
  mm = match(maternal$GENEWBID, annot$Gene.stable.ID)
  maternal$gene = annot$Gene.name[mm]
  sel = which((as.numeric(maternal$oocyte_rpkm)>=2 | as.numeric(maternal$`1cell_rpkm` >=2)))
  maternal = maternal[sel, ]
  
  ## compare philipp's 2cell stage with Stoeckium et al.
  load(file = "/Volumes/groups/cochella/jiwang/annotations/BioMart_WBcel235.Rdata")
  
  mm1 = match(rownames(res), maternal$GENEID)
  length(which(is.na(mm1)))
  
  mm2 = match(rownames(res), maternal$gene)
  length(which(is.na(mm2)))
  
  pdfname = paste0(outDir, "/MaternalRNA_candidates_Steockius", version.analysis, "_", Counts.to.Use, "_", compName, ".pdf")
  pdf(pdfname, width = 12, height = 8)
  par(cex = 1.0, las = 1, mgp = c(3,0.5,0), mar = c(5,4,4,2), tcl = -0.3)
  
  len = annot$Transcript.length..including.UTRs.and.CDS.[match(rownames(res), annot$Gene.name)]
  plot(apply(res[,c(1:3)], 1, mean)/len*10^3, maternal$`2cell_rpkm`[mm2], log='xy', cex=0.7, main = '2cell.stages.Philipp.vs.Steockius', 
       xlab = 'philip', ylab = 'steoeckius')
  #plot(res[,2]/len*10^3, maternal$`1cell_rpkm`[mm2], log='xy', cex=0.2)
  #plot(res[,1]/len*10^3, maternal$oocyte_rpkm[mm2], log='xy', cex=0.2)
  
  #maternal = data.frame(maternal, stringsAsFactors = FALSE)
  kk = which((as.numeric(maternal$oocyte_rpkm)>=2 | as.numeric(maternal$`1cell_rpkm` >=2)) & 
               as.numeric(maternal$log2FC.mRNA.1cell.ooc) < (-1) & 
               as.numeric(maternal$log2FC.mRNA.1cell.ooc_pv)<0.05)
  
  mats = maternal[kk, ]
  
  mm = match(mats$gene, rownames(res))
  mm = mm[which(!is.na(mm))]
  xx = data.frame(apply(res[mm, c(1:3)], 1, median), apply(res[mm, c(4:6)], 1, median))
  plot(xx[, c(1,2)], log='xy', xlab = 'mean of 2cell', ylab = 'mean of gastrulation')
  abline(0, 1)
  
  hist(res$log2FoldChange_wt_Gastrulation.vs.2cells[mm], main = 'log2FC ')
  
  dev.off()
  
}

# import targets from targetScan
if(Add.mirs.Targets){
  ff = read.table(file = '/Volumes/groups/cochella/jiwang/Databases/TargetScan_Worm_6.2/miR_Family_Info.txt', sep = "\t", header = TRUE)
  targets = read.delim(file = "/Volumes/groups/cochella/jiwang/Databases/TargetScan_Worm_6.2/Predicted_Targets_Info.txt", header = TRUE)
  targets.mir51 = unique(targets$Gene.Symbol[which(targets$miR.Family == 'miR-51/52/53/54/55/56')])
  targets.mir35 = unique(targets$Gene.Symbol[which(targets$miR.Family == 'miR-35/36/37/38/39/40/41-3p/42')])
}

#### Step0 : HERE define the maternal RNAs using N2
#### criterion: significantly decreased from 2cell to gastrulation
pdfname = paste0(outDir, "/Maternal_RNAs_decay_mir35.51", version.analysis, "_", Counts.to.Use, "_", compName, ".pdf")
pdf(pdfname, width = 18, height = 8)
par(cex = 1.0, las = 1, mgp = c(3,1,0), mar = c(5,4,4,2), tcl = -0.3)

head(res[index.wt, grep('.vs.', colnames(res))])
plot(res$log2FoldChange_wt_Gastrulation.vs.2cells[index.wt], -log10(res$pvalue_wt_Gastrulation.vs.2cells[index.wt]))
xx1 = apply(res[, c(1:3)], 1, mean)
xx2 = apply(res[, c(4:6)], 1, mean)

hist(log2(xx2), main = 'log2(N2.Gastrulation)');abline(v=3, col='red')

plot(xx1[index.wt], xx2[index.wt], log='xy', xlab = "N2.2cell", ylab='N2.Gastrulation', cex = 0.5)
abline(0, 1, lwd= 1.0, col='red')
abline(h=2^4, col='red')

length(which(xx2[index.wt]<2^3))
index.mat = index.wt[which(xx2[index.wt]<2^3)]


index.mir35 = match(targets.mir35, rownames(res))
index.mir35 = index.mir35[which(!is.na(index.mir35))]
index.mir51 = match(targets.mir51, rownames(res))
index.mir51 = index.mir51[which(!is.na(index.mir51))]

Save.Res.for.mir51.35.Targets = FALSE
if(Save.Res.for.mir51.35.Targets){
  
  write.csv(res[index.mir35, ], 
            file = paste0(outDir, "/NormalizedTable_DEanalysis_mir35_Targets_for_", Counts.to.Use, "_", compName, version.analysis, ".csv"), 
            row.names = TRUE)
  
  write.csv(res[index.mir51, ], 
            file = paste0(outDir, "/NormalizedTable_DEanalysis_mir51_Targets_for_", Counts.to.Use, "_", compName, version.analysis, ".csv"), 
            row.names = TRUE)
  
  write.csv(res[-c(index.mir35, index.mir51), ],
            file = paste0(outDir, "/NormalizedTable_DEanalysis_others_for_", Counts.to.Use, "_", compName, version.analysis, ".csv"), 
            row.names = TRUE)
  
  
  
}

names(index.all) = cond.sel
index.wt = index.all[[1]]
index.resccue = (intersect(index.all[[2]], index.all[[3]]))
index.mutant = (intersect(index.all[[4]], index.all[[5]]))

length(index.wt)
length(index.resccue)
length(index.mutant)
length(index.xx)

#index.mat = index.wt

## WT vs mutant in which all miRNAs were removed in principle
yy0 = res[, which(colnames(res) == paste0('log2FoldChange_', cond.sel[1], '_Gastrulation.vs.2cells'))]
xlims = c(-7, 6); ylims = c(-7, 6)

yy1 = res[, which(colnames(res) == paste0('log2FoldChange_', cond.sel[4], '_Gastrulation.vs.2cells'))]
yy2 = res[, which(colnames(res) == paste0('log2FoldChange_', cond.sel[2], '_Gastrulation.vs.2cells'))]

yy11 = res[, which(colnames(res) == paste0('log2FoldChange_', cond.sel[5], '_Gastrulation.vs.2cells'))]
yy22 = res[, which(colnames(res) == paste0('log2FoldChange_', cond.sel[3], '_Gastrulation.vs.2cells'))]

par(mfrow = c(1, 2))
plot(yy0[index.mat], yy11[index.mat], cex=0.7, xlim = xlims, ylim = ylims, main = "log2FC.drosha.pash1.aid.RNAi", xlab =  'wt', ylab = 'mutant')
points(yy0[index.mir35], yy11[index.mir35], col = 'blue', pch = 16)
points(yy0[index.mir51], yy11[index.mir51], col = 'darkgreen', pch = 15)
abline(0,1, lwd=1.5, col='darkred')
legend('topleft', legend = c('mir35', 'mir51'), col = c('blue', 'darkgreen'), bty = "n", pch = c(16, 15))

plot(yy0[index.mat], yy22[index.mat], cex=0.7, xlim = xlims, ylim = ylims,  main = "log2FC.drosha.pash1.aid.RNAi", xlab =  'wt', ylab = 'rescue')
points(yy0[index.mir35], yy22[index.mir35], col = 'blue', pch = 16)
points(yy0[index.mir51], yy22[index.mir51], col = 'darkgreen', pch = 15)
abline(0,1, lwd=1.5, col='darkred')
legend('topleft', legend = c('mir35', 'mir51'), col = c('blue', 'darkgreen'), bty = "n", pch = c(16, 15))

par(mfrow = c(1, 1))
plot((yy11 - yy0)[index.mat],  (yy22 - yy0)[index.mat], cex = 0.7, xlab = 'Diff.log2FC.mutant_wt', ylab = 'Diff.log2FC.rescue_wt', main = 'drosha.pash1.aid.RNAi')
points((yy11 - yy0)[index.mir35],  (yy22 - yy0)[index.mir35], col = 'blue', pch = 16)
points((yy11 - yy0)[index.mir51],  (yy22 - yy0)[index.mir51], col = 'darkgreen', pch = 15)
abline(0, 1, lwd = 1.5, col='darkred')
abline(h = 0, lwd =1.5, col='darkred')
abline(v = 0, lwd =1.5, col='darkred')
legend('topleft', legend = c('mir35', 'mir51'), col = c('blue', 'darkgreen'), bty = "n", pch = c(16, 15))

par(mfrow = c(1, 2))
plot(yy0[index.mat], yy1[index.mat], cex=0.7, xlim = xlims, ylim = ylims, main = "log2FC.pash1.ts", xlab =  'wt', ylab = 'mutant')
points(yy0[index.mir35], yy1[index.mir35], col = 'blue', pch = 16)
points(yy0[index.mir51], yy1[index.mir51], col = 'darkgreen', pch = 15)
abline(0,1, lwd=1.5, col='darkred')
legend('topleft', legend = c('mir35', 'mir51'), col = c('blue', 'darkgreen'), bty = "n", pch = c(16, 15))

plot(yy0[index.mat], yy2[index.mat], cex=0.7, xlim = xlims, ylim = ylims, main = "log2FC.pash1.ts", xlab =  'wt', ylab = 'rescue')
points(yy0[index.mir35], yy2[index.mir35], col = 'blue', pch = 16)
points(yy0[index.mir51], yy2[index.mir51], col = 'darkgreen', pch = 15)
abline(0,1, lwd=1.5, col='darkred')
legend('topleft', legend = c('mir35', 'mir51'), col = c('blue', 'darkgreen'), bty = "n", pch = c(16, 15))

par(mfrow = c(1, 1))
plot((yy1 - yy0)[index.mat],  (yy2 - yy0)[index.mat], cex = 0.7, main= 'pash1.ts',  xlab = 'Diff.log2FC.mutant_wt', ylab = 'Diff.log2FC.rescue_wt')
points((yy1 - yy0)[index.mir35],  (yy2 - yy0)[index.mir35], col = 'blue', pch = 16)
points((yy1 - yy0)[index.mir51],  (yy2 - yy0)[index.mir51], col = 'darkgreen', pch = 15)
abline(0, 1, lwd = 1.5, col='darkred')
abline(h = 0, lwd =1.5, col='darkred')
abline(v = 0, lwd =1.5, col='darkred')
legend('topleft', legend = c('mir35', 'mir51'), col = c('blue', 'darkgreen'), bty = "n", pch = c(16, 15))

dev.off()


########################################################
########################################################
# Section : check mir-35 and mir-51 target from modencode
# 
########################################################
########################################################
Add.mirs.Targets = TRUE
# import targets from targetScan
if(Add.mirs.Targets){
  ff = read.table(file = '/Volumes/groups/cochella/jiwang/Databases/TargetScan_Worm_6.2/miR_Family_Info.txt', sep = "\t", header = TRUE)
  targets = read.delim(file = "/Volumes/groups/cochella/jiwang/Databases/TargetScan_Worm_6.2/Predicted_Targets_Info.txt", header = TRUE)
  targets.mir51 = unique(targets$Gene.Symbol[which(targets$miR.Family == 'miR-51/52/53/54/55/56')])
  targets.mir35 = unique(targets$Gene.Symbol[which(targets$miR.Family == 'miR-35/36/37/38/39/40/41-3p/42')])
}

load('/Volumes/groups/cochella/jiwang/Databases/modEncode/modEncode_RNAseq_Supplemental_Table_Unified.Rdata') 

aa = encode[, c(33, 2:19)]

mm1 = match(targets.mir35, aa$gene)
mm2 = match(targets.mir51, aa$gene)
mm1 = mm1[!is.na(mm1)]
mm2 = mm2[!is.na(mm2)]
mm3 = unique(c(mm1, mm2))


write.csv(aa[mm1, ], file = paste0(tabDir, 'modEncode_embryo_mir35_targets.csv'), col.names = TRUE, row.names = FALSE)
write.csv(aa[mm2, ], file = paste0(tabDir, 'modEncode_embryo_mir51_targets.csv'), col.names = TRUE, row.names = FALSE)
write.csv(aa[-mm3, ], file = paste0(tabDir, 'modEncode_embryo_all_genes_other.than.mir35.mir51_targets.csv'), col.names = TRUE, row.names = FALSE)
