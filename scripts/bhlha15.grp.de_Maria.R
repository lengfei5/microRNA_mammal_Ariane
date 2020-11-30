################################################################################
# Maria Fischer
# Institute of Molecular Pathology (IMP), Vienna, Austria
# 2020/09/21
################################################################################
# Bhlha15 in-vitro stimulated plasmablasts vs wt and
# Xbp1 in-vitro stimulated plasmablasts vs wt

# ------------------------------------------------------------------------------
# setup
# ------------------------------------------------------------------------------
# packages
library(DESeq2)
library(readr)
library(stringr)
library(dplyr)
library(tidyr)
library(apeglm)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(enrichplot)

source("/groups/busslinger/Maria/Projects/scripts/pairs_funcs.R")
source("/groups/busslinger/Maria/Projects/scripts/DESeq2.readFeatureCounts.R")

prjDir <- "/groups/busslinger/Maria/Projects/Miriam/2020_0922/"
opf <- "bhlha15.grp.de" # output prefix

# these also include IG/TCR genes but RNA-Seq data did not consider them
xbp1.p2g <- "/groups/busslinger/jaritz/Busslinger/2018Nov12/results/peaks/macs2/full/29435ab.14951/0/p2g/p2g/p2g/refseq.mm9.2014_0110.imp.withIgTcr.merged.segmentedByTranscripts.2.5k_0k/assignment.txt"
bhlha15.atac_chip.file <- str_c(prjDir,"data/Bhlha15_ATAC_down_ChIP_bound.txt")

# colors
color.bhlha15 <- "darkgreen"
color.e2.2    <- "orange"
color.xbp1    <- "purple"
color.xbp1.unbound <- "#e2bcfa"
color.mix.bhlha15.xbp1<- "#cc9966"

# ------------------------------------------------------------------------------
# load data
# ------------------------------------------------------------------------------
# dataset
info <- read_tsv(paste0(prjDir, "/data/R/design.txt"),
                 col_types = "cf") %>%
  mutate(gt = factor(case_when(
    genotype == "wildtype" ~ "WT",
    str_detect(genotype, "Bhlha") ~ "Bhlha15",
    str_detect(genotype, "Xbp1") ~ "Xbp1",
    TRUE ~ ""),
    levels = c("WT", "Bhlha15", "Xbp1")))

sampleTable <- data.frame(sampleName=info$sample,
                          fileName = str_c(info$sample, ".star.uniq.sorted.bam"),
                          gt = info$gt,
                          genotype = info$genotype)

lst <- myDESeqDataSetFromFeatureCountv1_5(sampleTable = sampleTable,
                                          cfp = paste0(prjDir,"data/counts/counts.refseq.2014.txt"),
                                          design = ~ 1)
ds <- lst$ds

xbp1.binding <- read_tsv(xbp1.p2g,
                         comment = "#",
                         col_names = c("peak_chr",
                                       "peak_start",
                                       "peak_end",
                                       "gene",
                                       "gene_class",
                                       "distance"),
                         col_types = "fddcfc") # distance is a character because of unassigned peaks
bhlha15.atac_chip <- read_tsv(bhlha15.atac_chip.file,
                              comment = "#",
                              col_names = c("peak_chr",
                                            "peak_start",
                                            "peak_end",
                                            "gene",
                                            "gene_class",
                                            "distance"),
                              col_types = "fddcfc") # distance is a character because of unassigned peaks

# ------------------------------------------------------------------------------
# TPM
# ------------------------------------------------------------------------------
tpm <- counts(ds, normalized = FALSE)
tpm <- tpm/lst$geneLengths$length
tpm <- t(t(tpm)*10^6/rowSums(t(tpm)))

# alternative
tpm <- tibble(gene = rownames(tpm)) %>%
  bind_cols(as_tibble(tpm))

tpm_stat <- tpm %>%
  gather(sample, tpm, -gene) %>%
  inner_join(info %>% select(-genotype),
             by = "sample") %>%
  group_by(gt, gene) %>%
  summarize(mean = mean(tpm),
            sem = sd(tpm)/sqrt(n()))

tpm <- tpm %>%
  inner_join(tpm_stat %>%
               select(-sem) %>%
               spread(gt, mean),
             by = "gene") %>%
  inner_join(tpm_stat %>%
               select(-mean) %>%
               spread(gt, sem),
             by = "gene",
             suffix = c("", "_SEM"))

# ------------------------------------------------------------------------------
# transformations & first plots
# ------------------------------------------------------------------------------
# rlog transformations
rltVis <- rlogTransformation(ds, blind = TRUE)

# pairs plots
png(filename = paste(prjDir, "results/",opf,".pairs_rlt.png", sep = ""), width = 297*4, height = 210*4)
pairs(assay(rltVis),upper.panel=my_panel_smooth, lower.panel=my_panel_cor)
dev.off()

png(filename = paste(prjDir, "results/",opf,".pairs_log2cnts.png", sep = ""), width = 297*4, height = 210*4)
pairs(log2(counts(ds, normalize = FALSE)),upper.panel=my_panel_smooth, lower.panel=my_panel_cor)
dev.off()

pcaData <- plotPCA(rltVis, intgroup=c("genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pca <- ggplot(pcaData, aes(PC1, PC2, color = genotype, shape = genotype)) +
  geom_point(size = 2) +
  ggtitle("Bhlha15(cre/cre) vs Xbp1(F/F) Cd23-cre vs wildtype in-vitro plasmablasts\nPCA") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values = c("grey", color.bhlha15, color.xbp1)) +
  coord_fixed(ratio = 1) +
  theme_classic(base_size = 10)

pca +
  ggsave(paste0(prjDir, "results/",
                opf,".PCA.png"), width = 10, height = 10, dpi = "retina")

pca +
  ggsave(paste0(prjDir, "results/",
                opf,".PCA.pdf"), width = 10, height = 10, dpi = "retina")

# ------------------------------------------------------------------------------
# testing, group based approach
# ------------------------------------------------------------------------------
design(ds) <- formula(~ gt)
ds <- DESeq(ds, test="Wald")
padjThr <- 0.05

# dispersion
png(filename = paste0(prjDir, "results/",opf,".dispersion.png"))
plotDispEsts(ds)
dev.off()

# ------------------------------------
# Bhlha15 vs WT, results and shrinkage
# ------------------------------------
res.bhlha15 <- results(ds, name = "gt_Bhlha15_vs_WT", alpha = padjThr)

# shrinkage for potential MA plot
res.bhlha15.apeglm <- lfcShrink(ds,
                                coef="gt_Bhlha15_vs_WT",
                                res = res.bhlha15,
                                type = "apeglm")

res.bhlha15.apeglm <- bind_cols(tibble(gene = rownames(res.bhlha15.apeglm)),
                                as_tibble(res.bhlha15.apeglm))

res.bhlha15 <- tibble(gene = rownames(res.bhlha15)) %>%
  bind_cols(as_tibble(res.bhlha15)) %>%
  inner_join(tpm %>% 
               transmute(gene, 
                         TPM = case_when(
                           ((WT > 5) | (Bhlha15 > 5)) ~ "> 5",
                           TRUE ~ "<= 5"))) %>%
  inner_join(res.bhlha15.apeglm %>%
               transmute(gene,
                         shrunkenLog2FoldChange = log2FoldChange), 
             by = "gene")

# ------------------------------------
# Xbp1 vs WT, results and shrinkage
# ------------------------------------
res.xbp1 <- results(ds, name = "gt_Xbp1_vs_WT", alpha = padjThr)

# shrinkage for potential MA plot
res.xbp1.apeglm <- lfcShrink(ds,
                             coef="gt_Xbp1_vs_WT",
                             res = res.xbp1,
                             type = "apeglm")

res.xbp1.apeglm <- bind_cols(tibble(gene = rownames(res.xbp1.apeglm)),
                             as_tibble(res.xbp1.apeglm))

res.xbp1 <- tibble(gene = rownames(res.xbp1)) %>%
  bind_cols(as_tibble(res.xbp1)) %>%
  inner_join(tpm %>% 
               transmute(gene, 
                         TPM = case_when(
                           ((WT > 5) | (Xbp1 > 5)) ~ "> 5",
                           TRUE ~ "<= 5"))) %>%
  inner_join(res.xbp1.apeglm %>%
               transmute(gene,
                         shrunkenLog2FoldChange = log2FoldChange), 
             by = "gene")

# ------------------
# Xbp1 binding
# ------------------
xbp1.boundgenes <- xbp1.binding %>%
  filter(gene_class != "unassigned_peak") %>%
  group_by(gene) %>%
  summarise(peak_type = paste(unique(gene_class), collapse = ","))

res.xbp1 <- res.xbp1 %>%
  left_join(xbp1.boundgenes,
            by = "gene") %>%
  mutate(xbp1Bound = factor(if_else(is.na(peak_type), "No", "Yes")))

# ------------------
# Bhlha15 ACTAC down, ChIP binding
# ------------------
bhlha15.atac_chip_genes <- bhlha15.atac_chip  %>%
  filter(gene_class != "unassigned_peak") %>%
  group_by(gene) %>%
  summarise(peak_type = paste(unique(gene_class), collapse = ","))

res.bhlha15 <- res.bhlha15 %>%
  left_join(bhlha15.atac_chip_genes,
            by = "gene") %>%
  mutate(bhlha15.atac_chip = factor(if_else(is.na(peak_type), "No", "Yes")))
# ------------------------------------------------------------------------------
# output
# ------------------------------------------------------------------------------
# Spotfire, for ingenuity there was not enough going on to be worth the effort
rlt <- assay(rlogTransformation(ds, blind = FALSE))
rn <- rownames(rlt)

nrlog <- 2^rlt
# here, we intend to get something similar to RPMs, therefore, ALTHOUGH rlog values are already accounting for the size factor,
# we devide the genes by the total number of norm. counts per sample * 10^6 and put them on log10 scale (wished for by Meinrad)
nrlog <- t(t(nrlog*10^6)/colSums(nrlog)) # division by colSums is only possible in this way
nrlog <- log10(nrlog)

rlt <- as_tibble(rlt) %>%
  mutate(gene = rn) %>%
  rowwise() %>%
  transmute(gene = gene,
            rlog.wt      = mean(c(`124150`, `124151`)),
            rlog.bhlha15 = mean(c(`124152`, `124153`)),
            rlog.xbp1    = mean(c(`124154`, `124155`)))

nrlog <- as_tibble(nrlog) %>%
  mutate(gene = rn) %>%
  rowwise() %>%
  transmute(gene = gene,
            rRpms.wt      = mean(c(`124150`, `124151`)),
            rRpms.bhlha15 = mean(c(`124152`, `124153`)),
            rRpms.xbp1    = mean(c(`124154`, `124155`)))

spf_tbl <- res.bhlha15 %>%
  transmute(gene,
            baseMean,
            bhlha15.pvalue = pvalue,
            bhlha15.padj = padj,
            bhlha15.l2fc = log2FoldChange,
            bhlha15.sL2fc = shrunkenLog2FoldChange,
            bhlha15.stat = stat,
            bhlha15.sig = case_when(
              abs(log2FoldChange) > log2(3) & padj < 0.05 ~ "** |FC| 3",
              abs(log2FoldChange) > 1 & padj < 0.05 ~ "** |FC| 2-3",
              abs(log2FoldChange) > log2(3) & padj < 0.1 ~ "* |FC| 3",
              abs(log2FoldChange) > 1 & padj < 0.1 ~ "* |FC| 2-3",
              TRUE ~ "No"
            ),
            `bhlha15 fold change` = case_when(
              log2FoldChange > 0 ~ 2 ^ log2FoldChange,
              TRUE ~ -(2 ^ abs(log2FoldChange))),
            bhlha15.regulation = case_when(
              log2FoldChange > 0 ~ "repressed",
              log2FoldChange == 0 ~ "same",
              log2FoldChange < 0 ~ "activated",
              is.na(log2FoldChange) ~ "no counts",
              TRUE ~ "undefined"
            ),
            bhlha15.regulation.num = case_when(
              log2FoldChange < 0 ~ 1,
              log2FoldChange > 0 ~ -1,
              TRUE ~ 0
            ),
            `bhlha15 TPM filter` = TPM) %>%
  inner_join(res.xbp1 %>%
               transmute(gene,
                         xbp1.pvalue = pvalue,
                         xbp1.padj = padj,
                         xbp1.l2fc = log2FoldChange,
                         xbp1.sL2fc = shrunkenLog2FoldChange,
                         xbp1.stat = stat,
                         xbp1.sig = case_when(
                           abs(log2FoldChange) > log2(3) & padj < 0.05 ~ "** |FC| 3",
                           abs(log2FoldChange) > 1 & padj < 0.05 ~ "** |FC| 2-3",
                           abs(log2FoldChange) > log2(3) & padj < 0.1 ~ "* |FC| 3",
                           abs(log2FoldChange) > 1 & padj < 0.1 ~ "* |FC| 2-3",
                           TRUE ~ "No"
                         ),
                         `xbp1 fold change` = case_when(
                           log2FoldChange > 0 ~ 2 ^ log2FoldChange,
                           TRUE ~ -(2 ^ abs(log2FoldChange))),
                         xbp1.regulation = case_when(
                           log2FoldChange > 0 ~ "repressed",
                           log2FoldChange == 0 ~ "same",
                           log2FoldChange < 0 ~ "activated",
                           is.na(log2FoldChange) ~ "no counts",
                           TRUE ~ "undefined"
                         ),
                         xbp1.regulation.num = case_when(
                           log2FoldChange < 0 ~ 1,
                           log2FoldChange > 0 ~ -1,
                           TRUE ~ 0
                         ),
                         `xbp1 TPM filter` = TPM,
                         `xbp1 bound` = xbp1Bound),
             by = "gene") %>%
  inner_join(tpm,
             by = "gene") %>%
  inner_join(nrlog, by = "gene")

write_tsv(spf_tbl %>%
            rename(TPM_124150 = `124150`,
                   TPM_124151 = `124151`,
                   TPM_124152 = `124152`,
                   TPM_124153 = `124153`,
                   TPM_124154 = `124154`,
                   TPM_124155 = `124155`),
            na = "",
          path = str_c(prjDir, "results/", opf, ".txt"))

write_tsv(spf_tbl %>%
            filter(bhlha15.padj < 0.05 &
                     !is.na(bhlha15.padj) &
                     abs(bhlha15.l2fc) > 1 &
                     `bhlha15 TPM filter` == "> 5") %>%
            transmute(gene, bhlha15.padj, `bhlha15 fold change`,
                      bhlha15.l2fc,
                      `bhlha15 TPM filter`,
                      bhlha15.regulation,
                      bhlha15.sig,
                      TPM_WT = WT, TPM_Bhlha15KO = Bhlha15,
                      TPM_124150 = `124150`,
                      TPM_124151 = `124151`,
                      TPM_124152 = `124152`,
                      TPM_124153 = `124153`) %>%
            arrange(bhlha15.l2fc),
          na = "",
          path = str_c(prjDir, "results/", opf, ".filtered.txt"))

write_excel_csv(spf_tbl %>%
                  rename(TPM_124150 = `124150`,
                         TPM_124151 = `124151`,
                         TPM_124152 = `124152`,
                         TPM_124153 = `124153`,
                         TPM_124154 = `124154`,
                         TPM_124155 = `124155`),
                na = "",
                path = str_c(prjDir, "results/", opf, ".csv"))

write_tsv(res.bhlha15 %>%
            filter(!(is.na(padj))) %>%
            mutate(genename = str_remove(gene, "_IMP_.*")),
          str_c(prjDir, "results/bhlha15.atac_results_included.", opf, ".txt"))

write_tsv(res.bhlha15 %>%
            filter(!(is.na(padj)) &
                     bhlha15.atac_chip == "Yes") %>%
            mutate(genename = str_remove(gene, "_IMP_.*")),
          str_c(prjDir, "results/bhlha15.atac_chip_bound.", opf, ".txt"))
# ------------------------------------------------------------------------------
# plots, Bhlha15
# ------------------------------------------------------------------------------
# MA
ma.plot <- res.bhlha15 %>%
  filter(!is.na(shrunkenLog2FoldChange)) %>%
  ggplot(aes(x = baseMean, y = shrunkenLog2FoldChange, label = gene)) +
  geom_point(alpha = 0.1) +
  geom_point(data = res.bhlha15 %>%
               filter(padj < padjThr),
             col = "black", pch = 21, fill = color.bhlha15) +
  geom_hline(yintercept = c(-log2(3), -1, 1, log2(3)), lty = "21") +
  geom_text_repel(data = res.bhlha15 %>%
                    filter(padj < 0.05 &
                             abs(log2FoldChange) > log2(3) &
                             TPM == "> 5"), 
                  size = 2) +
  ggtitle("MA plot") +
  xlab("mean of normalized counts") +
  ylab("Bhlha15 dependent in-vitro plasmablasts\nshrunken log2 fold change") +
  scale_x_log10() + 
  annotation_logticks(sides = "b") +
  theme_classic(base_size = 10)

ma.plot +
  ggsave(paste0(prjDir, "results/",
                opf,".MA.png"), width = 10, height = 10, dpi = "retina")

ma.plot +
  ggsave(paste0(prjDir, "results/",
                opf,".MA.png"), width = 10, height = 10, dpi = "retina")

# Volcano
act <- res.bhlha15 %>%
  filter(log2FoldChange < -1 &
           padj < 0.05 & 
           TPM == "> 5")
rep <- res.bhlha15 %>%
  filter(log2FoldChange > 1 &
           padj < 0.05 & 
           TPM == "> 5")

volc.x.lims <- c(res.bhlha15 %>%
  filter(!is.na(padj)) %>%
    select(log2FoldChange) %>% 
    min() %>% floor(),
  res.bhlha15 %>%
    filter(!is.na(padj)) %>%
    select(log2FoldChange) %>% 
    max() %>% ceiling())

volc.plot <- res.bhlha15 %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), size = TPM, label = gene)) +
  geom_point(color = "grey", alpha = 0.4) +
  geom_point(data = rep,
             color = "red") +
  geom_text_repel(data = rep, 
                  size = 2) +
  geom_point(data = act,
             color = "blue") +
  geom_text_repel(data = act %>%
                    mutate(gene = str_remove(gene, "_IMP_.*")), 
                  size = 2) +
  geom_hline(yintercept = -log10(padjThr), lty = 2) +
  geom_vline(xintercept = c(-log2(3),-1,1,log2(3)), lty = 2) +
  xlab("log2 fold change") +
  ylab("-log10(padj)") +
  ggtitle("Bhlha15 vs wildtype") +
  annotate(geom = "text", 
           x = max(res.bhlha15$log2FoldChange, na.rm = TRUE) - 3, 
           y = -log10(min(res.bhlha15$padj, na.rm = TRUE)) + 3, 
           label = paste((rep %>% dim())[1], " repressed"),
           color = "red", size = 3) +
  annotate(geom = "text",
           x = min(res.bhlha15$log2FoldChange, na.rm = TRUE) + 0.5,
           y = -log10(min(res.bhlha15$padj, na.rm = TRUE)) + 3,
           label = paste((act %>% dim())[1], " activated"),
           color = "blue", size = 3) +
  annotate(geom = "text",
           x = c(-log2(3),log2(3)),
           y = -log10(min(res.bhlha15$padj, na.rm = TRUE)) + 3,
           label = "3x",
           color = "black", size = 3) +
  annotate(geom = "text",
           x = c(-1,1),
           y = -log10(min(res.bhlha15$padj, na.rm = TRUE)) + 3,
           label = "2x",
           color = "black", size = 3) +
  annotate(geom = "text",
           x = min(res.bhlha15$log2FoldChange, na.rm = TRUE) + 0.5,
           y = -log10(padjThr) + 1,
           label = paste("padj ", padjThr),
           color = "black", size = 3) +
  scale_size_manual(values = c(1,2)) +
  scale_x_continuous(breaks = volc.x.lims[1]:volc.x.lims[2],
                     limits = volc.x.lims) +
  theme_classic(base_size = 10)

volc.plot +
  ggsave(paste0(prjDir, "results/",
                opf,".volcano.png"), width = 10, height = 5, dpi = "retina")

volc.plot +
  ggsave(paste0(prjDir, "results/",
                opf,".volcano.pdf"), width = 10, height = 5, dpi = "retina")

limits <- c(floor(min(c(spf_tbl$rRpms.bhlha15, spf_tbl$rRpms.wt))),
            ceiling(max(c(spf_tbl$rRpms.bhlha15, spf_tbl$rRpms.wt))))
  
# scatter
scatter.plot <- spf_tbl %>%
  ggplot(aes(x = rRpms.bhlha15, y = rRpms.wt, size = -log10(bhlha15.padj))) +
  geom_point(color = "grey", alpha = 0.3, size = 1) +
  geom_point(data = spf_tbl %>%
               filter(bhlha15.l2fc > 1 &
                        bhlha15.padj < 0.05 & 
                        `bhlha15 TPM filter` == "> 5"),
             color = "red") +
  geom_point(data = spf_tbl %>%
               filter(bhlha15.l2fc < -1 &
                        bhlha15.padj < 0.05 & 
                        `bhlha15 TPM filter` == "> 5"),
             color = "blue") +
  ylab("wildtype\nExpression (nrom rlog)") +
  xlab("Bhlha15(cre/cre)\nExpression (nrom rlog)") +
  ggtitle("Bhlha15-dependence\nIn-vitro Plasmablasts") + 
  xlim(limits) + 
  ylim(limits) + 
  coord_fixed(ratio = 1) +
  annotate(geom = "text", 
           y = min(limits)+1, 
           x = max(limits)-1, 
           label = paste((rep %>% dim())[1], "repressed"),
           color = "red") +
  annotate(geom = "text",
           x = min(limits)+1,
           y = max(limits)-1,
           label = paste((act %>% dim())[1], "activated"),
           color = "blue") +
  theme_classic(base_size = 10) +
  scale_size_continuous(limits = c(-log10(0.05),-log10(min(spf_tbl$bhlha15.padj, na.rm = TRUE))), range = c(1,3))

scatter.plot +
  ggsave(paste0(prjDir, "results/",opf,".scatter.png"), width = 10, height = 10, dpi = "retina")

scatter.plot +
  ggsave(paste0(prjDir, "results/",opf,".scatter.pdf"), width = 10, height = 10, dpi = "retina")

scatter.plot.lab <- scatter.plot +
  geom_text_repel(data = spf_tbl %>%
                    filter(!is.na(bhlha15.padj) &
                             bhlha15.padj < 0.05 &
                             bhlha15.l2fc > 1 &
                             `bhlha15 TPM filter` == "> 5"),
                  aes(label = gene),
                  size = 2,
                  nudge_x = 0.25, nudge_y = -0.25,
                  segment.color = "black",
                  segment.alpha = 0.5,
                  segment.size = 0.1) +
  geom_text_repel(data = spf_tbl %>%
                    filter(!is.na(bhlha15.padj) &
                             bhlha15.padj < 0.05 &
                             bhlha15.l2fc < -1 &
                             `bhlha15 TPM filter` == "> 5") %>%
                    mutate(gene = str_remove(gene, "_IMP_.*")),
                  aes(label = gene),
                  size = 2,
                  nudge_x = -0.25, nudge_y = 0.25,
                  segment.color = "black",
                  segment.alpha = 0.5,
                  segment.size = 0.1)

scatter.plot.lab +
  ggsave(paste0(prjDir, "results/",opf,".scatter.labeled.png"), width = 10, height = 10, dpi = "retina")

scatter.plot.lab +
  ggsave(paste0(prjDir, "results/",opf,".scatter.labeled.pdf"), width = 10, height = 10, dpi = "retina")

# ------------------------------------------------------------------------------
# include ATAC&ChIP data
# ------------------------------------------------------------------------------
bhlha15.stats <- res.bhlha15 %>%
  filter(!is.na(padj)) %>%
  group_by(bhlha15.atac_chip) %>%
  summarise(mean = mean(log2FoldChange),
            median = median(log2FoldChange),
            n = n())

dens.bhlha15.atac.chip.plot <- res.bhlha15 %>%
  filter(!is.na(padj)) %>%
  ggplot(aes(x = log2FoldChange, color = bhlha15.atac_chip, fill = bhlha15.atac_chip)) +
  geom_density(alpha = 0.1) +
  geom_vline(xintercept = c(-log2(3),-1,1,log2(3)), lty = 2) +
  geom_vline(data = bhlha15.stats,
             aes(xintercept = median, color = bhlha15.atac_chip), lty = 2) +
  annotate(geom = "text", x = -8.4, y = 2.5, color = "grey", size = 3,
            label = paste("n =",
                          bhlha15.stats %>%
                            filter(bhlha15.atac_chip == "No") %>%
                            select(n))) +
  annotate(geom = "text", x = -8.5, y = 2.7, color = color.bhlha15, size = 3,
           label = paste("n =",
                         bhlha15.stats %>%
                           filter(bhlha15.atac_chip == "Yes") %>%
                           select(n))) +
  scale_color_manual(values = c("grey", color.bhlha15), name = "ATAC lost, Bhlha15 bound") +
  scale_fill_manual(values = c("grey", color.bhlha15), name = "ATAC lost, Bhlha15 bound") +
  scale_x_continuous(breaks = volc.x.lims[1]:volc.x.lims[2],
                     limits = volc.x.lims) +
  theme_classic(base_size = 10)

dens.bhlha15.atac.chip.plot +
  ggsave(paste0(prjDir, "results/",
                opf,".density.bhlha15.atac.chip.png"), width = 10, height = 2.5, dpi = "retina")

dens.bhlha15.atac.chip.plot +
  ggsave(paste0(prjDir, "results/",
                opf,".density.bhlha15.atac.chip.pdf"), width = 10, height = 2.5, dpi = "retina")

# ------------------------------------------------------------------------------
# check with E2-2 dependent genes
# ------------------------------------------------------------------------------
# THIS WAS CHANGED TO 2020_1020!
# e2.2 <- read_tsv(paste0(prjDir, "results/e2.2.de.txt"))
e2.2 <- read_tsv("/groups/busslinger/Maria/Projects/Miriam/2020_1020/results/e2.2.de.txt")

comb <- spf_tbl %>%
  select(gene, bhlha15.padj:bhlha15.stat, `bhlha15 TPM filter`,
         xbp1.padj:xbp1.stat, `xbp1 TPM filter`, `xbp1 bound`) %>%
  inner_join(e2.2 %>%
               transmute(gene, 
                         e2.2.padj = padj,
                         e2.2.l2fc = l2fc,
                         e2.2.sL2fc = sL2fc,
                         e2.2.stat = stat, 
                         `e2.2 TPM filter` = `TPM filter`),
             by = "gene",
             suffix = c("", ".e2.2")) %>%
  mutate(type = str_remove(str_c(if_else(!is.na(bhlha15.padj) & 
                                           bhlha15.padj < 0.05 &
                                           `bhlha15 TPM filter` == "> 5" &
                                           abs(bhlha15.l2fc) > 1, "Bhlha15;", ""),
                                 if_else(!is.na(xbp1.padj) & 
                                           xbp1.padj < 0.05 &
                                           `xbp1 TPM filter` == "> 5" &
                                           abs(xbp1.l2fc) > 1, "Xbp1;", ""),
                                 # Xbp1 targets (Xbp1 binding)
                                 if_else(!is.na(xbp1.padj) & 
                                           xbp1.padj < 0.05 &
                                           `xbp1 TPM filter` == "> 5" &
                                           abs(xbp1.l2fc) > 1 &
                                           `xbp1 bound` == "Yes", "|Xbp1|;", ""),
                                 if_else(!is.na(e2.2.padj) & 
                                           e2.2.padj < 0.05 &
                                           `e2.2 TPM filter` == "> 5" &
                                           abs(e2.2.l2fc) > 1, "E2-2", "")), ";$"),
         type1 = str_remove(str_remove(type, "\\|Xbp1\\|;*"), ";$"))

dens.plot <- comb %>%
  filter(str_detect(type,"E2-2|Xbp1")) %>%
  ggplot(aes(x = bhlha15.l2fc, color = type, fill = type)) +
  geom_density(data = comb %>%
                 filter(str_detect(type,"E2-2") &
                          !is.na(bhlha15.padj)) %>% # this is important because it belongs to data points in bhlha15 volcano plot!
                 mutate(type = "E2-2"),
               alpha = 0.1) +
  geom_density(data = comb %>%
                 filter(str_detect(type,"Xbp1") &
                          !is.na(bhlha15.padj)) %>% # this is important because it belongs to data points in bhlha15 volcano plot!
                 mutate(type = "Xbp1"),
               alpha = 0.4) +
  geom_density(data = comb %>%
                 filter(str_detect(type,"\\|Xbp1\\|") &
                          !is.na(bhlha15.padj)) %>% # this is important because it belongs to data points in bhlha15 volcano plot!
                 mutate(type = "Xbp1 bound"),
               alpha = 0.4) +
  geom_vline(xintercept = c(-log2(3),-1,1,log2(3)), lty = 2) +
  annotate(geom = "text", x = -8.5, y = 1.8, color = color.e2.2, size = 3,
           label = paste("n =",
                         (comb %>%
                           filter(str_detect(type,"E2-2") &
                                    !is.na(bhlha15.padj)) %>% dim())[1])) +
  annotate(geom = "text", x = -8.5, y = 2.2, color = color.xbp1.unbound, size = 3,
           label = paste("n =",
                         (comb %>%
                            filter(str_detect(type,"Xbp1") &
                                     !is.na(bhlha15.padj)) %>% dim())[1])) +
  annotate(geom = "text", x = -8.5, y = 2, color = color.xbp1, size = 3,
           label = paste("n =",
                         (comb %>%
                            filter(str_detect(type,"Xbp1") &
                                     !is.na(bhlha15.padj) &
                                     `xbp1 bound` == "Yes") %>% dim())[1])) +
  scale_color_manual(values = rep("black",3), name = "Deregulated in") +
  scale_fill_manual(values = c(color.e2.2, color.xbp1.unbound, color.xbp1), name = "Deregulated in") +
  scale_x_continuous(breaks = volc.x.lims[1]:volc.x.lims[2],
                     limits = volc.x.lims) +
  theme_classic(base_size = 10)

dens.plot +
  ggsave(paste0(prjDir, "results/",
                opf,".density.e2-2_xbp1.png"), width = 10, height = 2.5, dpi = "retina")

dens.plot +
  ggsave(paste0(prjDir, "results/",
                opf,".density.e2-2_xbp1.pdf"), width = 10, height = 2.5, dpi = "retina")

sl2fc.limit <- max(comb %>% select(contains("sL2fc")) %>% min(na.rm = TRUE) %>% abs(),
                   comb %>% select(contains("sL2fc")) %>% max(na.rm = TRUE)) %>% round()

bhlha15.sl2fc_vs_e2.2sl2fc <- comb %>%
  filter(!is.na(bhlha15.sL2fc) & !is.na(e2.2.sL2fc)) %>%
  mutate(type1 = str_remove(str_remove(type1, "Xbp1"),
                           "^;|;$")) %>%
  ggplot(aes(x = bhlha15.sL2fc, y = e2.2.sL2fc, fill = type1, label = gene)) +
  geom_point(alpha = 0.4, pch = 21, color = "black") +
  geom_hline(yintercept = c(-1,1), lty = 2) +
  geom_vline(xintercept = c(-1,1), lty = 2) +
  geom_point(data = comb %>%
               filter(str_detect(type, "Bhlha15|E2-2")) %>%
               mutate(type1 = str_remove(str_remove(type1, "Xbp1"),
                                        "^;|;$")),
             pch = 21, color = "black") +
  geom_text_repel(data = comb %>%
                    filter(str_detect(type, "Bhlha15|E2-2")) %>%
                    mutate(gene = str_remove(gene, "_IMP_.*"),
                           type1 = str_remove(str_remove(type1, "Xbp1"),
                                             "^;|;$")),
                  size = 2,
                  segment.size = 0.1) +
  ggtitle("shrunken log2 fold change vs shrunken log2 fold change\napeglm") +
  xlab("Bhlha15 dependence\n(shrunken log2 fold change)") +
  ylab("E2-2 dependence\n(shrunken log2 fold change)") +
  xlim(c(-sl2fc.limit, sl2fc.limit)) +
  ylim(c(-sl2fc.limit, sl2fc.limit)) +
  coord_fixed(ratio = 1) +
  scale_fill_manual(values = c("grey", color.bhlha15, color.e2.2),
                    name = "Deregulated in") +
  theme_classic(base_size = 10)

bhlha15.sl2fc_vs_e2.2sl2fc +
  ggsave(paste0(prjDir, "results/bhlha15.sl2fc_vs_e2.2.sl2fc.png"), width = 10, height = 10, dpi = "retina")

bhlha15.sl2fc_vs_e2.2sl2fc +
  ggsave(paste0(prjDir, "results/bhlha15.sl2fc_vs_e2.2.sl2fc.pdf"), width = 10, height = 10, dpi = "retina")

# here, we ignore Xbp1 binding
comb1 <- comb %>%
  filter(!is.na(bhlha15.sL2fc) & !is.na(xbp1.sL2fc)) %>%
  mutate(type = factor(str_replace(
    str_remove(str_remove(type1, "E2-2"),
               "^;|;$"),
    ";"," & "),
    levels = c("", "Bhlha15", "Xbp1", "Bhlha15 & Xbp1")))

bhlha15.sl2fc_vs_xbp1.sl2fc <- comb1 %>%
  ggplot(aes(x = bhlha15.sL2fc, y = xbp1.sL2fc, fill = type, label = gene)) +
  geom_point(alpha = 0.4, pch = 21, color = "black") +
  geom_hline(yintercept = c(-1,1), lty = 2) +
  geom_vline(xintercept = c(-1,1), lty = 2) +
  geom_point(data = comb1 %>%
               filter(type != ""),
             pch = 21, color = "black") +
    geom_text_repel(data = comb1 %>%
                    filter(type == "Bhlha15 & Xbp1") %>%
                    mutate(gene = str_remove(gene, "_IMP_.*")),
                  size = 2,
                  segment.size = 0.1) +
  ggtitle("shrunken log2 fold change vs shrunken log2 fold change\napeglm") +
  xlab("Bhlha15 dependence\n(shrunken log2 fold change)") +
  ylab("Xbp1 dependence\n(shrunken log2 fold change)") +
  xlim(c(-sl2fc.limit, sl2fc.limit)) +
  ylim(c(-sl2fc.limit, sl2fc.limit)) +
  coord_fixed(ratio = 1) +
  scale_fill_manual(values = c("grey", color.bhlha15, color.xbp1, color.mix.bhlha15.xbp1),
                    name = "Deregulated in") +
  theme_classic(base_size = 10)

bhlha15.sl2fc_vs_xbp1.sl2fc +
  ggsave(paste0(prjDir, "results/bhlha15.sl2fc_vs_xbp1.sl2fc.png"), width = 10, height = 10, dpi = "retina")

bhlha15.sl2fc_vs_xbp1.sl2fc +
  ggsave(paste0(prjDir, "results/bhlha15.sl2fc_vs_xbp1.sl2fc.pdf"), width = 10, height = 10, dpi = "retina")

comb2 <- comb %>%
  filter(!is.na(e2.2.sL2fc) & !is.na(xbp1.sL2fc)) %>%
  mutate(type = factor(str_replace(
    str_remove(str_remove(type1, "Bhlha15"),
               "^;|;$"),
    ";"," & "),
    levels = c("", "Xbp1", "E2-2", "Xbp1 & E2-2")))

e2.2.sl2fc_vs_xbp1.sl2fc <- comb2 %>%
  ggplot(aes(x = e2.2.sL2fc, y = xbp1.sL2fc, fill = type, label = gene)) +
  geom_point(alpha = 0.4, pch = 21, color = "black") +
  geom_hline(yintercept = c(-1,1), lty = 2) +
  geom_vline(xintercept = c(-1,1), lty = 2) +
  geom_point(data = comb2 %>%
               filter(type != ""),
             pch = 21, color = "black") +
  geom_text_repel(data = comb2 %>%
                    filter(type == "Xbp1 & E2-2") %>%
                    mutate(gene = str_remove(gene, "_IMP_.*")),
                  size = 2,
                  segment.size = 0.1) +
  ggtitle("shrunken log2 fold change vs shrunken log2 fold change\napeglm") +
  xlab("E2-2 dependence\n(shrunken log2 fold change)") +
  ylab("Xbp1 dependence\n(shrunken log2 fold change)") +
  xlim(c(-sl2fc.limit, sl2fc.limit)) +
  ylim(c(-sl2fc.limit, sl2fc.limit)) +
  coord_fixed(ratio = 1) +
  scale_fill_manual(values = c("grey", color.xbp1, color.e2.2, "#d06278"),
                    name = "Deregulated in") +
  theme_classic(base_size = 10)

e2.2.sl2fc_vs_xbp1.sl2fc +
  ggsave(paste0(prjDir, "results/e2.2.sl2fc_vs_xbp1.sl2fc.png"), width = 10, height = 10, dpi = "retina")

e2.2.sl2fc_vs_xbp1.sl2fc +
  ggsave(paste0(prjDir, "results/e2.2.sl2fc_vs_xbp1.sl2fc.pdf"), width = 10, height = 10, dpi = "retina")

# ------------------------------------------------------------------------------
# plots, Xbp1
# ------------------------------------------------------------------------------
# MA
xbp1.ma.plot <- res.xbp1 %>%
  filter(!is.na(shrunkenLog2FoldChange)) %>%
  ggplot(aes(x = baseMean, y = shrunkenLog2FoldChange, label = gene, pch = xbp1Bound, fill = xbp1Bound)) +
  geom_point(alpha = 0.1, col = "black", fill = "black") +
  geom_point(data = res.xbp1 %>%
               filter(padj < padjThr),
             col = "black") +
  geom_hline(yintercept = c(-log2(3), -1, 1, log2(3)), lty = "21") +
  # highlight top 20 genes
  geom_text_repel(data = res.xbp1 %>%
                    filter(padj < 0.05 &
                             abs(log2FoldChange) > log2(3) &
                             TPM == "> 5") %>%
                    arrange(desc(abs(log2FoldChange))) %>%
                    filter(row_number() < 21) %>%
                    mutate(gene = str_remove(gene, "_IMP.*")), 
                  size = 2) +
  ggtitle("MA plot") +
  xlab("mean of normalized counts") +
  ylab("Xbp1 dependent in-vitro plasmablasts\nshrunken log2 fold change") +
  scale_x_log10() +
  annotation_logticks(sides = "b") +
  scale_shape_manual(values = c(21,23), name = "Xbp1 binding") +
  scale_fill_manual(values = c(color.xbp1.unbound, color.xbp1), name = "Xbp1 binding") +
  theme_classic(base_size = 10)

xbp1.ma.plot +
  ggsave(paste0(prjDir, "results/xbp1.grp.de.MA.png"), width = 10, height = 10, dpi = "retina")

xbp1.ma.plot +
  ggsave(paste0(prjDir, "results/xbp1.grp.de.MA.pdf"), width = 10, height = 10, dpi = "retina")

# Volcano
xbp1.act <- res.xbp1 %>%
  filter(log2FoldChange < -1 &
           padj < 0.05 & 
           TPM == "> 5")
xbp1.rep <- res.xbp1 %>%
  filter(log2FoldChange > 1 &
           padj < 0.05 & 
           TPM == "> 5")

xbp1.volc.x.lims <- c(res.xbp1 %>%
                        filter(!is.na(padj)) %>%
                        select(log2FoldChange) %>% 
                        min() %>% floor(),
                      res.xbp1 %>%
                        filter(!is.na(padj)) %>%
                        select(log2FoldChange) %>% 
                        max() %>% ceiling())

xbp1.volc.plot <- res.xbp1 %>%
  filter(!is.na(padj) & padj != 0) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), size = TPM, label = gene, shape = xbp1Bound)) +
  geom_point(color = "grey", fill = "grey", alpha = 0.4) +
  geom_point(data = xbp1.act %>%
               filter(padj == 0),
             fill = "lightblue", color = "black") +
  geom_point(data = xbp1.rep %>%
               filter(xbp1Bound == "No"),
             color = "red") +
  geom_point(data = xbp1.rep %>%
               filter(xbp1Bound == "Yes"),
             fill = "red", color = "black") +
  geom_text_repel(data = bind_rows(res.xbp1 %>%
                                     filter(padj < 0.05 &
                                              log2FoldChange > log2(3) &
                                              TPM == "> 5") %>%
                                     arrange(desc(abs(log2FoldChange))) %>%
                                     mutate(rank.l2fc = row_number(),
                                            gene = str_remove(gene, "_IMP.*")) %>%
                                     arrange(padj) %>%
                                     mutate(rank.padj = row_number()) %>%
                                     filter(rank.l2fc < 11 | rank.padj < 11),
                                   res.xbp1 %>%
                                     filter(padj < 0.05 &
                                              log2FoldChange < -log2(3) &
                                              TPM == "> 5") %>%
                                     arrange(log2FoldChange) %>%
                                     mutate(rank.l2fc = row_number(),
                                            gene = str_remove(gene, "_IMP.*")) %>%
                                     arrange(padj) %>%
                                     mutate(rank.padj = row_number()) %>%
                                     filter(rank.l2fc < 11 | rank.padj < 11)), 
                  size = 2) +
  geom_point(data = xbp1.act %>%
               filter(padj != 0 & 
                        xbp1Bound == "No"),
             color = "blue") +
  geom_point(data = xbp1.act %>%
               filter(padj != 0 & 
                        xbp1Bound == "Yes"),
             fill = "blue", color = "black") +
  geom_hline(yintercept = -log10(padjThr), lty = 2) +
  geom_vline(xintercept = c(-log2(3),-1,1,log2(3)), lty = 2) +
  xlab("log2 fold change") +
  ylab("-log10(padj)") +
  ylim(c(0,400)) +
  ggtitle("Xbp1 vs wildtype",
          subtitle = "lightblue dots indicate VERY low adjust p-values (digit precision issue in R, wrongly shown as 0)") +
  annotate(geom = "text", 
           x = max(res.xbp1$log2FoldChange, na.rm = TRUE) - 3, 
           y = 375, 
           label = paste0((xbp1.rep %>% dim())[1], " repressed\n(",
                         (xbp1.rep %>% 
                            filter(xbp1Bound == "Yes") %>%
                            dim())[1], " bound)"),
           color = "red", size = 3) +
  annotate(geom = "text",
           x = min(res.xbp1$log2FoldChange, na.rm = TRUE) + 0.5,
           y = 375,
           label = paste0((xbp1.act %>% dim())[1], " activated\n(",
                          (xbp1.act %>% 
                             filter(xbp1Bound == "Yes") %>%
                             dim())[1], " bound)"),
           color = "blue", size = 3) +
  annotate(geom = "text",
           x = c(-log2(3),log2(3)), y = 400,
           label = "3x",
           color = "black", size = 3) +
  annotate(geom = "text",
           x = c(-1,1), y = 400,
           label = "2x",
           color = "black", size = 3) +
  annotate(geom = "text",
           x = min(res.xbp1$log2FoldChange, na.rm = TRUE) + 0.5,
           y = -log10(padjThr) + 1,
           label = paste("padj ", padjThr),
           color = "black", size = 3) +
  scale_size_manual(values = c(1,2)) +
  scale_shape_manual(values = c(19,23)) +
  scale_x_continuous(breaks = xbp1.volc.x.lims[1]:xbp1.volc.x.lims[2],
                     limits = xbp1.volc.x.lims) +
  theme_classic(base_size = 10)

xbp1.volc.plot +
  ggsave(paste0(prjDir, "results/xbp1.grp.de.volcano.png"), width = 10, height = 5, dpi = "retina")

xbp1.volc.plot +
  ggsave(paste0(prjDir, "results/xbp1.grp.de.volcano.pdf"), width = 10, height = 5, dpi = "retina")

xbp1.limits <- c(floor(min(c(spf_tbl$rRpms.xbp1, spf_tbl$rRpms.wt))),
                 ceiling(max(c(spf_tbl$rRpms.xbp1, spf_tbl$rRpms.wt))))

# scatter
xbp1.scatter.plot <- spf_tbl %>%
  ggplot(aes(x = rRpms.xbp1, y = rRpms.wt, size = -log10(xbp1.padj), shape = `xbp1 bound`)) +
  geom_point(color = "grey", fill = "grey", alpha = 0.3, size = 1) +
  geom_point(data = spf_tbl %>%
               filter(xbp1.l2fc > 1 &
                        xbp1.padj < 0.05 & 
                        `xbp1 TPM filter` == "> 5" &
                        `xbp1 bound` == "No"),
             color = "red") +
  geom_point(data = spf_tbl %>%
               filter(xbp1.l2fc > 1 &
                        xbp1.padj < 0.05 & 
                        `xbp1 TPM filter` == "> 5" &
                        `xbp1 bound` == "Yes"),
             fill = "red", color = "black") +
  geom_point(data = spf_tbl %>%
               filter(xbp1.l2fc < -1 &
                        xbp1.padj < 0.05 & 
                        `xbp1 TPM filter` == "> 5" &
                        `xbp1 bound` == "No"),
             color = "blue") +
  geom_point(data = spf_tbl %>%
               filter(xbp1.l2fc < -1 &
                        xbp1.padj < 0.05 & 
                        `xbp1 TPM filter` == "> 5" &
                        `xbp1 bound` == "Yes"),
             fill = "blue", color = "black") +
  ylab("wildtype\nExpression (nrom rlog)") +
  xlab("Xbp1(F/F) Cd23-cre\nExpression (nrom rlog)") +
  ggtitle("Xbp1-dependence\nIn-vitro Plasmablasts") +
  xlim(xbp1.limits) + 
  ylim(xbp1.limits) + 
  coord_fixed(ratio = 1) +
  annotate(geom = "text", 
           y = min(xbp1.limits)+1, 
           x = max(xbp1.limits)-1, 
           label = paste((xbp1.rep %>% dim())[1], " repressed\n(",
                         (xbp1.rep %>% 
                            filter(xbp1Bound == "Yes") %>%
                            dim())[1], " bound)"),
           color = "red") +
  annotate(geom = "text",
           x = min(xbp1.limits)+1,
           y = max(xbp1.limits)-1,
           label = paste((xbp1.act %>% dim())[1], " activated\n(",
                         (xbp1.act %>% 
                            filter(xbp1Bound == "Yes") %>%
                            dim())[1], " bound)"),
           color = "blue") +
  theme_classic(base_size = 10) +
  scale_size_continuous(limits = c(-log10(0.05),
                                   -log10(spf_tbl %>% filter(xbp1.padj != 0) %>% select(xbp1.padj) %>% min(na.rm = TRUE))),
                        range = c(1,3)) +
  scale_shape_manual(values = c(19,23))

xbp1.scatter.plot +
  ggsave(paste0(prjDir, "results/xbp1.grp.de.scatter.png"), width = 10, height = 10, dpi = "retina")

xbp1.scatter.plot +
  ggsave(paste0(prjDir, "results/xbp1.grp.de.scatter.pdf"), width = 10, height = 10, dpi = "retina")

xbp1.scatter.lab.plot <- spf_tbl %>%
  ggplot(aes(x = rRpms.xbp1, y = rRpms.wt, size = -log10(xbp1.padj))) +
  geom_point(color = "grey", fill = "grey", alpha = 0.3, size = 1) +
  geom_point(data = spf_tbl %>%
               filter(xbp1.l2fc > 1 &
                        !is.na(xbp1.padj) &
                        xbp1.padj < 0.05 & 
                        `xbp1 TPM filter` == "> 5"),
             color = "red") +
  geom_point(data = spf_tbl %>%
               filter(xbp1.l2fc < -1 &
                        !is.na(xbp1.padj) &
                        xbp1.padj < 0.05 & 
                        `xbp1 TPM filter` == "> 5") %>%
               mutate(xbp1.padj = if_else(xbp1.padj == 0, 10^-302, xbp1.padj)), # this needs to be done due to rounding errors otherwise
             color = "blue") +
  geom_text_repel(data = spf_tbl %>%
                    filter(!is.na(xbp1.padj) &
                             xbp1.padj < 0.05 &
                             xbp1.l2fc > 1 &
                             `xbp1 TPM filter` == "> 5") %>%
                    arrange(desc(xbp1.l2fc)) %>%
                    filter(row_number() < 16) %>%
                    mutate(gene = str_remove(gene, "_IMP_.*")),
                  aes(label = gene),
                  size = 2,
                  nudge_x = 0.25, nudge_y = -0.25,
                  segment.color = "black",
                  segment.alpha = 0.5,
                  segment.size = 0.1) +
  geom_text_repel(data = spf_tbl %>%
                    filter(!is.na(xbp1.padj) &
                             xbp1.padj < 0.05 &
                             xbp1.l2fc < -1 &
                             `xbp1 TPM filter` == "> 5") %>%
                    arrange(xbp1.l2fc) %>%
                    filter(row_number() < 16) %>%
                    mutate(gene = str_remove(gene, "_IMP_.*")),
                  aes(label = gene),
                  size = 2,
                  nudge_x = -0.25, nudge_y = 0.25,
                  segment.color = "black",
                  segment.alpha = 0.5,
                  segment.size = 0.1) +
  ylab("wildtype\nExpression (nrom rlog)") +
  xlab("Xbp1(F/F) Cd23-cre\nExpression (nrom rlog)") +
  ggtitle("Xbp1-dependence\nIn-vitro Plasmablasts") +
  xlim(xbp1.limits) + 
  ylim(xbp1.limits) + 
  coord_fixed(ratio = 1) +
  annotate(geom = "text", 
           y = min(xbp1.limits)+1, 
           x = max(xbp1.limits)-1, 
           label = paste((xbp1.rep %>% dim())[1], " repressed\n(",
                         (xbp1.rep %>% 
                            filter(xbp1Bound == "Yes") %>%
                            dim())[1], " bound)"),
           color = "red") +
  annotate(geom = "text",
           x = min(xbp1.limits)+1,
           y = max(xbp1.limits)-1,
           label = paste((xbp1.act %>% dim())[1], " activated\n(",
                         (xbp1.act %>% 
                            filter(xbp1Bound == "Yes") %>%
                            dim())[1], " bound)"),
           color = "blue") +
  theme_classic(base_size = 10) +
  scale_size_continuous(limits = c(-log10(0.05),
                                   -log10(spf_tbl %>% 
                                            mutate(xbp1.padj = if_else(xbp1.padj == 0, 10^-302, xbp1.padj)) %>%
                                            select(xbp1.padj) %>% min(na.rm = TRUE))),
                        range = c(1,3))

xbp1.scatter.lab.plot +
  ggsave(paste0(prjDir, "results/xbp1.grp.de.scatter.labeled.png"), width = 10, height = 10, dpi = "retina")

xbp1.scatter.lab.plot +
  ggsave(paste0(prjDir, "results/xbp1.grp.de.scatter.labeled.pdf"), width = 10, height = 10, dpi = "retina")

xbp1.scatter.lab3.plot <- spf_tbl %>%
  ggplot(aes(x = rRpms.xbp1, y = rRpms.wt, size = -log10(xbp1.padj))) +
  geom_point(color = "grey", fill = "grey", alpha = 0.3, size = 1) +
  geom_point(data = spf_tbl %>%
               filter(xbp1.l2fc > log2(3) &
                        !is.na(xbp1.padj) &
                        xbp1.padj < 0.05 & 
                        `xbp1 TPM filter` == "> 5"),
             color = "red") +
  geom_point(data = spf_tbl %>%
               filter(xbp1.l2fc < -log2(3) &
                        !is.na(xbp1.padj) &
                        xbp1.padj < 0.05 & 
                        `xbp1 TPM filter` == "> 5") %>%
               mutate(xbp1.padj = if_else(xbp1.padj == 0, 10^-302, xbp1.padj)), # this needs to be done due to rounding errors otherwise
             color = "blue") +
  geom_text_repel(data = spf_tbl %>%
                    filter(!is.na(xbp1.padj) &
                             xbp1.padj < 0.05 &
                             xbp1.l2fc > log2(3) &
                             `xbp1 TPM filter` == "> 5") %>%
                    arrange(desc(xbp1.l2fc)) %>%
                    filter(row_number() < 16) %>%
                    mutate(gene = str_remove(gene, "_IMP_.*")),
                  aes(label = gene),
                  size = 2,
                  nudge_x = 0.25, nudge_y = -0.25,
                  segment.color = "black",
                  segment.alpha = 0.5,
                  segment.size = 0.1) +
  geom_text_repel(data = spf_tbl %>%
                    filter(!is.na(xbp1.padj) &
                             xbp1.padj < 0.05 &
                             xbp1.l2fc < -log2(3) &
                             `xbp1 TPM filter` == "> 5") %>%
                    arrange(xbp1.l2fc) %>%
                    filter(row_number() < 16) %>%
                    mutate(gene = str_remove(gene, "_IMP_.*")),
                  aes(label = gene),
                  size = 2,
                  nudge_x = -0.25, nudge_y = 0.25,
                  segment.color = "black",
                  segment.alpha = 0.5,
                  segment.size = 0.1) +
  ylab("wildtype\nExpression (nrom rlog)") +
  xlab("Xbp1(F/F) Cd23-cre\nExpression (nrom rlog)") +
  ggtitle("Xbp1-dependence\nIn-vitro Plasmablasts") +
  xlim(xbp1.limits) + 
  ylim(xbp1.limits) + 
  coord_fixed(ratio = 1) +
  annotate(geom = "text", 
           y = min(xbp1.limits)+1, 
           x = max(xbp1.limits)-1, 
           label = paste((xbp1.rep %>% 
                            filter(log2FoldChange > log2(3)) %>%
                            dim())[1], " repressed\n(",
                         (xbp1.rep %>% 
                            filter(log2FoldChange > log2(3) &
                                     xbp1Bound == "Yes") %>%
                            dim())[1], " bound)"),
           color = "red") +
  annotate(geom = "text",
           x = min(xbp1.limits)+1,
           y = max(xbp1.limits)-1,
           label = paste((xbp1.act %>%
                            filter(log2FoldChange < -log2(3)) %>%
                            dim())[1], " activated\n(",
                         (xbp1.act %>% 
                            filter(log2FoldChange < -log2(3) &
                                     xbp1Bound == "Yes") %>%
                            dim())[1], " bound)"),
           color = "blue") +
  theme_classic(base_size = 10) +
  scale_size_continuous(limits = c(-log10(0.05),
                                   -log10(spf_tbl %>% 
                                            mutate(xbp1.padj = if_else(xbp1.padj == 0, 10^-302, xbp1.padj)) %>%
                                            select(xbp1.padj) %>% min(na.rm = TRUE))),
                        range = c(1,3))

xbp1.scatter.lab3.plot +
  ggsave(paste0(prjDir, "results/xbp1.grp.de.scatter.3fold_labeled.png"), width = 10, height = 10, dpi = "retina")

xbp1.scatter.lab3.plot +
  ggsave(paste0(prjDir, "results/xbp1.grp.de.scatter.3fold_labeled.pdf"), width = 10, height = 10, dpi = "retina")
# ------------------------------------------------------------------------------
# GSEA
# ------------------------------------------------------------------------------
# grp for GSEA Broad program (pobably no longer needed)
write_tsv(bind_rows(tibble(gene = c("mm9")),
                    spf_tbl %>%
                      filter(`bhlha15 TPM filter` == "> 5" &
                               bhlha15.padj < 0.05 &
                               bhlha15.l2fc > 1) %>%
                      select(gene)),
          path = str_c(prjDir, "results/gsea/bhlha15_rep.grp"),
          col_names = FALSE)

write_tsv(bind_rows(tibble(gene = c("mm9")),
                    spf_tbl %>%
                      filter(`bhlha15 TPM filter` == "> 5" &
                               bhlha15.padj < 0.05 &
                               bhlha15.l2fc < -1) %>%
                      select(gene)),
          path = str_c(prjDir, "results/gsea/bhlha15_act.grp"),
          col_names = FALSE)

write_tsv(bind_rows(tibble(gene = c("mm9")),
                    bhlha15.atac_chip_genes %>%
                      select(gene)),
          path = str_c(prjDir, "results/gsea/bhlha15_atac_chip.grp"),
          col_names = FALSE)

# rnk for GSEA
write_tsv(spf_tbl %>%
            transmute(`#gene` = gene, xbp1.sL2fc) %>%
            filter(!is.na(xbp1.sL2fc)),
          path = str_c(prjDir, "results/gsea/xbp1.rnk"))

write_tsv(spf_tbl %>%
            transmute(`#gene` = gene, bhlha15.sL2fc) %>%
            filter(!is.na(bhlha15.sL2fc)),
          path = str_c(prjDir, "results/gsea/bhlha15.rnk"))

##########################################################
# Bhlha15 ChIP ATAC vs Bhlha15 RNA-seq
##########################################################
geneList.bhlha15 <- res.bhlha15 %>%
  filter(!is.na(shrunkenLog2FoldChange)) %>%
  arrange(desc(shrunkenLog2FoldChange)) %>%
  pull(shrunkenLog2FoldChange)

names(geneList.bhlha15) <- res.bhlha15 %>%
  filter(!is.na(shrunkenLog2FoldChange)) %>%
  arrange(desc(shrunkenLog2FoldChange)) %>%
  pull(gene)

gsid2gene <- data.frame(term = rep("BHLHA15_ATAC_CHIP",
                                   bhlha15.atac_chip_genes %>%
                                     pull(gene) %>% length),
                        gene = bhlha15.atac_chip_genes %>%
                          pull(gene))

gsid2name <- data.frame(term = rep("BHLHA15_ATAC_CHIP",
                                   bhlha15.atac_chip_genes %>%
                                     pull(gene) %>% length),
                        name = rep("Bhlha15 ATAC lost, ChIP peak",
                                   bhlha15.atac_chip_genes %>%
                                     pull(gene) %>% length))

gsea.bhlha15ChIP <- GSEA(geneList.bhlha15, TERM2GENE = gsid2gene, TERM2NAME = gsid2name)

gseaplot2.bhlha15 <- gseaplot2(gsea.bhlha15ChIP, geneSetID = 1, 
                               title = str_c(gsea.bhlha15ChIP$Description[1],
                                             "\nNES: ", gsea.bhlha15ChIP$NES[1]), 
                               pvalue_table = TRUE, 
                               base_size = 10,
                               color = color.bhlha15,
                               subplots = 1:2,
                               rel_heights = c(8, 1))

gseaplot2.bhlha15 +
  ggsave(paste0(prjDir, "results/gsea/bhlha15.gsea.Rplot.png"), width = 10, height = 5, dpi = "retina")
gseaplot2.bhlha15 +
  ggsave(paste0(prjDir, "results/gsea/bhlha15.gsea.Rplot.pdf"), width = 10, height = 5, dpi = "retina")

gseaplot.bhlha15 <- gseaplot(gsea.bhlha15ChIP, 
                             "BHLHA15_ATAC_CHIP", 
                             color.line = color.bhlha15, 
                             color = "darkgrey",
                             color.vline = "darkgrey",
                             by = "runningScore",
                             base_size = 4,
                             title = str_c("Set: ", gsea.bhlha15ChIP$Description[1], 
                                           "; ranking: Bhlha15 RNA-seq (sL2fc)")) +
  annotate(geom = "text", 
           x = c(1000, 18500), y = c(0.1, 0.1), 
           label = c("Bhlha15 repression",
                     "Bhlha15 activation"),
           color = c("red", "blue"), size = 3) +
  annotate(geom = "text", 
           x = 3000, y = -0.6, 
           label = str_c("NES: ", format(gsea.bhlha15ChIP$NES[1], digits = 3),
                     "\npadj: ", format(gsea.bhlha15ChIP$p.adjust[1], digits = 3)),
           color = "black", size = 4) +
  annotate(geom = "text", 
           x = length(geneList.bhlha15) - gsea.bhlha15ChIP$rank[1], y = -0.1, 
           label = format(length(geneList.bhlha15) - gsea.bhlha15ChIP$rank[1], big.mark = ","),
           color = "black", size = 3)

gseaplot.bhlha15 +
  theme_classic(base_size = 10) +
  ggsave(paste0(prjDir, "results/gsea/bhlha15.gsea.Rplot.v2.png"), width = 10, height = 5, dpi = "retina")
gseaplot.bhlha15 +
  theme_classic(base_size = 10) +
  ggsave(paste0(prjDir, "results/gsea/bhlha15.gsea.Rplot.v2.pdf"), width = 10, height = 5, dpi = "retina")

as_tibble(gsea.bhlha15ChIP) %>%
  mutate(core_enrichment = str_replace_all(core_enrichment, "/", "\n")) %>%
  write_excel_csv(str_c(prjDir, "results/gsea/bhlha15.gsea.csv"))

##########################################################
# Bhlha15 RNA-seq
##########################################################
geneList.xbp1 <- res.xbp1 %>%
  filter(!is.na(shrunkenLog2FoldChange)) %>%
  arrange(desc(shrunkenLog2FoldChange)) %>%
  pull(shrunkenLog2FoldChange)

names(geneList.xbp1) <- res.xbp1 %>%
  filter(!is.na(shrunkenLog2FoldChange)) %>%
  arrange(desc(shrunkenLog2FoldChange)) %>%
  pull(gene)

bhlha15.gsid2gene <- data.frame(term = c(rep("BHLHA15_ATAC_CHIP",
                                             bhlha15.atac_chip_genes %>%
                                               pull(gene) %>% length),
                                         rep("BHLHA15_ATCTIVATED_RNASEQ",
                                             (act %>% dim())[1])),
                                gene = c(bhlha15.atac_chip_genes %>%
                                           pull(gene),
                                         act %>% pull(gene)))

bhlha15.gsid2name <- data.frame(term = c(rep("BHLHA15_ATAC_CHIP",
                                             bhlha15.atac_chip_genes %>%
                                               pull(gene) %>% length),
                                         rep("BHLHA15_ATCTIVATED_RNASEQ",
                                             (act %>% dim())[1])),
                                name = c(rep("Bhlha15 ATAC lost, ChIP peak",
                                             bhlha15.atac_chip_genes %>%
                                               pull(gene) %>% length),
                                         rep("Bhlha15 activated genes",
                                             (act %>% dim())[1])))

gsea.xbp1 <- GSEA(geneList.xbp1, 
                  TERM2GENE = bhlha15.gsid2gene, TERM2NAME = bhlha15.gsid2name)

t.gsea.xbp1 <- as_tibble(gsea.xbp1) %>%
  mutate(core_enrichment = str_replace_all(core_enrichment, "/", "\n"))

gseaplot.xbp1_vs_bhlha15act <- gseaplot(gsea.xbp1, 
                          "BHLHA15_ATCTIVATED_RNASEQ", 
                          color.line = color.xbp1, 
                          color = "darkgrey",
                          color.vline = "darkgrey",
                          by = "runningScore",
                          base_size = 4,
                          title = str_c("Set: ", gsea.xbp1$Description[2],
                                        "; ranking: Xbp1 RNA-seq (sL2fc)")) +
  annotate(geom = "text", 
           x = c(1000, 18500), y = c(0.1, 0.1), 
           label = c("Xbp1 repression",
                     "Xbp1 activation"),
           color = c("red", "blue"), size = 3) +
  annotate(geom = "text", 
           x = 3000, y = -0.6, 
           label = str_c("NES: ", format(gsea.xbp1$NES[2], digits = 3),
                         "\npadj: ", format(gsea.xbp1$p.adjust[2], digits = 3)),
           color = "black", size = 4)  +
  annotate(geom = "text", 
           x = length(geneList.xbp1) - gsea.xbp1$rank[2], y = -0.1, 
           label = format(length(geneList.xbp1) - gsea.xbp1$rank[2], big.mark = ","),
           color = "black", size = 3)

gseaplot.xbp1_vs_bhlha15act +
  theme_classic(base_size = 10) +
  ggsave(paste0(prjDir, "results/gsea/xbp1.gsea.Rplot.png"), width = 10, height = 5, dpi = "retina")

gseaplot.xbp1_vs_bhlha15act +
  theme_classic(base_size = 10) +
  ggsave(paste0(prjDir, "results/gsea/xbp1.gsea.Rplot.pdf"), width = 10, height = 5, dpi = "retina")

gseaplot.xbp1_vs_bhlha15ATAC <- gseaplot(gsea.xbp1, 
                                        "BHLHA15_ATAC_CHIP", 
                                        color.line = color.xbp1, 
                                        color = "darkgrey",
                                        color.vline = "darkgrey",
                                        by = "runningScore",
                                        base_size = 4,
                                        title = str_c("Set: ", gsea.xbp1$Description[1],
                                                      "; ranking: Xbp1 RNA-seq (sL2fc)")) +
  annotate(geom = "text", 
           x = c(1000, 18500), y = c(0.1, 0.1), 
           label = c("Xbp1 repression",
                     "Xbp1 activation"),
           color = c("red", "blue"), size = 3) +
  annotate(geom = "text", 
           x = 3000, y = -0.6, 
           label = str_c("NES: ", format(gsea.xbp1$NES[1], digits = 3),
                         "\npadj: ", format(gsea.xbp1$p.adjust[1], digits = 3)),
           color = "black", size = 4)  +
  annotate(geom = "text", 
           x = length(geneList.xbp1) - gsea.xbp1$rank[1], y = -0.1, 
           label = format(length(geneList.xbp1) - gsea.xbp1$rank[1], big.mark = ","),
           color = "black", size = 3)

gseaplot.xbp1_vs_bhlha15ATAC +
  theme_classic(base_size = 10) +
  ggsave(paste0(prjDir, "results/gsea/xbp1_vs_Bhlha15_ATAC_ChIP.gsea.Rplot.png"), width = 10, height = 5, dpi = "retina")

gseaplot.xbp1_vs_bhlha15ATAC +
  theme_classic(base_size = 10) +
  ggsave(paste0(prjDir, "results/gsea/xbp1_vs_Bhlha15_ATAC_ChIP.gsea.Rplot.pdf"), width = 10, height = 5, dpi = "retina")

# this is just playing around
gseaplot2(gsea.xbp1, geneSetID = 1:2, 
          title = str_c("Xbp1 RNA-seq ranks",
                        "\n",
                        as_tibble(gsea.xbp1) %>% transmute(desc = str_c(Description, ": NES ", NES)) %>% pull() %>% paste0(collapse = "\n")), 
          pvalue_table = TRUE, 
          base_size = 10,
          color = c(color.bhlha15, "green"),
          subplots = 1,
          rel_heights = c(8, 1)) +
  theme_classic(base_size = 10)  +
  guides(color = FALSE)  +
  geom_vline(xintercept = length(geneList.xbp1) - t.gsea.xbp1$rank,
             color = c("green", color.bhlha15),
             lty = 2) +
  ggsave(paste0(prjDir, "results/gsea/playground/xbp1.gsea.Rplot.v4.png"), width = 10, height = 5, dpi = "retina")

gseaplot2(gsea.xbp1, geneSetID = 1:2, 
          title = str_c("Xbp1 RNA-seq ranks",
                        "\n",
                        as_tibble(gsea.xbp1) %>% transmute(desc = str_c(Description, ": NES ", NES)) %>% pull() %>% paste0(collapse = "\n")), 
          pvalue_table = TRUE, 
          base_size = 10,
          color = c(color.bhlha15, "green"),
          subplots = 1:2,
          rel_heights = c(8, 1)) +
  ggsave(paste0(prjDir, "results/gsea/playground/xbp1.gsea.Rplot.v3.png"), width = 10, height = 5, dpi = "retina")
  
t.gsea.xbp1 %>%
  write_excel_csv(str_c(prjDir, "results/gsea/xbp1.gsea.csv"))
# ------------------------------------------------------------------------------
# log session info
# ------------------------------------------------------------------------------
sink(paste0(prjDir,"data/R/",opf,".sessionInfo.DESeq2.txt"))
sessionInfo()
sink()
save.image(file=paste0(prjDir,"data/R/",opf,".RData"))
