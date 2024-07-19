#Delfi Dorussen
#Differential expression analysis for MET1 mutant lines
#This should be run after import_kallisto_quantifications.R

#BiocManager::install("DESeq2")
library("DESeq2")
library("dplyr")
library("tidyr")
library(ggplot2)

setwd("U:/Year1/met1 Mutants Project/RNA_seq_analysis/Output/")

#Import tpm data
abundances <- read.table("MET1_tpm.tsv")
head(abundances)

#Remove low confidence (LC) genes 
abundances <- cbind(rownames(abundances), data.frame(abundances, row.names=NULL))
colnames(abundances)[1] <- "gene"
head(abundances)
abundances <- abundances[!grepl("LC",abundances$gene),]

#Import count data and subset in the same way as the tpm data
counts <- read.table("MET1_count.tsv")
head(counts)
counts <- counts[rownames(counts) %in% abundances$gene,]

#Provide metadata to the count matrices
setwd("U:/Year1/met1 Mutants Project/RNA_seq_analysis/")
samples <- read.csv("sample_index.csv")
setwd("U:/Year1/met1 Mutants Project/RNA_seq_analysis/Output")
rownames(samples) <- samples$Name
samples <- samples[,c("Genotype","Type")]
samples$Genotype <- factor(samples$Genotype)
samples$Type <- factor(samples$Type)

#Check if metadata is in the correct order - this command should give output 'TRUE'
all(rownames(samples) == colnames(counts))

#DESeq2 Analysis
dds1 <- DESeqDataSetFromMatrix(countData=round(counts), colData=samples, design=~Genotype)

#Producing RangedSummarizedExperiment object & plotting PCA
RSE <- rlog(dds1, blind=FALSE)
plotPCA(RSE, intgroup="Genotype", ntop=1000, pcsToUse=1:2) +
  theme_bw()
ggsave("PCA_by_genotype.png", width=6, height=6)
plotPCA(RSE, intgroup="Type", ntop=1000, pcsToUse=1:2) +
  theme_bw()
ggsave("PCA_by_category.png", width=6, height=6)

#Comparing 5/6 mutant to WT segregant
dds1_DE <- DESeq(dds1)

res_AABBDD_vs_5out6 <- results(dds1_DE, contrast=c("Genotype","Aabbdd","AABBDD"))
AABBDD_vs_5out6_sig <- subset(res_AABBDD_vs_5out6, padj<0.01) #6444 DEGs
AABBDD_vs_5out6_sig_up <- subset(AABBDD_vs_5out6_sig, log2FoldChange>=1) #2004 upregulated DEGs
AABBDD_vs_5out6_sig_down <- subset(AABBDD_vs_5out6_sig, log2FoldChange<=-1) #3900 downregulated DEGs
write.csv(AABBDD_vs_5out6_sig, "AABBDD_vs_5out6_sig.csv")

#Single Mutants
res_AABBDD_vs_aaBBDD <- results(dds1_DE, contrast=c("Genotype","aaBBDD","AABBDD"))
AABBDD_vs_aaBBDD_sig <- subset(res_AABBDD_vs_aaBBDD, padj<0.01) #4 DEGs
AABBDD_vs_aaBBDD_sig_up <- subset(AABBDD_vs_aaBBDD_sig, log2FoldChange>=1) #1 upregulated DEGs
AABBDD_vs_aaBBDD_sig_down <- subset(AABBDD_vs_aaBBDD_sig, log2FoldChange<=-1) #3 downregulated DEGs
write.csv(AABBDD_vs_aaBBDD_sig, "AABBDD_vs_A_single_sig.csv")

res_AABBDD_vs_AAbbDD <- results(dds1_DE, contrast=c("Genotype","AAbbDD","AABBDD"))
AABBDD_vs_AAbbDD_sig <- subset(res_AABBDD_vs_AAbbDD, padj<0.01) #20 DEGs
AABBDD_vs_AAbbDD_sig_up <- subset(AABBDD_vs_AAbbDD_sig, log2FoldChange>=1) #11 upregulated DEGs
AABBDD_vs_AAbbDD_sig_down <- subset(AABBDD_vs_AAbbDD_sig, log2FoldChange<=-1) #8 downregulated DEGs
write.csv(AABBDD_vs_AAbbDD_sig, "AABBDD_vs_B_single_sig.csv")

res_AABBDD_vs_AABBdd <- results(dds1_DE, contrast=c("Genotype","AABBdd","AABBDD"))
AABBDD_vs_AABBdd_sig <- subset(res_AABBDD_vs_AABBdd, padj<0.01) #4 DEGs
AABBDD_vs_AABBdd_sig_up <- subset(AABBDD_vs_AABBdd_sig, log2FoldChange>=1) #3 upregulated DEGs
AABBDD_vs_AABBdd_sig_down <- subset(AABBDD_vs_AABBdd_sig, log2FoldChange<=-1) #0 downregulated DEGs
write.csv(AABBDD_vs_AABBdd_sig, "AABBDD_vs_D_single_sig.csv")

#Double Mutants
res_AABBDD_vs_aabbDD <- results(dds1_DE, contrast=c("Genotype","aabbDD","AABBDD"))
AABBDD_vs_aabbDD_sig <- subset(res_AABBDD_vs_aabbDD, padj<0.01) #146 DEGs
AABBDD_vs_aabbDD_sig_up <- subset(AABBDD_vs_aabbDD_sig, log2FoldChange>=1) #118 upregulated DEGs
AABBDD_vs_aabbDD_sig_down <- subset(AABBDD_vs_aabbDD_sig, log2FoldChange<=-1) #18 downregulated DEGs
write.csv(AABBDD_vs_aabbDD_sig, "AABBDD_vs_AB_double_sig.csv")

res_AABBDD_vs_aaBBdd <- results(dds1_DE, contrast=c("Genotype","aaBBdd","AABBDD"))
AABBDD_vs_aaBBdd_sig <- subset(res_AABBDD_vs_aaBBdd, padj<0.01) #46 DEGs
AABBDD_vs_aaBBdd_sig_up <- subset(AABBDD_vs_aaBBdd_sig, log2FoldChange>=1) #39 upregulated DEGs
AABBDD_vs_aaBBdd_sig_down <- subset(AABBDD_vs_aaBBdd_sig, log2FoldChange<=-1) #5 downregulated DEGs
write.csv(AABBDD_vs_aaBBdd_sig, "AABBDD_vs_AD_double_sig.csv")

res_AABBDD_vs_AAbbdd <- results(dds1_DE, contrast=c("Genotype","AAbbdd","AABBDD"))
AABBDD_vs_AAbbdd_sig <- subset(res_AABBDD_vs_AAbbdd, padj<0.01) #59 DEGs
AABBDD_vs_AAbbdd_sig_up <- subset(AABBDD_vs_AAbbdd_sig, log2FoldChange>=1) #46 upregulated DEGs
AABBDD_vs_AAbbdd_sig_down <- subset(AABBDD_vs_AAbbdd_sig, log2FoldChange<=-1) #11 downregulated DEGs
write.csv(AABBDD_vs_AAbbdd_sig, "AABBDD_vs_BD_double_sig.csv")

double_overlap_1 <- intersect(rownames(AABBDD_vs_aabbDD_sig),rownames(AABBDD_vs_aaBBdd_sig)) #41 DEGs
double_overlap_2 <- intersect(double_overlap_1, rownames(AABBDD_vs_AAbbdd_sig)) #26 DEGs
write.csv(double_overlap_2, "double_mutants_vs_AABBDD_DEGs.csv")


#Overlap between 5/6 and at least one double mutant
aabbDD_5out6_overlap_up <- intersect(rownames(AABBDD_vs_aabbDD_sig_up), rownames(AABBDD_vs_5out6_sig_up)) #73 DEGS
aaBBdd_5out6_overlap_up <- intersect(rownames(AABBDD_vs_aaBBdd_sig_up), rownames(AABBDD_vs_5out6_sig_up)) #26 DEG
AAbbdd_5out6_overlap_up <- intersect(rownames(AABBDD_vs_AAbbdd_sig_up), rownames(AABBDD_vs_5out6_sig_up)) #30 DEGS
double_up_1 <- append(aabbDD_5out6_overlap_up, aaBBdd_5out6_overlap_up)
double_up_2 <- as.data.frame(append(double_up_1, AAbbdd_5out6_overlap_up))

double_5out6_overlap_up <- distinct(double_up_2) #76 DEGs
write.csv(double_5out6_overlap_up, "double_5out6_overlap_up.csv")

aabbDD_5out6_overlap_down <- intersect(rownames(AABBDD_vs_aabbDD_sig_down), rownames(AABBDD_vs_5out6_sig_down)) #9 DEGS
aaBBdd_5out6_overlap_down <- intersect(rownames(AABBDD_vs_aaBBdd_sig_down), rownames(AABBDD_vs_5out6_sig_down)) #1 DEG
AAbbdd_5out6_overlap_down <- intersect(rownames(AABBDD_vs_AAbbdd_sig_down), rownames(AABBDD_vs_5out6_sig_down)) #6 DEGS
double_down_1 <- append(aabbDD_5out6_overlap_down, aaBBdd_5out6_overlap_down)
double_down_2 <- as.data.frame(append(double_down_1, AAbbdd_5out6_overlap_down))

double_5out6_overlap_down <- distinct(double_down_2) #10 DEGs
write.csv(double_5out6_overlap_down, "double_5out6_overlap_down.csv")
