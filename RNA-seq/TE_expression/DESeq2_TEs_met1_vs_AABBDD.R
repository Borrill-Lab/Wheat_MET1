#Delfi Dorussen
#Aim to carry out differential expression analysis for transposable elements in
#the met1 mutants, compared to the WT segregant (AABBDD)

library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

setwd("U:/Year1/met1 Mutants Project/TE_RNA_seq_analysis/")

TE_counts <- read.csv("TE_counts.csv", row.names="Feature")

#Count number of TEs that have uniquely mapped reads in at least one sample
TE_counts_unique <- subset(TE_counts, AABBDD_rep1!=0 | AABBDD_rep2!=0 | AABBDD_rep3!=0 |
                             aaBBDD_rep1!=0 | aaBBDD_rep2!=0 | aaBBDD_rep3!=0 |
                             AAbbDD_rep1!=0 | AAbbDD_rep2!=0 | AAbbDD_rep3!=0 |
                             AABBdd_rep1!=0 | AABBdd_rep2!=0 | AABBdd_rep3!=0 |
                             aabbDD_rep1!=0 | aabbDD_rep2!=0 | aabbDD_rep3!=0 |
                             aaBBdd_rep1!=0 | aaBBdd_rep2!=0 | aaBBdd_rep3!=0 |
                             AAbbdd_rep1!=0 | AAbbdd_rep2!=0 | AAbbdd_rep3!=0 |
                             Aabbdd_rep1!=0)
#Keeps 79,878 TEs with uniquely mapped reads


samples <- read.csv("sample_index.csv")
rownames(samples) <- samples$Name
samples <- samples[,c("Genotype","Type")]
samples$Genotype <- factor(samples$Genotype)
samples$Type <- factor(samples$Type)

all(rownames(samples) == colnames(TE_counts_unique)) #Should be TRUE

#DESeq2 Analysis
dds1 <- DESeqDataSetFromMatrix(countData=TE_counts_unique, colData=samples, design=~Genotype)

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
AABBDD_vs_5out6_sig <- subset(res_AABBDD_vs_5out6, padj<0.01) #91 DEGs, 143 if use TE_counts_unique
AABBDD_vs_5out6_sig_up <- subset(AABBDD_vs_5out6_sig, log2FoldChange>=1) #71 upregulated DEGs, 103 if use TE_counts_unique
AABBDD_vs_5out6_sig_down <- subset(AABBDD_vs_5out6_sig, log2FoldChange<=-1) #20 downregulated DEGs, 38 if use TE_counts_unique
write.csv(AABBDD_vs_5out6_sig, "AABBDD_vs_5out6_sig.csv")

#Single Mutants
res_AABBDD_vs_aaBBDD <- results(dds1_DE, contrast=c("Genotype","aaBBDD","AABBDD"))
AABBDD_vs_aaBBDD_sig <- subset(res_AABBDD_vs_aaBBDD, padj<0.01) #0 DEGs

res_AABBDD_vs_AAbbDD <- results(dds1_DE, contrast=c("Genotype","AAbbDD","AABBDD"))
AABBDD_vs_AAbbDD_sig <- subset(res_AABBDD_vs_AAbbDD, padj<0.01) #0 DEGs

res_AABBDD_vs_AABBdd <- results(dds1_DE, contrast=c("Genotype","AABBdd","AABBDD"))
AABBDD_vs_AABBdd_sig <- subset(res_AABBDD_vs_AABBdd, padj<0.01) #0 DEGs

#Double Mutants
res_AABBDD_vs_aabbDD <- results(dds1_DE, contrast=c("Genotype","aabbDD","AABBDD"))
AABBDD_vs_aabbDD_sig <- subset(res_AABBDD_vs_aabbDD, padj<0.01) #49 DEGs, 69 if use TE_counts_unique
AABBDD_vs_aabbDD_sig_up <- subset(AABBDD_vs_aabbDD_sig, log2FoldChange>=1) #46 upregulated DEGs, 65 if use TE_counts_unique
AABBDD_vs_aabbDD_sig_down <- subset(AABBDD_vs_aabbDD_sig, log2FoldChange<=-1) #3 downregulated DEGs, 4 if use TE_counts_unique
write.csv(AABBDD_vs_aabbDD_sig, "AABBDD_vs_AB_double_sig.csv")

res_AABBDD_vs_aaBBdd <- results(dds1_DE, contrast=c("Genotype","aaBBdd","AABBDD"))
AABBDD_vs_aaBBdd_sig <- subset(res_AABBDD_vs_aaBBdd, padj<0.01) #14 DEGs, 18 if use TE_counts_unique
AABBDD_vs_aaBBdd_sig_up <- subset(AABBDD_vs_aaBBdd_sig, log2FoldChange>=1) #14 upregulated DEGs, 18 if use TE_counts_unique
AABBDD_vs_aaBBdd_sig_down <- subset(AABBDD_vs_aaBBdd_sig, log2FoldChange<=-1) #0 downregulated DEGs
write.csv(AABBDD_vs_aaBBdd_sig, "AABBDD_vs_AD_double_sig.csv")

res_AABBDD_vs_AAbbdd <- results(dds1_DE, contrast=c("Genotype","AAbbdd","AABBDD"))
AABBDD_vs_AAbbdd_sig <- subset(res_AABBDD_vs_AAbbdd, padj<0.01) #12 DEGs, 22 if use TE_counts_unique
AABBDD_vs_AAbbdd_sig_up <- subset(AABBDD_vs_AAbbdd_sig, log2FoldChange>=1) #12 upregulated DEGs, 22 if use TE_counts_unique
AABBDD_vs_AAbbdd_sig_down <- subset(AABBDD_vs_AAbbdd_sig, log2FoldChange<=-1) #0 downregulated DEGs
write.csv(AABBDD_vs_AAbbdd_sig, "AABBDD_vs_BD_double_sig.csv")


#Make plot with number of DETs 
Aabbdd <- c(nrow(AABBDD_vs_5out6_sig_up), nrow(AABBDD_vs_5out6_sig_down))
aabbDD <- c(nrow(AABBDD_vs_aabbDD_sig_up), nrow(AABBDD_vs_aabbDD_sig_down))
aaBBdd <- c(nrow(AABBDD_vs_aaBBdd_sig_up), nrow(AABBDD_vs_aaBBdd_sig_down))
AAbbdd <- c(nrow(AABBDD_vs_AAbbdd_sig_up), nrow(AABBDD_vs_AAbbdd_sig_down))

number_DETs <- as.data.frame(rbind(Aabbdd, aabbDD, aaBBdd, AAbbdd))
number_DETs$Genotype <- rownames(number_DETs)
colnames(number_DETs) <- c("Up.DETs","Down.DETs","Genotype")
number_DETs <- pivot_longer(number_DETs, Up.DETs:Down.DETs, names_to="Class", values_to="Number")
number_DETs$Genotype <- factor(number_DETs$Genotype, levels=c("aabbDD","aaBBdd","AAbbdd","Aabbdd"))

ggplot(number_DETs, aes(x=Genotype, y=Number, fill=Class)) +
  geom_bar(position="stack", stat="identity", colour="black") +
  theme_bw() +
  scale_fill_manual(values=c("#a4c7ff","#ff6643"), labels=c("Down-DETs","Up-DETs")) +
  ylab("Number of DETs") +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
ggsave("number_of_DETs_AABBDD_FC1.pdf", height=6, width=6, units="in")
