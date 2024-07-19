#Delfi Dorussen
#Aim to make plots comparing (normalised) levels of expression of particular genes
#of interest in the met1 mutants

#BiocManager::install("DESeq2")
library("DESeq2")
library("dplyr")
library("tidyr")
library(ggplot2)
library(rio)

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
dds1$Type <- relevel(dds1$Genotype,"AABBDD") #Sets WT as reference level
dds1 <- DESeq(dds1)

counts_dds1 <- counts(dds1, normalized=TRUE)

samples$Sample <- rownames(samples)


#Plot individual genes 

#Histone-lysine methyltransferase
counts_dds1_HLM <- counts_dds1[row.names(counts_dds1) %in% c("TraesCS2B02G503200"),]
counts_dds1_HLM <- as.data.frame(counts_dds1_HLM)
counts_dds1_HLM$Sample <- rownames(counts_dds1_HLM)
counts_dds1_HLM <- merge(samples, counts_dds1_HLM)

HLM_summary <- counts_dds1_HLM %>%
  group_by(Genotype) %>%
  summarise(mean=mean(counts_dds1_HLM), SE=sd(counts_dds1_HLM)/sqrt(3))

HLM_summary_min <- subset(HLM_summary, Genotype!="CWT")

counts_dds1_HLM_min <- subset(counts_dds1_HLM, Genotype!="CWT")

HLM_summary_min$Genotype <- factor(HLM_summary_min$Genotype, levels=c("AABBDD","aaBBDD","AAbbDD","AABBdd",
                                                                      "aabbDD","aaBBdd","AAbbdd","Aabbdd"))

colours=c("darkgray","#e9df00","#e9df00","#e9df00","#4056f4","#4056f4","#4056f4","#cc2936")

HLM_bar <- ggplot(HLM_summary_min, aes(x=Genotype, y=mean, fill=Genotype)) +
  geom_bar(stat="identity", alpha=0.8, colour="black") +
  geom_jitter(data=counts_dds1_HLM_min, aes(x=Genotype, y=counts_dds1_HLM), width=0.1, size=2) +
  scale_fill_manual(values=colours) +
  #facet_grid(rows=NULL, cols=vars(gene)) +
  geom_errorbar(aes(ymin=mean-(1.96*SE), ymax=mean+(1.96*SE)), width=.2) +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
  ylab("Normalised Counts")
ggsave("HLM_individual_plot.pdf", height=6, width=4, units="in")

#RdDM gene (B homoeologue as example)
counts_dds1_sub4 <- counts_dds1[row.names(counts_dds1) %in% c("TraesCS6B02G315200"),]
counts_dds1_sub4 <- as.data.frame(counts_dds1_sub4)
counts_dds1_sub4$Sample <- rownames(counts_dds1_sub4)
counts_dds1_sub4 <- merge(samples, counts_dds1_sub4)

sub4_summary <- counts_dds1_sub4 %>%
  group_by(Genotype) %>%
  summarise(mean=mean(counts_dds1_sub4), SE=sd(counts_dds1_sub4)/sqrt(3))

sub4_summary_min <- subset(sub4_summary, Genotype!="CWT")

counts_dds1_sub4_min <- subset(counts_dds1_sub4, Genotype!="CWT")

sub4_summary_min$Genotype <- factor(sub4_summary_min$Genotype, levels=c("AABBDD","aaBBDD","AAbbDD","AABBdd",
                                                                        "aabbDD","aaBBdd","AAbbdd","Aabbdd"))

colours=c("darkgray","#e9df00","#e9df00","#e9df00","#4056f4","#4056f4","#4056f4","#cc2936")

sub4_bar <- ggplot(sub4_summary_min, aes(x=Genotype, y=mean, fill=Genotype)) +
  geom_bar(stat="identity", alpha=0.8, colour="black") +
  geom_jitter(data=counts_dds1_sub4_min, aes(x=Genotype, y=counts_dds1_sub4), width=0.1, size=2) +
  scale_fill_manual(values=colours) +
  #facet_grid(rows=NULL, cols=vars(gene)) +
  geom_errorbar(aes(ymin=mean-(1.96*SE), ymax=mean+(1.96*SE)), width=.2) +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
  ylab("Normalised Counts")
ggsave("sub4_individual_plot.pdf", height=6, width=4, units="in")

#DEMETER-like DNA Glycosylase (DME-like)
counts_dds1_DME <- counts_dds1[row.names(counts_dds1) %in% c("TraesCS3B02G023200"),]
counts_dds1_DME <- as.data.frame(counts_dds1_DME)
counts_dds1_DME$Sample <- rownames(counts_dds1_DME)
counts_dds1_DME <- merge(samples, counts_dds1_DME)

DME_summary <- counts_dds1_DME %>%
  group_by(Genotype) %>%
  summarise(mean=mean(counts_dds1_DME), SE=sd(counts_dds1_DME)/sqrt(3))

DME_summary_min <- subset(DME_summary, Genotype!="CWT")

counts_dds1_DME_min <- subset(counts_dds1_DME, Genotype!="CWT")

DME_summary_min$Genotype <- factor(DME_summary_min$Genotype, levels=c("AABBDD","aaBBDD","AAbbDD","AABBdd",
                                                                        "aabbDD","aaBBdd","AAbbdd","Aabbdd"))

colours=c("darkgray","#e9df00","#e9df00","#e9df00","#4056f4","#4056f4","#4056f4","#cc2936")

DME_bar <- ggplot(DME_summary_min, aes(x=Genotype, y=mean, fill=Genotype)) +
  geom_bar(stat="identity", alpha=0.8, colour="black") +
  geom_jitter(data=counts_dds1_DME_min, aes(x=Genotype, y=counts_dds1_DME), width=0.1, size=2) +
  scale_fill_manual(values=colours) +
  #facet_grid(rows=NULL, cols=vars(gene)) +
  geom_errorbar(aes(ymin=mean-(1.96*SE), ymax=mean+(1.96*SE)), width=.2) +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
  ylab("Normalised Counts")
ggsave("DME_individual_plot.pdf", height=6, width=4, units="in")


#HEN1
counts_dds1_HEN1 <- counts_dds1[row.names(counts_dds1) %in% c("TraesCS2D02G114600"),]
counts_dds1_HEN1 <- as.data.frame(counts_dds1_HEN1)
counts_dds1_HEN1$Sample <- rownames(counts_dds1_HEN1)
counts_dds1_HEN1 <- merge(samples, counts_dds1_HEN1)

HEN1_summary <- counts_dds1_HEN1 %>%
  group_by(Genotype) %>%
  summarise(mean=mean(counts_dds1_HEN1), SE=sd(counts_dds1_HEN1)/sqrt(3))

HEN1_summary_min <- subset(HEN1_summary, Genotype!="CWT")

counts_dds1_HEN1_min <- subset(counts_dds1_HEN1, Genotype!="CWT")

HEN1_summary_min$Genotype <- factor(HEN1_summary_min$Genotype, levels=c("AABBDD","aaBBDD","AAbbDD","AABBdd",
                                                                      "aabbDD","aaBBdd","AAbbdd","Aabbdd"))

colours=c("darkgray","#e9df00","#e9df00","#e9df00","#4056f4","#4056f4","#4056f4","#cc2936")

HEN1_bar <- ggplot(HEN1_summary_min, aes(x=Genotype, y=mean, fill=Genotype)) +
  geom_bar(stat="identity", alpha=0.8, colour="black") +
  geom_jitter(data=counts_dds1_HEN1_min, aes(x=Genotype, y=counts_dds1_HEN1), width=0.1, size=2) +
  scale_fill_manual(values=colours) +
  #facet_grid(rows=NULL, cols=vars(gene)) +
  geom_errorbar(aes(ymin=mean-(1.96*SE), ymax=mean+(1.96*SE)), width=.2) +
  theme_bw() +
  theme(legend.position="none", axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
  ylab("Normalised Counts")
ggsave("HEN1_individual_plot.pdf", height=6, width=4, units="in")
