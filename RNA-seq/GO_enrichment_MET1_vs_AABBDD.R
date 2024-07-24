#Delfi Dorussen
#Aim to perform GO enrichment analysis for DEGs in MET1 mutants
#Used https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2020/RNAseq/extended_html/06_Gene_set_testing.html
#for GO enrichment plots

setwd("U:/Year1/met1 Mutants Project/RNA_seq_analysis/Output/")
mydir <- getwd()

#BiocManager::install("goseq")
library(goseq)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

#supportedOrganisms() #T. aestivum not in pre-loaded database so will have to load in manually

#Load in all genes with GO terms 
all_go <- read.csv("IWGSC_stress_GO.csv") #file from https://github.com/Borrill-Lab/WheatFlagLeafSenescence/blob/master/data/IWGSC_stress_GO.csv
head(all_go) #Gene IDs are v1.0, so need to convert to v1.1

v1.0_v1.1 <- read.csv("genes_to_transfer_qcov90_pident99_same_ID.csv") #file from https://github.com/Borrill-Lab/WheatFlagLeafSenescence/blob/master/data/genes_to_transfer_qcov90_pident99_same_ID.csv
all_go <- merge(all_go, v1.0_v1.1, by.x="Gene", by.y="gene_v1.0")
all_go <- dplyr::select(all_go, ID:gene_v1.1)

#Load in gene lengths (produced in import_kallisto_quantifications.R)
gene_lengths <- read.csv("MET1_gene_lengths.csv")
gene_lengths <- gene_lengths[!grepl("LC",gene_lengths$X),]
colnames(gene_lengths) <- c("Gene", "Length")

#Ensure all genes with GO terms have gene lengths and vice versa
all_go <- subset(all_go, gene_v1.1 %in% gene_lengths$Gene)
colnames(all_go) <- c("GO","Gene")
all_go <- all_go[,c(2,1)]
gene_lengths <- subset(gene_lengths, Gene %in% all_go$Gene)
gene_lengths_vector <- as.vector(gene_lengths$Length)

#write.csv(all_go, "GO_gene_IDs.csv")

#Looking at DEGs in the 5/6 mutant compared to AABBDD
DEGs_5_6 <- read.csv("AABBDD_vs_5out6_sig.csv")
DEGs_5_6_up <- subset(DEGs_5_6, log2FoldChange>=1)[,1]
DEGs_5_6_down <- subset(DEGs_5_6, log2FoldChange<=1)[,1]

#Up-regulated genes
#Create list of all genes, specifying with are DEGs and which are not
gene_vector_up_5_6 <- as.integer(gene_lengths$Gene %in% DEGs_5_6_up)
table(gene_vector_up_5_6)

gene_vector_down_5_6 <- as.integer(gene_lengths$Gene %in% DEGs_5_6_down)
table(gene_vector_down_5_6)

pwf_up_5_6 <- nullp(DEgenes=gene_vector_up_5_6, bias.data=gene_lengths_vector, plot.fit=TRUE)
rownames(pwf_up_5_6) <- gene_lengths$Gene
GO.up.5_6 <- goseq(pwf_up_5_6, gene2cat=all_go) 
GO.up.BP.5_6 <- subset(GO.up.5_6, ontology=="BP")
GO.up.BP.5_6 <- subset(GO.up.BP.5_6, category!="GO:0034402") #Remove obsolete GO term
GO.up.BP.5_6$p.adjust <- p.adjust(GO.up.BP.5_6$over_represented_pvalue, method="BH")

##Filter to keep only GO terms with more than 15 DEGs in cat & p.adj<0.05
GO.up.BP.5_6.15DEGs <- subset(GO.up.BP.5_6, numDEInCat>15&p.adjust<0.05)

GO.up.BP.5_6.15DEGs <- GO.up.BP.5_6.15DEGs %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>%
  top_n(15, wt=-p.adjust) %>%
  arrange(hitsPerc)

GO.up.BP.5_6.15DEGs$term <- factor(GO.up.BP.5_6.15DEGs$term,levels=GO.up.BP.5_6.15DEGs$term)

ggplot(GO.up.BP.5_6.15DEGs, aes(x=hitsPerc, y=term, size=numDEInCat)) +
  geom_point() +
  expand_limits(x=10) +
  labs(x="Hits (%)", y="GO term", size="DEGs in Category") +
  theme_bw() +
  xlim(0,30) +
  geom_segment( aes(x=0, xend=hitsPerc, y=term, yend=term), linewidth=1)
ggsave("5_6_Up_GO.pdf", width=7.5, height=6, units="in")

#Down-regulated genes
pwf_down_5_6 <- nullp(DEgenes=gene_vector_down_5_6, bias.data=gene_lengths_vector, plot.fit=TRUE)
rownames(pwf_down_5_6) <- gene_lengths$Gene
GO.down.5_6 <- goseq(pwf_down_5_6, gene2cat=all_go) 
GO.down.BP.5_6 <- subset(GO.down.5_6, ontology=="BP")
GO.down.BP.5_6$p.adjust <- p.adjust(GO.down.BP.5_6$over_represented_pvalue, method="BH")

##Filter to keep only GO terms with more than 15 DEGs in cat & p.adj<0.05
GO.down.BP.5_6.15DEGs <- subset(GO.down.BP.5_6, numDEInCat>15&p.adjust<0.05)

GO.down.BP.5_6.15DEGs <- GO.down.BP.5_6.15DEGs %>% 
  mutate(hitsPerc=numDEInCat*100/numInCat) %>%
  top_n(15, wt=-p.adjust) %>%
  arrange(hitsPerc)

GO.down.BP.5_6.15DEGs$term <- factor(GO.down.BP.5_6.15DEGs$term,levels=GO.down.BP.5_6.15DEGs$term)

ggplot(GO.down.BP.5_6.15DEGs, aes(x=hitsPerc, y=term, size=numDEInCat)) +
  geom_point() +
  expand_limits(x=10) +
  labs(x="Hits (%)", y="GO term", size="DEGs in Category") +
  theme_bw() +
  xlim(0,70) +
  geom_segment( aes(x=0, xend=hitsPerc, y=term, yend=term), size=1)
ggsave("5_6_down_GO.pdf", width=9.5, height=6, units="in")
