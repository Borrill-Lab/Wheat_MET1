library(ggplot2)
library(dplyr)
library(tidyr)
library(ggforce)

setwd("U:/Year1/met1 Mutants Project/RNA_seq_analysis/Output/")

number_DEGs <- read.csv("number_of_DEGs_AABBDD.csv")
number_DEGs <- pivot_longer(number_DEGs, Down.DEGs.FC1:Up.DEGs.FC1, names_to="Class", values_to="Number")
#number_DEGs <- subset(number_DEGs, Class!="DEGs")
number_DEGs$Genotype <- factor(number_DEGs$Genotype, levels=c("aaBBDD","AAbbDD","AABBdd","aabbDD","aaBBdd","AAbbdd","Aabbdd"))
number_DEGs$Class <- factor(number_DEGs$Class, levels=c("Down.DEGs.FC1","Up.DEGs.FC1"))

ggplot(number_DEGs, aes(x=Genotype, y=Number, fill=Class)) +
  geom_bar(position="stack", stat="identity", colour="black") +
  theme_bw() +
  scale_fill_manual(values=c("#a4c7ff","#ff6643"), labels=c("Down-DEGs","Up-DEGs")) +
  ylab("Number of DEGs") +
  facet_zoom(ylim = c(0, 200)) +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
ggsave("number_of_DEGs_AABBDD_FC1.pdf", height=6, width=9, units="in")


number_DEGs <- read.csv("number_of_DEGs_AABBDD.csv")
number_DEGs <- pivot_longer(number_DEGs, Down.DEGs.FC1:DEGs.FC1, names_to="Class", values_to="Number")
number_DEGs <- subset(number_DEGs, Class!="Up.DEGs.FC1"&Class!="Down.DEGs.FC1")
number_DEGs$log.Number <- log10(number_DEGs$Number)

ggplot(number_DEGs, aes(x=X..of.AABBDD, y=log.Number, colour=Genotype)) +
  geom_point()

#number_DEGs <- subset(number_DEGs, Genotype!="Aabbdd")
ggplot(number_DEGs, aes(x=X..of.AABBDD, y=log.Number, colour=Type)) +
  geom_point(size=5) +
  theme_bw() +
  scale_colour_manual(values=c("#CC2936","#4056f4","#e9df00")) +
  xlab("% CG Methylation Relative to AABBDD") +
  ylab("log(Number of DEGs)")
ggsave("DEGs_vs_methylation_FC1.pdf",height=6, width=7, units="in")
