#Delfi Dorussen
#Find out how many genes intersect with DMRs in the met1 mutants

library(dplyr)
library(tidyr)
library(ggplot2)

setwd("Z:/Delfi/MET1_WGBS/DMR_overlap")

mutants <- c("A_single","B_single","D_single","AB_double","AD_double","BD_double","Aabbdd")

number_gene_DMRs <- NULL

for(i in 1:7){
  mutant <- mutants[i]
  
  mutant_overlap_hypo <- read.table(paste(mutant,"_overlap_hypoDMRs.txt",sep=""), sep="\t")
  
  mutant_overlap_hypo <- separate_wider_delim(mutant_overlap_hypo, cols=V9, delim=";",
                                              names=c("ID","previous_id","primconf","Name"),
                                              too_few="align_start", too_many="drop")
  mutant_overlap_hypo <- separate_wider_delim(mutant_overlap_hypo, cols=ID, delim="=",
                                              names=c("x","TraesID"))
  
  mutant_overlap_hypo$ID <- factor(mutant_overlap_hypo$TraesID)
  mutant_overlap_hypo_unique <- mutant_overlap_hypo %>%
    group_by(TraesID) %>%
    filter(row_number()==1)
  
  mutant_overlap_hyper <- read.table(paste(mutant,"_overlap_hyperDMRs.txt",sep=""), sep="\t")
  
  mutant_overlap_hyper <- separate_wider_delim(mutant_overlap_hyper, cols=V9, delim=";",
                                               names=c("ID","previous_id","primconf","Name"),
                                               too_few="align_start", too_many="drop")
  mutant_overlap_hyper <- separate_wider_delim(mutant_overlap_hyper, cols=ID, delim="=",
                                               names=c("x","TraesID"))
  
  mutant_overlap_hyper$ID <- factor(mutant_overlap_hyper$TraesID)
  mutant_overlap_hyper_unique <- mutant_overlap_hyper %>%
    group_by(TraesID) %>%
    filter(row_number()==1)
  
  mutant_overlap_unique <- rbind(mutant_overlap_hyper_unique, mutant_overlap_hypo_unique) %>%
    group_by(TraesID) %>%
    filter(row_number()==1)
  
  mutant_number_gene_DMRs <- cbind(mutant, nrow(mutant_overlap_unique), nrow(mutant_overlap_hypo_unique), nrow(mutant_overlap_hyper_unique))
  number_gene_DMRs <- rbind(number_gene_DMRs, mutant_number_gene_DMRs)
}

colnames(number_gene_DMRs) <- c("Mutant","Total gene-DMRs","gene_Hypo-DMRs","gene_Hyper-DMRs")

#Total number of genes (excluding chrUn) is 105200
number_gene_DMRs <- data.frame(number_gene_DMRs)
number_gene_DMRs$Total.gene.DMRs <- as.numeric(number_gene_DMRs$Total.gene.DMRs)
number_gene_DMRs$gene_Hypo.DMRs <- as.numeric(number_gene_DMRs$gene_Hypo.DMRs)
number_gene_DMRs$gene_Hyper.DMRs <- as.numeric(number_gene_DMRs$gene_Hyper.DMRs)

number_gene_DMRs$Total_gene_percent <- (number_gene_DMRs$Total.gene.DMRs/105200)*100
number_gene_DMRs$Total_Hypogene_percent <- (number_gene_DMRs$gene_Hypo.DMRs/105200)*100
number_gene_DMRs$Total_Hypegene_percent <- (number_gene_DMRs$gene_Hyper.DMRs/105200)*100

number_gene_DMRs <- pivot_longer(number_gene_DMRs, cols=Total_Hypogene_percent:Total_Hypegene_percent,
                               names_to="Type",values_to="Percentage")
number_gene_DMRs$Mutant <- factor(number_gene_DMRs$Mutant, levels=c("A_single","B_single","D_single",
                                                                "AB_double","AD_double","BD_double",
                                                                "Aabbdd"))

ggplot(number_gene_DMRs, aes(x=Mutant, y=Percentage, fill=Type)) +
  geom_bar(position="stack",stat="identity", colour="black") +
  theme_bw() +
  theme(legend.position="none")

setwd("U:/Year1/met1 Mutants Project/WGBS_analysis/")
ggsave("percentage_gene_DMRs.pdf", width=5, height=5, units="in")


##Incorporating DEGs into the analysis
mutants <- c("AB_double","AD_double","BD_double","Aabbdd")
mutants2 <- c("AB_double","AD_double","BD_double","5out6")

number_DEG_DMRs <- NULL

for(i in 1:4){
  mutant <- mutants[i]
  mutant2 <- mutants2[i]
  
  setwd("Z:/Delfi/MET1_WGBS/DMR_overlap")
  
  mutant_overlap_hypo <- read.table(paste(mutant,"_overlap_hypoDMRs.txt",sep=""), sep="\t")
  
  mutant_overlap_hypo <- separate_wider_delim(mutant_overlap_hypo, cols=V9, delim=";",
                                              names=c("ID","previous_id","primconf","Name"),
                                              too_few="align_start", too_many="drop")
  mutant_overlap_hypo <- separate_wider_delim(mutant_overlap_hypo, cols=ID, delim="=",
                                              names=c("x","TraesID"))
  
  mutant_overlap_hypo$ID <- factor(mutant_overlap_hypo$TraesID)
  mutant_overlap_hypo_unique <- mutant_overlap_hypo %>%
    group_by(TraesID) %>%
    filter(row_number()==1)
  
  setwd("U:/Year1/met1 Mutants Project/RNA_seq_analysis/Output/")
  mutant_DEGs <- read.csv(paste("AABBDD_vs_",mutant2,"_sig.csv", sep=""))
  mutant_DEGs <- subset(mutant_DEGs, log2FoldChange>=1|log2FoldChange<=-1)
  
  mutant_DEGs_hypoDMRs <- merge(mutant_overlap_hypo_unique, mutant_DEGs, by.x="TraesID", by.y="X")
  mutant_DEGs_up_hypo <- subset(mutant_DEGs_hypoDMRs, log2FoldChange>=1) 
  mutant_DEGs_down_hypo <- subset(mutant_DEGs_hypoDMRs, log2FoldChange<=-1) 
  
  setwd("Z:/Delfi/MET1_WGBS/DMR_overlap")
  
  mutant_overlap_hyper <- read.table(paste(mutant,"_overlap_hyperDMRs.txt",sep=""), sep="\t")
  
  mutant_overlap_hyper <- separate_wider_delim(mutant_overlap_hyper, cols=V9, delim=";",
                                               names=c("ID","previous_id","primconf","Name"),
                                               too_few="align_start", too_many="drop")
  mutant_overlap_hyper <- separate_wider_delim(mutant_overlap_hyper, cols=ID, delim="=",
                                               names=c("x","TraesID"))
  
  mutant_overlap_hyper$ID <- factor(mutant_overlap_hyper$TraesID)
  mutant_overlap_hyper_unique <- mutant_overlap_hyper %>%
    group_by(TraesID) %>%
    filter(row_number()==1)
  
  mutant_DEGs_hyperDMRs <- merge(mutant_overlap_hyper_unique, mutant_DEGs, by.x="TraesID", by.y="X")
  mutant_DEGs_up_hyper <- subset(mutant_DEGs_hyperDMRs, log2FoldChange>=1) 
  mutant_DEGs_down_hyper <- subset(mutant_DEGs_hyperDMRs, log2FoldChange<=-1)
  
  setwd("U:/Year1/met1 Mutants Project/WGBS_analysis/Output/")
  write.csv(mutant_DEGs_hypoDMRs, paste(mutant,"_DEGs_hypoDMRs.csv",sep=""))
  write.csv(mutant_DEGs_hyperDMRs, paste(mutant,"_DEGs_hyperDMRs.csv",sep=""))
  
  #Get total number of DEGs in the mutant (remove ChrUn)
  mutant_DEGs <- subset(mutant_DEGs, !grepl("CSU",X))
  
  mutant_number_DEG_DMRs <- cbind(mutant, nrow(mutant_DEGs), nrow(mutant_DEGs_hypoDMRs), nrow(mutant_DEGs_up_hypo), nrow(mutant_DEGs_down_hypo),
                                  nrow(mutant_DEGs_hyperDMRs), nrow(mutant_DEGs_up_hyper), nrow(mutant_DEGs_down_hyper))
  number_DEG_DMRs <- rbind(number_DEG_DMRs, mutant_number_DEG_DMRs)
}

colnames(number_DEG_DMRs) <- c("Mutant","Number.DEGs","DEG.hypoDMRs","Up.DEG.hypoDMRs","Down.DEG.hypoDMRs",
                               "DEG.hyperDMRs","Up.DEG.hyperDMRs","Down.DEG.hyperDMRs")

number_DEG_DMRs <- data.frame(number_DEG_DMRs)
number_DEG_DMRs$Number.DEGs <- as.numeric(number_DEG_DMRs$Number.DEGs)
number_DEG_DMRs$DEG.hypoDMRs <- as.numeric(number_DEG_DMRs$DEG.hypoDMRs)
number_DEG_DMRs$DEG.hyperDMRs <- as.numeric(number_DEG_DMRs$DEG.hyperDMRs)

number_DEG_DMRs$Total_HypoDEG_percent <- (number_DEG_DMRs$DEG.hypoDMRs/number_DEG_DMRs$Number.DEGs)*100
number_DEG_DMRs$Total_HyperDEG_percent <- (number_DEG_DMRs$DEG.hyperDMRs/number_DEG_DMRs$Number.DEGs)*100

number_DEG_DMRs <- pivot_longer(number_DEG_DMRs, cols=Total_HypoDEG_percent:Total_HyperDEG_percent,
                                 names_to="Type",values_to="Percentage")
number_DEG_DMRs$Mutant <- factor(number_DEG_DMRs$Mutant, levels=c("AB_double","AD_double","BD_double",
                                                                    "Aabbdd"))

ggplot(number_DEG_DMRs, aes(x=Mutant, y=Percentage, fill=Type)) +
  geom_bar(position="stack",stat="identity", colour="black") +
  theme_bw() +
  theme(legend.position="none")

setwd("U:/Year1/met1 Mutants Project/WGBS_analysis/")
ggsave("percentage_DEG_DMRs.pdf", width=3, height=5, units="in")


#Plot for just hypo-DMRs
number_DEG_hypoDMRs <- select(number_DEG_DMRs, Mutant, DEG.hypoDMRs:Down.DEG.hypoDMRs)
number_DEG_hypoDMRs <- number_DEG_hypoDMRs[c(1,3,5,7),]
number_DEG_hypoDMRs <- pivot_longer(number_DEG_hypoDMRs, cols=Up.DEG.hypoDMRs:Down.DEG.hypoDMRs,
                                    names_to="Type",values_to="Number")
number_DEG_hypoDMRs$Number <- as.numeric(number_DEG_hypoDMRs$Number)
number_DEG_hypoDMRs$Percentage <- (number_DEG_hypoDMRs$Number/number_DEG_hypoDMRs$DEG.hypoDMRs)*100

ggplot(number_DEG_hypoDMRs, aes(x=Mutant, y=Percentage, fill=Type)) +
  geom_bar(position="stack",stat="identity", colour="black") +
  theme_bw() +
  theme(legend.position="none")

setwd("U:/Year1/met1 Mutants Project/WGBS_analysis/")
ggsave("percentage_DEG_hypoDMRs.pdf", width=3, height=5, units="in")
