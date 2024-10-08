#Delfi Dorussen
#Find out how many TEs intersect with DMRs in the met1 mutants

library(dplyr)
library(tidyr)
library(ggplot2)

setwd("Z:/Delfi/MET1_WGBS/DMR_overlap")

##Transposons
mutants <- c("AB_double","AD_double","BD_double","Aabbdd")
mutants2 <- c("AB_double","AD_double","BD_double","5out6")
  
number_DET_DMRs <- NULL

    
for(i in 1:4){
    mutant <- mutants[i]
    mutant2 <- mutants2[i]
    
    setwd("Z:/Delfi/MET1_WGBS/DMR_overlap")
    
    mutant_overlap_hypo <- read.table(paste(mutant,"_overlap_hypoDMRs_TEs.txt",sep=""), sep="\t")
    
    mutant_overlap_hypo <- separate_wider_delim(mutant_overlap_hypo, cols=V9, delim=";",
                                                names=c("ID","Name","Ontology_term","Compo","Status"),
                                                too_few="align_start", too_many="drop")
    mutant_overlap_hypo <- separate_wider_delim(mutant_overlap_hypo, cols=ID, delim="=",
                                                names=c("x","ID"))
    
    mutant_overlap_hypo$ID <- factor(mutant_overlap_hypo$ID)
    mutant_overlap_hypo_unique <- mutant_overlap_hypo %>%
      group_by(ID) %>%
      filter(row_number()==1)
    
    setwd("U:/Year1/met1 Mutants Project/TE_RNA_seq_analysis/")
    mutant_DETs <- read.csv(paste("AABBDD_vs_",mutant2,"_sig.csv", sep=""))
    mutant_DETs <- subset(mutant_DETs, log2FoldChange>=1|log2FoldChange<=-1)
    
    mutant_DETs_hypoDMRs <- merge(mutant_overlap_hypo_unique, mutant_DETs, by.x="ID", by.y="X")
    mutant_DETs_up_hypo <- subset(mutant_DETs_hypoDMRs, log2FoldChange>=1) 
    mutant_DETs_down_hypo <- subset(mutant_DETs_hypoDMRs, log2FoldChange<=-1) 
    
    setwd("Z:/Delfi/MET1_WGBS/DMR_overlap")
    
    mutant_overlap_hyper <- read.table(paste(mutant,"_overlap_hyperDMRs.txt",sep=""), sep="\t")
    
    mutant_overlap_hyper <- separate_wider_delim(mutant_overlap_hyper, cols=V9, delim=";",
                                                 names=c("ID","Name","Ontology_term","Compo","Status"),
                                                 too_few="align_start", too_many="drop")
    mutant_overlap_hyper <- separate_wider_delim(mutant_overlap_hyper, cols=ID, delim="=",
                                                 names=c("x","ID"))
    
    mutant_overlap_hyper$ID <- factor(mutant_overlap_hyper$ID)
    mutant_overlap_hyper_unique <- mutant_overlap_hyper %>%
      group_by(ID) %>%
      filter(row_number()==1)
    
    mutant_DETs_hyperDMRs <- merge(mutant_overlap_hyper_unique, mutant_DETs, by.x="ID", by.y="X")
    mutant_DETs_up_hyper <- subset(mutant_DETs_hyperDMRs, log2FoldChange>=1) 
    mutant_DETs_down_hyper <- subset(mutant_DETs_hyperDMRs, log2FoldChange<=-1)
    
    setwd("U:/Year1/met1 Mutants Project/TE_RNA_seq_analysis/")
    write.csv(mutant_DETs_hypoDMRs, paste(mutant,"_DETs_hypoDMRs.csv",sep=""))
    write.csv(mutant_DETs_hyperDMRs, paste(mutant,"_DETs_hyperDMRs.csv",sep=""))
    
    #Get total number of DETs in the mutant (remove ChrUn)
    mutant_number_DET_DMRs <- cbind(mutant, nrow(mutant_DETs), nrow(mutant_DETs_hypoDMRs), nrow(mutant_DETs_up_hypo), nrow(mutant_DETs_down_hypo),
                                    nrow(mutant_DETs_hyperDMRs), nrow(mutant_DETs_up_hyper), nrow(mutant_DETs_down_hyper))
    number_DET_DMRs <- rbind(number_DET_DMRs, mutant_number_DET_DMRs)
}

colnames(number_DET_DMRs) <- c("Mutant","Number.DETs","DET.hypoDMRs","Up.DET.hypoDMRs","Down.DET.hypoDMRs",
                               "DET.hyperDMRs","Up.DET.hyperDMRs","Down.DET.hyperDMRs")

number_DET_DMRs <- data.frame(number_DET_DMRs)
number_DET_DMRs$Number.DETs <- as.numeric(number_DET_DMRs$Number.DETs)
number_DET_DMRs$DET.hypoDMRs <- as.numeric(number_DET_DMRs$DET.hypoDMRs)
number_DET_DMRs$DET.hyperDMRs <- as.numeric(number_DET_DMRs$DET.hyperDMRs)

number_DET_DMRs$Total_HypoDET_percent <- (number_DET_DMRs$DET.hypoDMRs/number_DET_DMRs$Number.DETs)*100
number_DET_DMRs$Total_HyperDET_percent <- (number_DET_DMRs$DET.hyperDMRs/number_DET_DMRs$Number.DETs)*100

number_DET_DMRs <- pivot_longer(number_DET_DMRs, cols=Total_HypoDET_percent:Total_HyperDET_percent,
                                names_to="Type",values_to="Percentage")
number_DET_DMRs$Mutant <- factor(number_DET_DMRs$Mutant, levels=c("AB_double","AD_double","BD_double",
                                                                  "Aabbdd"))

ggplot(number_DET_DMRs, aes(x=Mutant, y=Percentage, fill=Type)) +
  geom_bar(position="stack",stat="identity", colour="black") +
  theme_bw() +
  theme(legend.position="none")

setwd("U:/Year1/met1 Mutants Project/TE_RNA_seq_analysis/")
ggsave("percentage_DET_DMRs.pdf", width=3, height=5, units="in")


#Plot for just hypo-DMRs
number_DET_hypoDMRs <- select(number_DET_DMRs, Mutant, DET.hypoDMRs:Down.DET.hypoDMRs)
number_DET_hypoDMRs <- number_DET_hypoDMRs[c(1,3,5,7),]
number_DET_hypoDMRs <- pivot_longer(number_DET_hypoDMRs, cols=Up.DET.hypoDMRs:Down.DET.hypoDMRs,
                                    names_to="Type",values_to="Number")
number_DET_hypoDMRs$Number <- as.numeric(number_DET_hypoDMRs$Number)
number_DET_hypoDMRs$Percentage <- (number_DET_hypoDMRs$Number/number_DET_hypoDMRs$DET.hypoDMRs)*100

ggplot(number_DET_hypoDMRs, aes(x=Mutant, y=Percentage, fill=Type)) +
  geom_bar(position="stack",stat="identity", colour="black") +
  theme_bw() +
  theme(legend.position="none")

setwd("U:/Year1/met1 Mutants Project/TE_RNA_seq_analysis/")
ggsave("percentage_DET_hypoDMRs.pdf", width=3, height=5, units="in")


##Genes
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


#Combined plot for hypoDMRs for genes and TEs
number_DET_hypoDMRs <- separate_wider_delim(number_DET_hypoDMRs, cols=Type, delim=".",
                                            names=c("Type","Feature",NA))
colnames(number_DET_hypoDMRs) <- c("Mutant","DE.hypoDMRs","Type","Feature","Number","Percentage")

number_DEG_hypoDMRs <- separate_wider_delim(number_DEG_hypoDMRs, cols=Type, delim=".",
                                            names=c("Type","Feature",NA))
colnames(number_DEG_hypoDMRs) <- c("Mutant","DE.hypoDMRs","Type","Feature","Number","Percentage")

number_DE_hypoDMRs <- rbind(number_DEG_hypoDMRs, number_DET_hypoDMRs)


setwd("U:/Year1/met1 Mutants Project/TE_RNA_seq_analysis/")
ggplot(number_DE_hypoDMRs, aes(x=Feature, y=Percentage, fill=Type)) +
  geom_bar(position="stack", stat="identity", colour="black") +
  facet_wrap(vars(Mutant), nrow=1) +
  theme_bw()
ggsave("hypoDMRs_TEs_DEGs_all.pdf", width=9, height=4, units="in", dpi=300)

#only keeping AB_double and Aabbdd
ggplot(subset(number_DE_hypoDMRs, Mutant=="AB_double"|Mutant=="Aabbdd"), aes(x=Feature, y=Percentage, fill=Type)) +
  geom_bar(position="stack", stat="identity", colour="black") +
  facet_wrap(vars(Mutant), nrow=1) +
  theme_bw()
ggsave("hypoDMRs_TEs_DEGs_Aabbdd.pdf", width=5, height=6, units="in", dpi=300)
