#Delfi Dorussen
#Find out how many TEs intersect with DMRs in the met1 mutants

library(dplyr)
library(tidyr)
library(ggplot2)

setwd("Z:/Delfi/MET1_WGBS/DMR_overlap")

mutants <- c("A_single","B_single","D_single","AB_double","AD_double","BD_double","Aabbdd")

number_TE_DMRs <- NULL
number_family_DMRs_hypo <- NULL
number_family_DMRs_hyper <- NULL
number_subfamily_DMRs_hypo <- NULL
number_subfamily_DMRs_hyper <- NULL

for(i in 1:7){
  mutant <- mutants[i]
  
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
  
  #Separate out by TE family
  mutant_overlap_hypo_LTR <- subset(mutant_overlap_hypo_unique, grepl("compo=RL",Compo))
  mutant_overlap_hypo_Copia <- subset(mutant_overlap_hypo_unique, grepl("RLC",Compo))
  mutant_overlap_hypo_Gypsy <- subset(mutant_overlap_hypo_unique, grepl("RLG",Compo))
  mutant_overlap_hypo_TIR <- subset(mutant_overlap_hypo_unique, grepl("compo=DT",Compo))
  
  #Separate into 6 most common TE families
  mutant_overlap_hypo_Angela <- subset(mutant_overlap_hypo_unique, grepl("RLC_famc1",Compo))
  mutant_overlap_hypo_Jorge <- subset(mutant_overlap_hypo_unique, grepl("DTC_famc2",Compo))
  mutant_overlap_hypo_Sabrina <- subset(mutant_overlap_hypo_unique, grepl("RLG_famc2",Compo))
  mutant_overlap_hypo_Fatima <- subset(mutant_overlap_hypo_unique, grepl("RLG_famc1",Compo))
  mutant_overlap_hypo_Sumana <- subset(mutant_overlap_hypo_unique, grepl("RLG_famc7",Compo))
  mutant_overlap_hypo_WHAM <- subset(mutant_overlap_hypo_unique, grepl("RLG_famc5",Compo))
  
  mutant_overlap_hyper <- read.table(paste(mutant,"_overlap_hyperDMRs_TEs.txt",sep=""), sep="\t")
  
  mutant_overlap_hyper <- separate_wider_delim(mutant_overlap_hyper, cols=V9, delim=";",
                                               names=c("ID","Name","Ontology_term","Compo","Status"),
                                               too_few="align_start", too_many="drop")
  mutant_overlap_hyper <- separate_wider_delim(mutant_overlap_hyper, cols=ID, delim="=",
                                               names=c("x","ID"))
  
  mutant_overlap_hyper$ID <- factor(mutant_overlap_hyper$ID)
  mutant_overlap_hyper_unique <- mutant_overlap_hyper %>%
    group_by(ID) %>%
    filter(row_number()==1)
  
  #Separate out by TE family
  mutant_overlap_hyper_LTR <- subset(mutant_overlap_hyper_unique, grepl("compo=RL",Compo))
  mutant_overlap_hyper_Copia <- subset(mutant_overlap_hyper_unique, grepl("RLC",Compo))
  mutant_overlap_hyper_Gypsy <- subset(mutant_overlap_hyper_unique, grepl("RLG",Compo))
  mutant_overlap_hyper_TIR <- subset(mutant_overlap_hyper_unique, grepl("compo=DT",Compo))
  
  #Separate into 6 most common TE families
  mutant_overlap_hyper_Angela <- subset(mutant_overlap_hyper_unique, grepl("RLC_famc1",Compo))
  mutant_overlap_hyper_Jorge <- subset(mutant_overlap_hyper_unique, grepl("DTC_famc2",Compo))
  mutant_overlap_hyper_Sabrina <- subset(mutant_overlap_hyper_unique, grepl("RLG_famc2",Compo))
  mutant_overlap_hyper_Fatima <- subset(mutant_overlap_hyper_unique, grepl("RLG_famc1",Compo))
  mutant_overlap_hyper_Sumana <- subset(mutant_overlap_hyper_unique, grepl("RLG_famc7",Compo))
  mutant_overlap_hyper_WHAM <- subset(mutant_overlap_hyper_unique, grepl("RLG_famc5",Compo))
  
  mutant_overlap_unique <- rbind(mutant_overlap_hyper_unique, mutant_overlap_hypo_unique) %>%
    group_by(ID) %>%
    filter(row_number()==1)
  
  mutant_number_TE_DMRs <- cbind(mutant, nrow(mutant_overlap_unique), nrow(mutant_overlap_hypo_unique), nrow(mutant_overlap_hyper_unique))
  number_TE_DMRs <- rbind(number_TE_DMRs, mutant_number_TE_DMRs)
  
  mutant_number_family_DMRs_hypo <- cbind(mutant, nrow(mutant_overlap_hypo_unique), nrow(mutant_overlap_hypo_LTR), 
                                          nrow(mutant_overlap_hypo_Copia), nrow(mutant_overlap_hypo_Gypsy), nrow(mutant_overlap_hypo_TIR))
  number_family_DMRs_hypo <- rbind(number_family_DMRs_hypo, mutant_number_family_DMRs_hypo)
  
  mutant_number_family_DMRs_hyper <- cbind(mutant, nrow(mutant_overlap_hyper_unique), nrow(mutant_overlap_hyper_LTR),
                                           nrow(mutant_overlap_hyper_Copia), nrow(mutant_overlap_hyper_Gypsy),nrow(mutant_overlap_hyper_TIR))
  number_family_DMRs_hyper <- rbind(number_family_DMRs_hyper, mutant_number_family_DMRs_hyper)
  
  mutant_number_subfamily_DMRs_hypo <- cbind(mutant, nrow(mutant_overlap_hypo_unique), 
                                             nrow(mutant_overlap_hypo_Angela), nrow(mutant_overlap_hypo_Jorge),
                                             nrow(mutant_overlap_hypo_Sabrina), nrow(mutant_overlap_hypo_Fatima),
                                             nrow(mutant_overlap_hypo_Sumana), nrow(mutant_overlap_hypo_WHAM))
  number_subfamily_DMRs_hypo <- rbind(number_subfamily_DMRs_hypo, mutant_number_subfamily_DMRs_hypo)
  
  mutant_number_subfamily_DMRs_hyper <- cbind(mutant, nrow(mutant_overlap_hyper_unique), 
                                             nrow(mutant_overlap_hyper_Angela), nrow(mutant_overlap_hyper_Jorge),
                                             nrow(mutant_overlap_hyper_Sabrina), nrow(mutant_overlap_hyper_Fatima),
                                             nrow(mutant_overlap_hyper_Sumana), nrow(mutant_overlap_hyper_WHAM))
  number_subfamily_DMRs_hyper <- rbind(number_subfamily_DMRs_hyper, mutant_number_subfamily_DMRs_hyper)
}

colnames(number_TE_DMRs) <- c("Mutant","Total TE-DMRs","TE_Hypo-DMRs","TE_Hyper-DMRs")

#Total number of TEs (excluding chrUn) is 853523
number_TE_DMRs <- data.frame(number_TE_DMRs)
number_TE_DMRs$Total.TE.DMRs <- as.numeric(number_TE_DMRs$Total.TE.DMRs)
number_TE_DMRs$TE_Hypo.DMRs <- as.numeric(number_TE_DMRs$TE_Hypo.DMRs)
number_TE_DMRs$TE_Hyper.DMRs <- as.numeric(number_TE_DMRs$TE_Hyper.DMRs)

number_TE_DMRs$Total_TE_percent <- (number_TE_DMRs$Total.TE.DMRs/853523)*100
number_TE_DMRs$Total_HypoTE_percent <- (number_TE_DMRs$TE_Hypo.DMRs/853523)*100
number_TE_DMRs$Total_HyperTE_percent <- (number_TE_DMRs$TE_Hyper.DMRs/853523)*100

number_TE_DMRs <- pivot_longer(number_TE_DMRs, cols=Total_HypoTE_percent:Total_HyperTE_percent,
                               names_to="Type",values_to="Percentage")
number_TE_DMRs$Mutant <- factor(number_TE_DMRs$Mutant, levels=c("A_single","B_single","D_single",
                                                                "AB_double","AD_double","BD_double",
                                                                "Aabbdd"))

ggplot(number_TE_DMRs, aes(x=Mutant, y=Percentage, fill=Type)) +
  geom_bar(position="stack",stat="identity", colour="black") +
  theme_bw() +
  theme(legend.position="none")

setwd("U:/Year1/met1 Mutants Project/WGBS_analysis/")
ggsave("percentage_TE_DMRs.pdf", width=5, height=5, units="in")
