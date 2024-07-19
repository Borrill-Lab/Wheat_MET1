#Delfi Dorussen
#Aim to combine DMR files for all chromosomes for each mutant and compare number
#of DMRs in each (and whether they are hypo- or hypermethylated)
#Also plot the number of DMRs for each of the mutants

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggforce)

setwd("Z:/Delfi/MET1_WGBS/DMR_calling/")

mutants <- c("A_single","B_single","D_single","AB_double","AD_double","BD_double","Aabbdd")

all_mutants_by_chrom <- NULL
all_mutants <- NULL

for(i in 1:7){
  mutant <- mutants[i]
  
  WT_mutant_DMRs_Chr1A <- read.table(paste("WT_",mutant,"_DMRs_Chr1A.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr1B <- read.table(paste("WT_",mutant,"_DMRs_Chr1B.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr1D <- read.table(paste("WT_",mutant,"_DMRs_Chr1D.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr2A <- read.table(paste("WT_",mutant,"_DMRs_Chr2A.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr2B <- read.table(paste("WT_",mutant,"_DMRs_Chr2B.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr2D <- read.table(paste("WT_",mutant,"_DMRs_Chr2D.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr3A <- read.table(paste("WT_",mutant,"_DMRs_Chr3A.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr3B <- read.table(paste("WT_",mutant,"_DMRs_Chr3B.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr3D <- read.table(paste("WT_",mutant,"_DMRs_Chr3D.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr4A <- read.table(paste("WT_",mutant,"_DMRs_Chr4A.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr4B <- read.table(paste("WT_",mutant,"_DMRs_Chr4B.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr4D <- read.table(paste("WT_",mutant,"_DMRs_Chr4D.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr5A <- read.table(paste("WT_",mutant,"_DMRs_Chr5A.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr5B <- read.table(paste("WT_",mutant,"_DMRs_Chr5B.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr5D <- read.table(paste("WT_",mutant,"_DMRs_Chr5D.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr6A <- read.table(paste("WT_",mutant,"_DMRs_Chr6A.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr6B <- read.table(paste("WT_",mutant,"_DMRs_Chr6B.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr6D <- read.table(paste("WT_",mutant,"_DMRs_Chr6D.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr7A <- read.table(paste("WT_",mutant,"_DMRs_Chr7A.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr7B <- read.table(paste("WT_",mutant,"_DMRs_Chr7B.txt",sep=""), header=TRUE)
  WT_mutant_DMRs_Chr7D <- read.table(paste("WT_",mutant,"_DMRs_Chr7D.txt",sep=""), header=TRUE)
  
  
  WT_mutant_DMRs <- rbind(WT_mutant_DMRs_Chr1A,WT_mutant_DMRs_Chr1B,WT_mutant_DMRs_Chr1D,
                            WT_mutant_DMRs_Chr2A,WT_mutant_DMRs_Chr2B,WT_mutant_DMRs_Chr2D,
                            WT_mutant_DMRs_Chr3A,WT_mutant_DMRs_Chr3B,WT_mutant_DMRs_Chr3D,
                            WT_mutant_DMRs_Chr4A,WT_mutant_DMRs_Chr4B,WT_mutant_DMRs_Chr4D,
                            WT_mutant_DMRs_Chr5A,WT_mutant_DMRs_Chr5B,WT_mutant_DMRs_Chr5D,
                            WT_mutant_DMRs_Chr6A,WT_mutant_DMRs_Chr6B,WT_mutant_DMRs_Chr6D,
                            WT_mutant_DMRs_Chr7A,WT_mutant_DMRs_Chr7B,WT_mutant_DMRs_Chr7D)
  
  WT_mutant_DMRs_sum_chrom <- WT_mutant_DMRs %>%
    group_by(seqnames,regionType) %>%
    summarise(count=n())
  WT_mutant_DMRs_sum_chrom$Sample <- mutant
  all_mutants_by_chrom <- rbind(all_mutants_by_chrom,WT_mutant_DMRs_sum_chrom)
  
  WT_mutant_DMRs_sum <- WT_mutant_DMRs %>%
    group_by(regionType) %>%
    summarise(count=n())
  WT_mutant_DMRs_sum$Sample <- mutant
  all_mutants <- rbind(all_mutants,WT_mutant_DMRs_sum)
}


all_mutants$Sample <- factor(all_mutants$Sample, levels=c("A_single","B_single","D_single",
                                                          "AB_double","AD_double","BD_double","Aabbdd"))
all_mutants$regionType <- factor(all_mutants$regionType, levels=c("gain","loss"))

ggplot(all_mutants, aes(x=Sample, y=count, fill=regionType)) +
  geom_bar(position="stack", stat="identity", colour="black") +
  theme_bw() +
  scale_fill_manual(values=c("#a4c7ff","#ff6643"), labels=c("Gain","Loss")) +
  ylab("Number of DMRs") +
  facet_zoom(ylim = c(0, 300000)) +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
ggsave("number_of_DMRs_AABBDD.pdf", height=6, width=9, units="in")



###Calculating number of DMRs in each chromosome region
chromosome_partitions <- read.csv("region_partition.csv")

chromosome_partitions$Length <- chromosome_partitions$Length * 10^6
chromosome_partitions$R1_R2a <- chromosome_partitions$R1_R2a * 10^6
chromosome_partitions$R2a_C <- chromosome_partitions$R2a_C * 10^6
chromosome_partitions$C_R2b <- chromosome_partitions$C_R2b * 10^6
chromosome_partitions$R2b_R3 <- chromosome_partitions$R2b_R3 * 10^6

chromosomes <- c("Chr1A","Chr1B","Chr1D","Chr2A","Chr2B","Chr2D","Chr3A","Chr3B","Chr3D",
                 "Chr4A","Chr4B","Chr4D","Chr5A","Chr5B","Chr5D","Chr6A","Chr6B","Chr6D",
                 "Chr7A","Chr7B","Chr7D")

chromosomes_lc <- c("chr1A","chr1B","chr1D","chr2A","chr2B","chr2D","chr3A","chr3B","chr3D",
                    "chr4A","chr4B","chr4D","chr5A","chr5B","chr5D","chr6A","chr6B","chr6D",
                    "chr7A","chr7B","chr7D")


WT_mutant_DMRs_regions_all <- NULL

for(k in 1:7) {
  WT_mutant_DMRs_regions <- NULL
  
  mutant <- mutants[k]
  
  for(j in 1:21){
    chromosome <- chromosomes[j]
    chromosome_lc <- chromosomes_lc[j]
    
    chromosome_partitions_chr <- subset(chromosome_partitions, Chr==chromosome_lc)
    
    WT_mutant_DMRs <- read.table(paste("WT_",mutant,"_DMRs_",chromosome,".txt",sep=""), header=TRUE)
    for(i in 1:nrow(WT_mutant_DMRs)){
      if(WT_mutant_DMRs$start[i]<chromosome_partitions_chr$R1_R2a[1]){
        WT_mutant_DMRs$partition[i] <- "R1"
      } else if(WT_mutant_DMRs$start[i]>=chromosome_partitions_chr$R1_R2a[1]&WT_mutant_DMRs$start[i]<chromosome_partitions_chr$R2a_C[1]){
        WT_mutant_DMRs$partition[i] <- "R2a"
      } else if(WT_mutant_DMRs$start[i]>=chromosome_partitions_chr$R2a_C[1]&WT_mutant_DMRs$start[i]<chromosome_partitions_chr$C_R2b[1]){
        WT_mutant_DMRs$partition[i] <- "C"
      } else if(WT_mutant_DMRs$start[i]>=chromosome_partitions_chr$C_R2b[1]&WT_mutant_DMRs$start[i]<chromosome_partitions_chr$R2b_R3[1]){
        WT_mutant_DMRs$partition[i] <- "R2b"
      } else if(WT_mutant_DMRs$start[i]>=chromosome_partitions_chr$R2b_R3[1]){
        WT_mutant_DMRs$partition[i] <- "R3"
      } else{
        WT_mutant_DMRs$partition[i] <- "NA"
      }
    }
    
    WT_mutant_DMRs_chr_sum <- WT_mutant_DMRs %>%
      group_by(partition, regionType) %>%
      summarise(count=n())
    
    WT_mutant_DMRs_chr_sum$Chr <- chromosome
    
    WT_mutant_DMRs_regions <- rbind(WT_mutant_DMRs_regions, WT_mutant_DMRs_chr_sum)
  }
  
  WT_mutant_DMRs_regions_sum <- WT_mutant_DMRs_regions %>%
    group_by(partition, regionType) %>%
    summarise(total=sum(count))
  WT_mutant_DMRs_regions_sum$Genotype <- mutant
  
  WT_mutant_DMRs_regions_all <- rbind(WT_mutant_DMRs_regions_all, WT_mutant_DMRs_regions_sum)
  
}

WT_mutant_DMRs_regions_all$partition <- factor(WT_mutant_DMRs_regions_all$partition,
                                               levels=c("R1","R2a","C","R2b","R3"))

ggplot(subset(WT_mutant_DMRs_regions_all, regionType=="loss"), aes(x=Genotype, y=total)) +
  geom_point() +
  facet_wrap(vars(partition), nrow=1)

ggplot(subset(WT_mutant_DMRs_regions_all, regionType=="loss"), aes(x=Genotype, y=total, fill=partition)) +
  geom_bar(position="dodge", stat="identity")

#Correct for length of the regions
chromosome_partitions$R2a_length <- chromosome_partitions$R2a_C - chromosome_partitions$R1_R2a
chromosome_partitions$C_length <- chromosome_partitions$C_R2b - chromosome_partitions$R2a_C
chromosome_partitions$R2b_length <- chromosome_partitions$R2b_R3 - chromosome_partitions$C_R2b
chromosome_partitions$R3_length <- chromosome_partitions$Length - chromosome_partitions$R2b_R3

R1_total <- sum(chromosome_partitions$R1_R2a)
R2a_total <- sum(chromosome_partitions$R2a_length)
C_total <- sum(chromosome_partitions$C_length)
R2b_total <- sum(chromosome_partitions$R2b_length)
R3_total <- sum(chromosome_partitions$R3_length)

WT_mutant_DMRs_regions_all$total_corrected <- ifelse(WT_mutant_DMRs_regions_all$partition=="R1", WT_mutant_DMRs_regions_all$total/R1_total, NA)
WT_mutant_DMRs_regions_all$total_corrected <- ifelse(WT_mutant_DMRs_regions_all$partition=="R2a", WT_mutant_DMRs_regions_all$total/R2a_total, WT_mutant_DMRs_regions_all$total_corrected)
WT_mutant_DMRs_regions_all$total_corrected <- ifelse(WT_mutant_DMRs_regions_all$partition=="C", WT_mutant_DMRs_regions_all$total/C_total, WT_mutant_DMRs_regions_all$total_corrected)
WT_mutant_DMRs_regions_all$total_corrected <- ifelse(WT_mutant_DMRs_regions_all$partition=="R2b", WT_mutant_DMRs_regions_all$total/R2b_total, WT_mutant_DMRs_regions_all$total_corrected)
WT_mutant_DMRs_regions_all$total_corrected <- ifelse(WT_mutant_DMRs_regions_all$partition=="R3", WT_mutant_DMRs_regions_all$total/R3_total, WT_mutant_DMRs_regions_all$total_corrected)


WT_mutant_DMRs_regions_all$Genotype <- factor(WT_mutant_DMRs_regions_all$Genotype, 
                                              levels=c("A_single","B_single","D_single",
                                                       "AB_double","AD_double","BD_double","Aabbdd"))
WT_mutant_DMRs_regions_all$partition <- factor(WT_mutant_DMRs_regions_all$partition,
                                               levels=c("R1","R2a","C","R2b","R3"))

colours <- c("#1d4c89","#36b6c1","#fbf9d0","#36b6c1","#1d4c89")

ggplot(subset(WT_mutant_DMRs_regions_all, regionType=="loss"), aes(x=Genotype, y=total_corrected, fill=partition)) +
  geom_bar(position="dodge", stat="identity", colour="black") +
  theme_bw() +
  scale_fill_manual(values=colours)
ggsave("Hypo-DMRs_by_partition.pdf",height=6, width=8, units="in")

WT_mutant_DMRs_higher <- subset(WT_mutant_DMRs_regions_all, Genotype!="A_single"&Genotype!="B_single"&Genotype!="D_single")
ggplot(subset(WT_mutant_DMRs_higher, regionType=="loss"), aes(x=Genotype, y=total_corrected, fill=partition)) +
  geom_bar(position="dodge", stat="identity", colour="black") +
  theme_bw() +
  scale_fill_manual(values=colours)
ggsave("Hypo-DMRs_by_partition_double_Aabbdd.pdf",height=6, width=6, units="in")


ggplot(subset(WT_mutant_DMRs_regions_all, regionType=="gain"), aes(x=Genotype, y=total_corrected, fill=partition)) +
  geom_bar(position="dodge", stat="identity",colour="black") +
  theme_bw() +
  scale_fill_manual(values=colours)
ggsave("Hyper-DMRs_by_partition.pdf",height=6, width=8, units="in")

ggplot(subset(WT_mutant_DMRs_higher, regionType=="gain"), aes(x=Genotype, y=total_corrected, fill=partition)) +
  geom_bar(position="dodge", stat="identity",colour="black") +
  theme_bw() +
  scale_fill_manual(values=colours)
ggsave("Hyper-DMRs_by_partition_double_Aabbdd.pdf",height=6, width=6, units="in")


##Plotting number of DMRs for each chromosome
all_mutants_by_chrom_single <- subset(all_mutants_by_chrom,Sample=="A_single"|Sample=="B_single"|Sample=="D_single")

all_mutants_by_chrom_single$Sample <- factor(all_mutants_by_chrom_single$Sample, 
                                              levels=c("A_single","B_single","D_single"))


ggplot(all_mutants_by_chrom_single, aes(x=Sample, y=count, fill=regionType)) +
  geom_bar(position="stack", stat="identity", colour="black") +
  theme_bw() +
  scale_fill_manual(values=c("#a4c7ff","#ff6643"), labels=c("Gain","Loss")) +
  ylab("Number of DMRs") +
  facet_wrap(vars(seqnames),ncol=3) +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))

#Just looking at Chromosome 2A,2B,2D
ggplot(subset(all_mutants_by_chrom_single,seqnames=="chr2A"|seqnames=="chr2B"|seqnames=="chr2D"), 
       aes(x=Sample, y=count, fill=regionType)) +
  geom_bar(position="stack", stat="identity", colour="black") +
  theme_bw() +
  scale_fill_manual(values=c("#a4c7ff","#ff6643"), labels=c("Gain","Loss")) +
  ylab("Number of DMRs") +
  facet_wrap(vars(seqnames),ncol=3) +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
ggsave("single_mutants_group2.pdf", height=4, width=6, units="in")
