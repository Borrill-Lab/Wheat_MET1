#Delfi Dorussen
#Aim to load in the TE count files for each of the met1 samples and combine into 
#one file

library(dplyr)

setwd("Z:/Delfi/MET1_TE/htseq_count/slurm_output/")

AABBDD_rep1 <- read.table("htseq-Unknown_BU702-001T0001-slurm-2776305.out")
AABBDD_rep1 <- AABBDD_rep1[1:911917,]
colnames(AABBDD_rep1) <- c("Feature","AABBDD_rep1")

AABBDD_rep2 <- read.table("htseq-Unknown_BU702-001T0002-slurm-2776306.out")
AABBDD_rep2 <- AABBDD_rep2[1:911917,]
colnames(AABBDD_rep2) <- c("Feature","AABBDD_rep2")

AABBDD_rep3 <- read.table("htseq-Unknown_BU702-001T0003-slurm-2776307.out")
AABBDD_rep3 <- AABBDD_rep3[1:911917,]
colnames(AABBDD_rep3) <- c("Feature","AABBDD_rep3")

aaBBDD_rep1 <- read.table("htseq-Unknown_BU702-001T0004-slurm-2776308.out")
aaBBDD_rep1 <- aaBBDD_rep1[1:911917,]
colnames(aaBBDD_rep1) <- c("Feature","aaBBDD_rep1")

aaBBDD_rep2 <- read.table("htseq-Unknown_BU702-001T0005-slurm-2776309.out")
aaBBDD_rep2 <- aaBBDD_rep2[1:911917,]
colnames(aaBBDD_rep2) <- c("Feature","aaBBDD_rep2")

aaBBDD_rep3 <- read.table("htseq-Unknown_BU702-001T0006-slurm-2776310.out")
aaBBDD_rep3 <- aaBBDD_rep3[1:911917,]
colnames(aaBBDD_rep3) <- c("Feature","aaBBDD_rep3")

AAbbDD_rep1 <- read.table("htseq-Unknown_BU702-001T0007-slurm-2776311.out")
AAbbDD_rep1 <- AAbbDD_rep1[1:911917,]
colnames(AAbbDD_rep1) <- c("Feature","AAbbDD_rep1")

AAbbDD_rep2 <- read.table("htseq-Unknown_BU702-001T0008-slurm-2776312.out")
AAbbDD_rep2 <- AAbbDD_rep2[1:911917,]
colnames(AAbbDD_rep2) <- c("Feature","AAbbDD_rep2")

AAbbDD_rep3 <- read.table("htseq-Unknown_BU702-001T0009-slurm-2776313.out")
AAbbDD_rep3 <- AAbbDD_rep3[1:911917,]
colnames(AAbbDD_rep3) <- c("Feature","AAbbDD_rep3")

AABBdd_rep1 <- read.table("htseq-Unknown_BU702-001T0010-slurm-2776314.out")
AABBdd_rep1 <- AABBdd_rep1[1:911917,]
colnames(AABBdd_rep1) <- c("Feature","AABBdd_rep1")

AABBdd_rep2 <- read.table("htseq-Unknown_BU702-001T0011-slurm-2776315.out")
AABBdd_rep2 <- AABBdd_rep2[1:911917,]
colnames(AABBdd_rep2) <- c("Feature","AABBdd_rep2")

AABBdd_rep3 <- read.table("htseq-Unknown_BU702-001T0012-slurm-2776316.out")
AABBdd_rep3 <- AABBdd_rep3[1:911917,]
colnames(AABBdd_rep3) <- c("Feature","AABBdd_rep3")

aabbDD_rep1 <- read.table("htseq-Unknown_BU702-001T0013-slurm-2776317.out")
aabbDD_rep1 <- aabbDD_rep1[1:911917,]
colnames(aabbDD_rep1) <- c("Feature","aabbDD_rep1")

aabbDD_rep2 <- read.table("htseq-Unknown_BU702-001T0014-slurm-2776318.out")
aabbDD_rep2 <- aabbDD_rep2[1:911917,]
colnames(aabbDD_rep2) <- c("Feature","aabbDD_rep2")

aabbDD_rep3 <- read.table("htseq-Unknown_BU702-001T0015-slurm-2776319.out")
aabbDD_rep3 <- aabbDD_rep3[1:911917,]
colnames(aabbDD_rep3) <- c("Feature","aabbDD_rep3")

aaBBdd_rep1 <- read.table("htseq-Unknown_BU702-001T0016-slurm-2776320.out")
aaBBdd_rep1 <- aaBBdd_rep1[1:911917,]
colnames(aaBBdd_rep1) <- c("Feature","aaBBdd_rep1")

aaBBdd_rep2 <- read.table("htseq-Unknown_BU702-001T0017-slurm-2776321.out")
aaBBdd_rep2 <- aaBBdd_rep2[1:911917,]
colnames(aaBBdd_rep2) <- c("Feature","aaBBdd_rep2")

aaBBdd_rep3 <- read.table("htseq-Unknown_BU702-001T0018-slurm-2776322.out")
aaBBdd_rep3 <- aaBBdd_rep3[1:911917,]
colnames(aaBBdd_rep3) <- c("Feature","aaBBdd_rep3")

AAbbdd_rep1 <- read.table("htseq-Unknown_BU702-001T0019-slurm-2776323.out")
AAbbdd_rep1 <- AAbbdd_rep1[1:911917,]
colnames(AAbbdd_rep1) <- c("Feature","AAbbdd_rep1")

AAbbdd_rep2 <- read.table("htseq-Unknown_BU702-001T0020-slurm-2776324.out")
AAbbdd_rep2 <- AAbbdd_rep2[1:911917,]
colnames(AAbbdd_rep2) <- c("Feature","AAbbdd_rep2")

AAbbdd_rep3 <- read.table("htseq-Unknown_BU702-001T0021-slurm-2776325.out")
AAbbdd_rep3 <- AAbbdd_rep3[1:911917,]
colnames(AAbbdd_rep3) <- c("Feature","AAbbdd_rep3")

Aabbdd_rep1 <- read.table("htseq-Unknown_BU702-001T0025-slurm-2776329.out")
Aabbdd_rep1 <- Aabbdd_rep1[1:911917,]
colnames(Aabbdd_rep1) <- c("Feature","Aabbdd_rep1")


all_counts <- AABBDD_rep1 %>%
  left_join(AABBDD_rep2, by="Feature") %>%
  left_join(AABBDD_rep3, by="Feature") %>%
  left_join(aaBBDD_rep1, by="Feature") %>%
  left_join(aaBBDD_rep2, by="Feature") %>%
  left_join(aaBBDD_rep3, by="Feature") %>%
  left_join(AAbbDD_rep1, by="Feature") %>%
  left_join(AAbbDD_rep2, by="Feature") %>%
  left_join(AAbbDD_rep3, by="Feature") %>%
  left_join(AABBdd_rep1, by="Feature") %>%
  left_join(AABBdd_rep2, by="Feature") %>%
  left_join(AABBdd_rep3, by="Feature") %>%
  left_join(aabbDD_rep1, by="Feature") %>%
  left_join(aabbDD_rep2, by="Feature") %>%
  left_join(aabbDD_rep3, by="Feature") %>%
  left_join(aaBBdd_rep1, by="Feature") %>%
  left_join(aaBBdd_rep2, by="Feature") %>%
  left_join(aaBBdd_rep3, by="Feature") %>%
  left_join(AAbbdd_rep1, by="Feature") %>%
  left_join(AAbbdd_rep2, by="Feature") %>%
  left_join(AAbbdd_rep3, by="Feature") %>%
  left_join(Aabbdd_rep1, by="Feature")

setwd("U:/Year1/met1 Mutants Project/TE_RNA_seq_analysis/")
write.csv(all_counts, "TE_counts.csv", quote=FALSE, row.names=FALSE)

