#Delfi Dorussen
#Edit gff file with TE positions and filter to keep only those with status=complete
#and only keep repeat_region TEs (not nested repeats as these overlap one another)

setwd("Z:/Delfi/MET1_WGBS/DMR_overlap/")

library(microseq)

TE_gff <- readGFF("iwgsc_refseqv1.0_TransposableElements_2017Mar13.gff3")
head(TE_gff)

TE_gff_complete <- subset(TE_gff, grepl("status=complete",Attributes))
TE_gff_complete <- subset(TE_gff_complete, Type=="repeat_region")

writeGFF(TE_gff_complete, "iwgsc_refseqv1.0_TransposableElements_2017Mar13_complete.gff3")