#Delfi Dorussen
#Edit gff file with all HC genes (v1.1) to contain only reads annotated as 'gene'
#and change 'start' and 'end' to include +/- 1kb up/downstream

setwd("Z:/References/iwgsc_ref_seq_1.1/iwgsc_refseqv1.1_genes_2017July06/")

library(microseq)

HC_gff <- readGFF("IWGSC_v1.1_HC_20170706.gff3")
head(HC_gff)

HC_gff_genes <- subset(HC_gff, Type=="gene")
HC_gff_genes$Start <- HC_gff_genes$Start - 1000
HC_gff_genes$End <- HC_gff_genes$End + 1000
HC_gff_genes$Start <- ifelse(HC_gff_genes$Start<0, 1, HC_gff_genes$Start)

writeGFF(HC_gff_genes, "IWGSC_v1.1_HC_20170706_genes1kb.gff3")