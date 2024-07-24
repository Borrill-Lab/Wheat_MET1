# Wheat_MET1

## Data
Raw read files for RNA-seq and Whole Genome Bisulfite Sequencing (WGBS) data can be found in the ENA with accession numbers PRJEB77425 and PRJEB77426 respectively.  

## Analysis of RNA-seq Data  

All scripts for the RNA-seq analysis are found in the [RNA-seq](https://github.com/Borrill-Lab/Wheat_MET1/tree/main/RNA-seq) folder.  

- **Trim**: fastp.pl
- **Pseudoalignment & quantification with kallisto**: kallisto.pl
- **Import quantifications and summarise at gene level**: import_kallisto_quantifications.R
- **Identify differentially expressed genes (DEGs) in *met1* mutants**: DESeq_analysis_MET1_vs_AABBDD.R
- **Plot the number of DEGs in each of the *met1* mutants**: number_of_DEGs_plot.R (produces Figures 4A-B)
- **Plot particularly interesting genes**: plot_genes_of_interest.R (produces Figure 4C)
- **Gene Ontology enrichment for DEGs in Aabbdd mutant**: GO_enrichment_MET1_vs_AABBDD.R (produces Figures 4D-E)

## Analysis of Whole Genome Bisulfite Sequencing (WGBS) Data

All scripts for the WGBS analysis are found in the [WGBS](https://github.com/Borrill-Lab/Wheat_MET1/tree/main/WGBS) folder.  
  
- **Trim**: fastp.pl
- **Create bisulfite-converted reference**: bismark_index.sh, requires IWGSC v1.0 genome, available [here](https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Assemblies/v1.0/)  
- **Map to bisulfite-converted reference**: bismark.pl  
- **Extract methylation data and produce CX_report**: bismark_methylation_extraction.pl  
- **Split CX_report files into separate chromosomes**: filter_all_chromosomes.sh
- **Split CX_report files by sequence context and convert to bigwig files**: CXreport_to_bigwig.R run via CXreport_to_bigwig.sh 
- **Convert gene and transposon gff3 files to gtf files and filter for chromosome of interest**: gff_to_gtf.sh and gff_to_gtf_TE.sh, requires IWGSC v1.1 HC genes gff3 file, available [here](https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v1.1/) and IWGSC Transposable Elements gff3 file, available [here](https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v1.0/)
- **Plot average methylation across genes and TEs**: computeMatrix_5D.sh and plotProfile_5D.sh (produces Figures 3B-C and Supplementary Figure 4) 
- **Produce methylation profile plots**: plot_all_chrom.R run via plot_all_chrom.sh (produces Figures 3D-F and Supplementary Figures 6-12)  
- **Plot absolute difference in methylation between WT and *met1* mutants**: plot_all_chrom_difference.R run via plot_all_chrom_difference.sh (produces Supplementary Figure 5)  
- **Plot methylation profile plots for all single and double mutants**: plot_all_chrom_singles.R run via plot_all_chrom_singles.sh, and plot_all_chrom_doubles.R run via plot_all_chrom_doubles.sh (produces Supplementary Figure 12)  
- **Call CG-DMRs**: call_DMRs_all_chrom.R run via call_DMRs.sh  
- **Combine DMRs from all chromosomes, plot number of DMRs, and density of DMRs in chromosome regions**: combine_all_DMRs.R (produces Figures 3G-H and Supplementary Figure 13)  
- **Convert DMR files into bed files**: DMRs_to_bed.R  
- **Filter gff file with genes and add 1kb up- and downstream**: gff_to_genes1kb.R, requires IWGSC v1.1 HC genes gff3 file, available [here](https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v1.1/)   
- **Find overlap between genes and DMRs**: DMRs_intersect_genes.sh  
- **Plot number of genes overlapping DMRs and identify overlap with DEGs**: DMR_gene_overlap.R (produces Figures 5A,C-D)  
- **Filter gff file with TEs**: TEgff_to_completeTE.R, requires IWGSC Transposable Elements gff3 file, available [here](https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Annotations/v1.0/)  
- **Find overlap between TEs and DMRs**: DMRs_intersect_TEs.sh  
- **Plot the number of TEs overlapping DMRs**: DMR_TE_overlap.R (produces Figure 5B)
