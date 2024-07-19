# Wheat_MET1

## Analysis of Whole Genome Bisulfite Sequencing (WGBS) Data

All scripts for the WGBS analysis are found in the 'WGBS' folder.  
  
- **Trim**: fastp.pl  
- **Map to bisulfite-converted reference**: bismark.pl  
- **Extract methylation data and produce CX_report**: bismark_methylation_extraction.pl  
- **Split CX_report files into separate chromosomes**: filter_all_chromosomes.sh  
- **Produce methylation profile plots**: plot_all_chrom.R run via plot_all_chrom.sh (produces Figures 3D-F and Supplementary Figures 5-11)  
- **Plot absolute difference in methylation between WT and *met1* mutants**: plot_all_chrom_difference.R run via plot_all_chrom_difference.sh (produces Supplementary Figure 4)  
- **Plot methylation profile plots for all single and double mutants**: plot_all_chrom_singles.R run via plot_all_chrom_singles.sh, and plot_all_chrom_doubles.R run via plot_all_chrom_doubles.sh (produces Supplementary Figure 12)  
- **Call CG-DMRs**: call_DMRs_all_chrom.R run via call_DMRs.sh  
- **Combine DMRs from all chromosomes, plot number of DMRs, and density of DMRs in chromosome regions**: combine_all_DMRs.R (produces Figures 3G-H and Supplementary Figure 13)  
- **Convert DMR files into bed files**: DMRs_to_bed.R  
- **Filter gff file with genes and add 1kb up- and downstream**: gff_to_genes1kb.R  
- **Find overlap between genes and DMRs**: DMRs_intersect_genes.sh  
- **Plot number of genes overlapping DMRs and identify overlap with DEGs**: DMR_gene_overlap.R (produces Figures 5A,C-D)  
- **Filter gff file with TEs**: TEgff_to_completeTE.R  
- **Find overlap between TEs and DMRs**: DMRs_intersect_TEs.sh  
- **Plot the number of TEs overlapping DMRs**: DMR_TE_overlap.R (produces Figure 5B)
