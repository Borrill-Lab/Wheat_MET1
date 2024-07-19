# Wheat_MET1

## Analysis of Whole Genome Bisulfite Sequencing (WGBS) Data

All scripts for the WGBS analysis are found in the 'WGBS' folder.  
  
**Trim**: fastp.pl  
**Map to bisulfite-converted reference**: bismark.pl  
**Extract methylation data and produce CX_report**: bismark_methylation_extraction.pl  
**Split CX_report files into separate chromosomes**: filter_all_chromosomes.sh  
**Produce methylation profile plots**: plot_all_chrom.R run via plot_all_chrom.sh  
- Produces Figures 3D-F and Supplementary Figures 5-11
**Plot absolute difference in methylation between WT and *met1* mutants**:
