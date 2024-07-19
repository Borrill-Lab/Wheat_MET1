#!/bin/bash
#SBATCH -p jic-medium
#SBATCH -t 2-00:00
#SBATCH -c 4
#SBATCH --mem=200000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=delfi.dorussen@jic.ac.uk
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/filter_by_chrom/chr_plots_difference/slurm_output/plot-slurm-%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/filter_by_chrom/chr_plots_difference/slurm_output/plot-slurm-%j.err

singularity exec --overlay /jic/scratch/groups/Philippa-Borrill/scripts/R-3.6-all-2.img /jic/scratch/groups/Philippa-Borrill/scripts/r-upd.img Rscript  plot_all_chrom_difference.R