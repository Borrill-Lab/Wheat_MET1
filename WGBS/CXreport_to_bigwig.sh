#!/bin/bash
#SBATCH -p jic-medium
#SBATCH -t 1-00:00
#SBATCH -c 4
#SBATCH --mem=120000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=delfi.dorussen@jic.ac.uk
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/filter_by_chrom/bigwig_files/slurm_output/bigwig-slurm-%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/filter_by_chrom/bigwig_files/slurm_output/bigwig-slurm-%j.err

singularity exec --overlay R-3.6-all-2.img r-upd.img Rscript  CXreport_to_bigwig.R