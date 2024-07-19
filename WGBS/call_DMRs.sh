#!/bin/bash
#SBATCH -p jic-medium
#SBATCH -t 1-00:00
#SBATCH -c 4
#SBATCH --mem=120000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=delfi.dorussen@jic.ac.uk
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/DMR_calling/slurm_output/DMR-slurm-%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/DMR_calling/slurm_output/DMR-slurm-%j.err

singularity exec --overlay /jic/scratch/groups/Philippa-Borrill/scripts/R-3.6-all.img /jic/scratch/groups/Philippa-Borrill/scripts/r-upd.img Rscript  call_DMRs_all_chrom.R