#!/bin/bash
#SBATCH -p jic-medium
#SBATCH -t 1-00:00
#SBATCH -c 4
#SBATCH --mem=200000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=delfi.dorussen@jic.ac.uk
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/methylation_extraction/slurm_output/gffread-slurm-%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/methylation_extraction/slurm_output/gffread-slurm-%j.err

source package  /tgac/software/testing/bin/gffread-0.11.4

cd /jic/scratch/groups/Philippa-Borrill/References/iwgsc_ref_seq_1.1/iwgsc_refseqv1.1_genes_2017July06/

gffread IWGSC_v1.1_HC_20170706.gff3 -T -o IWGSC_v1.1_HC_20170706_v2.gtf

grep "chr5D" IWGSC_v1.1_HC_20170706_v2.gtf > IWGSC_v1.1_HC_20170706_v2_5D.gtf
