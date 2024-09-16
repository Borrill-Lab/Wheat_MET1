#!/bin/bash
#SBATCH -p jic-medium
#SBATCH -t 1-00:00
#SBATCH -c 4
#SBATCH --mem=200000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=delfi.dorussen@jic.ac.uk
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_TE/hisat_index/slurm_output/index-slurm-%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_TE/hisat_index/slurm_output/index-slurm-%j.err

source package f9c1e0c5-d0e8-4ba0-9edd-88235400fa13

cd /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_TE/hisat_index/

hisat2-build /jic/scratch/groups/Philippa-Borrill/References/iwgsc_ref_seq_1.1/iwgsc_refseqv1.0_all_chromosomes/161010_Chinese_Spring_v1.0_pseudomolecules.fasta \
161010_Chinese_Spring_v1.0_pseudomolecules