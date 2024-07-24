#!/bin/bash
#SBATCH -p jic-medium
#SBATCH -t 2-00:00
#SBATCH -c 4
#SBATCH --mem=120000
#SBATCH --mail-type=none
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Marek/wgbs_Paragon_Charger/bismark-mapping/bis-slurm-%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Marek/wgbs_Paragon_Charger/bismark-mapping/bis-slurm-%j.err

source package 33c48798-0827-4add-8153-909c1bd83e89


bismark_genome_preparation /jic/scratch/groups/Philippa-Borrill/iwgsc_ref_seq_1.1/BSseq-ref


