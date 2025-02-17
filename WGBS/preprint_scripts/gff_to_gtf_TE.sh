#!/bin/bash
#SBATCH -p jic-medium
#SBATCH -t 1-00:00
#SBATCH -c 4
#SBATCH --mem=200000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=delfi.dorussen@jic.ac.uk
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/DMR_overlap/slurm_output/gffread-slurm-%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/DMR_overlap/slurm_output/gffread-slurm-%j.err

source package  /tgac/software/testing/bin/gffread-0.11.4

cd /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/DMR_overlap/

grep "chr5D" iwgsc_refseqv1.0_TransposableElements_2017Mar13.gff3 > iwgsc_refseqv1.0_TransposableElements_2017Mar13_5D.gff3
grep "status=complete" iwgsc_refseqv1.0_TransposableElements_2017Mar13_5D.gff3 > iwgsc_refseqv1.0_TransposableElements_2017Mar13_complete_5D.gff3
grep "repeat_region" iwgsc_refseqv1.0_TransposableElements_2017Mar13_complete_5D.gff3 > iwgsc_refseqv1.0_TransposableElements_2017Mar13_complete_5D_repeatregion.gff3
sed -i 's/repeat_region/mRNA/g' iwgsc_refseqv1.0_TransposableElements_2017Mar13_complete_5D_repeatregion.gff3

gffread iwgsc_refseqv1.0_TransposableElements_2017Mar13_complete_5D_repeatregion.gff3 -T -o iwgsc_refseqv1.0_TransposableElements_2017Mar13_complete_5D.gtf
