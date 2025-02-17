#!/bin/bash
#SBATCH -p jic-medium
#SBATCH -t 1-00:00
#SBATCH -c 4
#SBATCH --mem=200000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=delfi.dorussen@jic.ac.uk
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/filter_by_chrom/bigwig_files/slurm_output/bigwigmerge-slurm-%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/filter_by_chrom/bigwig_files/slurm_output/bigwigmerge-slurm-%j.err

source package 71de5a7a-135b-417a-8de1-ede16dc52660

cd /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/filter_by_chrom/bigwig_files/

bigWigMerge \
Z0001_Chr1A_CHH.bw Z0001_Chr1B_CHH.bw Z0001_Chr1D_CHH.bw \
Z0001_Chr2A_CHH.bw Z0001_Chr2B_CHH.bw Z0001_Chr2D_CHH.bw \
Z0001_Chr3A_CHH.bw Z0001_Chr3B_CHH.bw Z0001_Chr3D_CHH.bw \
Z0001_Chr4A_CHH.bw Z0001_Chr4B_CHH.bw Z0001_Chr4D_CHH.bw \
Z0001_Chr5A_CHH.bw Z0001_Chr5B_CHH.bw Z0001_Chr5D_CHH.bw \
Z0001_Chr6A_CHH.bw Z0001_Chr6B_CHH.bw Z0001_Chr6D_CHH.bw \
Z0001_Chr7A_CHH.bw Z0001_Chr7B_CHH.bw Z0001_Chr7D_CHH.bw \
Z0001_all_chrom_CHH.bg

bigWigMerge \
Z0002_Chr1A_CHH.bw Z0002_Chr1B_CHH.bw Z0002_Chr1D_CHH.bw \
Z0002_Chr2A_CHH.bw Z0002_Chr2B_CHH.bw Z0002_Chr2D_CHH.bw \
Z0002_Chr3A_CHH.bw Z0002_Chr3B_CHH.bw Z0002_Chr3D_CHH.bw \
Z0002_Chr4A_CHH.bw Z0002_Chr4B_CHH.bw Z0002_Chr4D_CHH.bw \
Z0002_Chr5A_CHH.bw Z0002_Chr5B_CHH.bw Z0002_Chr5D_CHH.bw \
Z0002_Chr6A_CHH.bw Z0002_Chr6B_CHH.bw Z0002_Chr6D_CHH.bw \
Z0002_Chr7A_CHH.bw Z0002_Chr7B_CHH.bw Z0002_Chr7D_CHH.bw \
Z0002_all_chrom_CHH.bg

bigWigMerge \
Z0003_Chr1A_CHH.bw Z0003_Chr1B_CHH.bw Z0003_Chr1D_CHH.bw \
Z0003_Chr2A_CHH.bw Z0003_Chr2B_CHH.bw Z0003_Chr2D_CHH.bw \
Z0003_Chr3A_CHH.bw Z0003_Chr3B_CHH.bw Z0003_Chr3D_CHH.bw \
Z0003_Chr4A_CHH.bw Z0003_Chr4B_CHH.bw Z0003_Chr4D_CHH.bw \
Z0003_Chr5A_CHH.bw Z0003_Chr5B_CHH.bw Z0003_Chr5D_CHH.bw \
Z0003_Chr6A_CHH.bw Z0003_Chr6B_CHH.bw Z0003_Chr6D_CHH.bw \
Z0003_Chr7A_CHH.bw Z0003_Chr7B_CHH.bw Z0003_Chr7D_CHH.bw \
Z0003_all_chrom_CHH.bg

bigWigMerge \
Z0004_Chr1A_CHH.bw Z0004_Chr1B_CHH.bw Z0004_Chr1D_CHH.bw \
Z0004_Chr2A_CHH.bw Z0004_Chr2B_CHH.bw Z0004_Chr2D_CHH.bw \
Z0004_Chr3A_CHH.bw Z0004_Chr3B_CHH.bw Z0004_Chr3D_CHH.bw \
Z0004_Chr4A_CHH.bw Z0004_Chr4B_CHH.bw Z0004_Chr4D_CHH.bw \
Z0004_Chr5A_CHH.bw Z0004_Chr5B_CHH.bw Z0004_Chr5D_CHH.bw \
Z0004_Chr6A_CHH.bw Z0004_Chr6B_CHH.bw Z0004_Chr6D_CHH.bw \
Z0004_Chr7A_CHH.bw Z0004_Chr7B_CHH.bw Z0004_Chr7D_CHH.bw \
Z0004_all_chrom_CHH.bg

bigWigMerge \
Z0005_Chr1A_CHH.bw Z0005_Chr1B_CHH.bw Z0005_Chr1D_CHH.bw \
Z0005_Chr2A_CHH.bw Z0005_Chr2B_CHH.bw Z0005_Chr2D_CHH.bw \
Z0005_Chr3A_CHH.bw Z0005_Chr3B_CHH.bw Z0005_Chr3D_CHH.bw \
Z0005_Chr4A_CHH.bw Z0005_Chr4B_CHH.bw Z0005_Chr4D_CHH.bw \
Z0005_Chr5A_CHH.bw Z0005_Chr5B_CHH.bw Z0005_Chr5D_CHH.bw \
Z0005_Chr6A_CHH.bw Z0005_Chr6B_CHH.bw Z0005_Chr6D_CHH.bw \
Z0005_Chr7A_CHH.bw Z0005_Chr7B_CHH.bw Z0005_Chr7D_CHH.bw \
Z0005_all_chrom_CHH.bg

bigWigMerge \
Z0006_Chr1A_CHH.bw Z0006_Chr1B_CHH.bw Z0006_Chr1D_CHH.bw \
Z0006_Chr2A_CHH.bw Z0006_Chr2B_CHH.bw Z0006_Chr2D_CHH.bw \
Z0006_Chr3A_CHH.bw Z0006_Chr3B_CHH.bw Z0006_Chr3D_CHH.bw \
Z0006_Chr4A_CHH.bw Z0006_Chr4B_CHH.bw Z0006_Chr4D_CHH.bw \
Z0006_Chr5A_CHH.bw Z0006_Chr5B_CHH.bw Z0006_Chr5D_CHH.bw \
Z0006_Chr6A_CHH.bw Z0006_Chr6B_CHH.bw Z0006_Chr6D_CHH.bw \
Z0006_Chr7A_CHH.bw Z0006_Chr7B_CHH.bw Z0006_Chr7D_CHH.bw \
Z0006_all_chrom_CHH.bg

bigWigMerge \
Z0007_Chr1A_CHH.bw Z0007_Chr1B_CHH.bw Z0007_Chr1D_CHH.bw \
Z0007_Chr2A_CHH.bw Z0007_Chr2B_CHH.bw Z0007_Chr2D_CHH.bw \
Z0007_Chr3A_CHH.bw Z0007_Chr3B_CHH.bw Z0007_Chr3D_CHH.bw \
Z0007_Chr4A_CHH.bw Z0007_Chr4B_CHH.bw Z0007_Chr4D_CHH.bw \
Z0007_Chr5A_CHH.bw Z0007_Chr5B_CHH.bw Z0007_Chr5D_CHH.bw \
Z0007_Chr6A_CHH.bw Z0007_Chr6B_CHH.bw Z0007_Chr6D_CHH.bw \
Z0007_Chr7A_CHH.bw Z0007_Chr7B_CHH.bw Z0007_Chr7D_CHH.bw \
Z0007_all_chrom_CHH.bg

bigWigMerge \
Z0008_Chr1A_CHH.bw Z0008_Chr1B_CHH.bw Z0008_Chr1D_CHH.bw \
Z0008_Chr2A_CHH.bw Z0008_Chr2B_CHH.bw Z0008_Chr2D_CHH.bw \
Z0008_Chr3A_CHH.bw Z0008_Chr3B_CHH.bw Z0008_Chr3D_CHH.bw \
Z0008_Chr4A_CHH.bw Z0008_Chr4B_CHH.bw Z0008_Chr4D_CHH.bw \
Z0008_Chr5A_CHH.bw Z0008_Chr5B_CHH.bw Z0008_Chr5D_CHH.bw \
Z0008_Chr6A_CHH.bw Z0008_Chr6B_CHH.bw Z0008_Chr6D_CHH.bw \
Z0008_Chr7A_CHH.bw Z0008_Chr7B_CHH.bw Z0008_Chr7D_CHH.bw \
Z0008_all_chrom_CHH.bg

bigWigMerge \
Z0009_Chr1A_CHH.bw Z0009_Chr1B_CHH.bw Z0009_Chr1D_CHH.bw \
Z0009_Chr2A_CHH.bw Z0009_Chr2B_CHH.bw Z0009_Chr2D_CHH.bw \
Z0009_Chr3A_CHH.bw Z0009_Chr3B_CHH.bw Z0009_Chr3D_CHH.bw \
Z0009_Chr4A_CHH.bw Z0009_Chr4B_CHH.bw Z0009_Chr4D_CHH.bw \
Z0009_Chr5A_CHH.bw Z0009_Chr5B_CHH.bw Z0009_Chr5D_CHH.bw \
Z0009_Chr6A_CHH.bw Z0009_Chr6B_CHH.bw Z0009_Chr6D_CHH.bw \
Z0009_Chr7A_CHH.bw Z0009_Chr7B_CHH.bw Z0009_Chr7D_CHH.bw \
Z0009_all_chrom_CHH.bg

bedGraphToBigWig Z0001_all_chrom_CHH.bg chrom_sizes.txt Z0001_all_chrom_CHH.bw
bedGraphToBigWig Z0002_all_chrom_CHH.bg chrom_sizes.txt Z0002_all_chrom_CHH.bw
bedGraphToBigWig Z0003_all_chrom_CHH.bg chrom_sizes.txt Z0003_all_chrom_CHH.bw
bedGraphToBigWig Z0004_all_chrom_CHH.bg chrom_sizes.txt Z0004_all_chrom_CHH.bw
bedGraphToBigWig Z0005_all_chrom_CHH.bg chrom_sizes.txt Z0005_all_chrom_CHH.bw
bedGraphToBigWig Z0006_all_chrom_CHH.bg chrom_sizes.txt Z0006_all_chrom_CHH.bw
bedGraphToBigWig Z0007_all_chrom_CHH.bg chrom_sizes.txt Z0007_all_chrom_CHH.bw
bedGraphToBigWig Z0008_all_chrom_CHH.bg chrom_sizes.txt Z0008_all_chrom_CHH.bw
bedGraphToBigWig Z0009_all_chrom_CHH.bg chrom_sizes.txt Z0009_all_chrom_CHH.bw