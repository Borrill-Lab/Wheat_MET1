#!/bin/bash
#SBATCH -p jic-medium
#SBATCH -t 00-12:00
#SBATCH -c 4
#SBATCH --mem=120000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=delfi.dorussen@jic.ac.uk
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/filter_by_chrom/slurm_output/filter-slurm-%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/filter_by_chrom/slurm_output/filter-slurm-%j.err

cd /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/filter_by_chrom/

sed -n '1,268261344p;268261345q' /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/methylation_extraction/Unknown_BU703-001Z0008_1.trimmed_bismark_bt2_pe.CX_report.txt > Unknown_BU703-001Z0008_1.Chr1A.CX_report.txt
sed -n '268261345,581167332p;581167333q' /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/methylation_extraction/Unknown_BU703-001Z0008_1.trimmed_bismark_bt2_pe.CX_report.txt > Unknown_BU703-001Z0008_1.Chr1B.CX_report.txt
sed -n '581167333,806387160p;806387161q' /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/methylation_extraction/Unknown_BU703-001Z0008_1.trimmed_bismark_bt2_pe.CX_report.txt > Unknown_BU703-001Z0008_1.Chr1D.CX_report.txt
sed -n '806387161,1159663434p;1159663435q' /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/methylation_extraction/Unknown_BU703-001Z0008_1.trimmed_bismark_bt2_pe.CX_report.txt > Unknown_BU703-001Z0008_1.Chr2A.CX_report.txt
sed -n '1159663435,1523453572p;1523453573q' /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/methylation_extraction/Unknown_BU703-001Z0008_1.trimmed_bismark_bt2_pe.CX_report.txt > Unknown_BU703-001Z0008_1.Chr2B.CX_report.txt
sed -n '1523453573,1820105918p;1820105919q' /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/methylation_extraction/Unknown_BU703-001Z0008_1.trimmed_bismark_bt2_pe.CX_report.txt > Unknown_BU703-001Z0008_1.Chr2D.CX_report.txt
sed -n '1820105919,2159234672p;2159234673q' /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/methylation_extraction/Unknown_BU703-001Z0008_1.trimmed_bismark_bt2_pe.CX_report.txt > Unknown_BU703-001Z0008_1.Chr3A.CX_report.txt
sed -n '2159234673,2536449719p;2536449720q' /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/methylation_extraction/Unknown_BU703-001Z0008_1.trimmed_bismark_bt2_pe.CX_report.txt > Unknown_BU703-001Z0008_1.Chr3B.CX_report.txt
sed -n '2536449720,2816384911p;2816384912q' /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/methylation_extraction/Unknown_BU703-001Z0008_1.trimmed_bismark_bt2_pe.CX_report.txt > Unknown_BU703-001Z0008_1.Chr3D.CX_report.txt
sed -n '2816384912,3152997001p;3152997002q' /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/methylation_extraction/Unknown_BU703-001Z0008_1.trimmed_bismark_bt2_pe.CX_report.txt > Unknown_BU703-001Z0008_1.Chr4A.CX_report.txt
sed -n '3152997002,3459787450p;3459787451q' /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/methylation_extraction/Unknown_BU703-001Z0008_1.trimmed_bismark_bt2_pe.CX_report.txt > Unknown_BU703-001Z0008_1.Chr4B.CX_report.txt
sed -n '3459787451,3693011655p;3693011656q' /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/methylation_extraction/Unknown_BU703-001Z0008_1.trimmed_bismark_bt2_pe.CX_report.txt > Unknown_BU703-001Z0008_1.Chr4D.CX_report.txt
sed -n '3693011656,4013418154p;4013418155q' /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/methylation_extraction/Unknown_BU703-001Z0008_1.trimmed_bismark_bt2_pe.CX_report.txt > Unknown_BU703-001Z0008_1.Chr5A.CX_report.txt
sed -n '4013418155,4337032235p;4337032236q' /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/methylation_extraction/Unknown_BU703-001Z0008_1.trimmed_bismark_bt2_pe.CX_report.txt > Unknown_BU703-001Z0008_1.Chr5B.CX_report.txt
sed -n '4337032236,4593735789p;4593735790q' /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/methylation_extraction/Unknown_BU703-001Z0008_1.trimmed_bismark_bt2_pe.CX_report.txt > Unknown_BU703-001Z0008_1.Chr5D.CX_report.txt
sed -n '4593735790,4873095882p;4873095883q' /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/methylation_extraction/Unknown_BU703-001Z0008_1.trimmed_bismark_bt2_pe.CX_report.txt > Unknown_BU703-001Z0008_1.Chr6A.CX_report.txt
sed -n '4873095883,5200864306p;5200864307q' /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/methylation_extraction/Unknown_BU703-001Z0008_1.trimmed_bismark_bt2_pe.CX_report.txt > Unknown_BU703-001Z0008_1.Chr6B.CX_report.txt
sed -n '5200864307,5416638604p;5416638605q' /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/methylation_extraction/Unknown_BU703-001Z0008_1.trimmed_bismark_bt2_pe.CX_report.txt > Unknown_BU703-001Z0008_1.Chr6D.CX_report.txt
sed -n '5416638605,5747932628p;5747932629q' /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/methylation_extraction/Unknown_BU703-001Z0008_1.trimmed_bismark_bt2_pe.CX_report.txt > Unknown_BU703-001Z0008_1.Chr7A.CX_report.txt
sed -n '5747932629,6087964118p;6087964119q' /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/methylation_extraction/Unknown_BU703-001Z0008_1.trimmed_bismark_bt2_pe.CX_report.txt > Unknown_BU703-001Z0008_1.Chr7B.CX_report.txt
sed -n '6087964119,6377460089p;6377460090q' /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/methylation_extraction/Unknown_BU703-001Z0008_1.trimmed_bismark_bt2_pe.CX_report.txt > Unknown_BU703-001Z0008_1.Chr7D.CX_report.txt