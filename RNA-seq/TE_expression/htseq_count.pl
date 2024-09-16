#!/usr/bin/perl -w

# Delfi Dorussen
#
# Aim of script is to run htseq for multiple samples of mapped met1 RNA-seq data

#### paths and references:
my $read_path_triticum = "/jic/scratch/groups/Philippa-Borrill/Delfi/MET1_TE/hisat_mapping/";
my $output_dir = "/jic/scratch/groups/Philippa-Borrill/Delfi/MET1_TE/htseq_count/";

### lists of samples (text file containing directory/subdirectory with .fastq to map e.g. each line should look like: ERP004505/ERR392073/ in these subdirectories are the fastq.gz - text file must be in $output_dir):
my $input_list_dir = "/jic/scratch/groups/Philippa-Borrill/Delfi/MET1_RNA/trim/";
my $input_for_trim = "/jic/scratch/groups/Philippa-Borrill/Delfi/MET1_RNA/trim/input_for_kallisto.txt";

#############################

#open the input file and go through the lines one by one so go to each directories where the fastq.gz should be located
chdir("$input_list_dir") or die "couldn't move to input directory";

open (INPUT_FILE, "$input_for_trim") || die "couldn't open the input file $input_for_trim!";
		    while (my $line = <INPUT_FILE>) {
			chomp $line;
my @array = split(/\t/,$line);
#print "\nmy line was: $line\n";

#print "\nmy array: @array\n";
#print "\narray element 1: @array[0]\n";

my $sample = $array[0];
my $f1 = $array[1];
my $f2 = $array[2];


chdir("$read_path_triticum") or die "couldn't move to specific read directory $read_path_triticum";


my $SLURM_header = <<"SLURM";
#!/bin/bash
#
# SLURM batch script to launch parallel htseq-count tasks
#
#SBATCH -p jic-medium
#SBATCH -t 1-00:00
#SBATCH -c 4
#SBATCH --mem=200000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=delfi.dorussen@jic.ac.uk
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_TE/htseq_count/slurm_output/htseq-$sample-slurm-%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_TE/htseq_count/slurm_output/htseq-$sample-slurm-%j.err
SLURM

my $tmp_file = "$output_dir/tmp/htseq.$sample";


  open (SLURM, ">$tmp_file") or die "Couldn't open temp file\n";
  $SLURM_header = $SLURM_header;
  print SLURM "$SLURM_header\n\n";
  print SLURM "\ncd $read_path_triticum\n";

  print SLURM "set -e\n";

  print SLURM "source package c8476f5b-da1a-45a2-a27c-c0bb32763998\n";

	### this part all needs editing!!!! ###
	print SLURM "htseq-count -f bam -m union --nonunique none -t repeat_region -i ID -o $output_dir/$sample.sam $sample.sorted.bam /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/DMR_overlap/iwgsc_refseqv1.0_TransposableElements_2017Mar13_complete.gff3\n";

	close SLURM;
  	system("sbatch $tmp_file");
 # unlink $tmp_file;

}

            close(INPUT_FILE);






