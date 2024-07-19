#!/usr/bin/perl -w

# script obtained from Philippa Borrill, modified by Delfi Dorussen for Bismark methylation extraction
#
# Aim of script is to run bismark_methylation_extractor on WGBS data mapped with bismark for multiple samples

#### paths and references:
my $index = "/jic/scratch/groups/Philippa-Borrill/References/iwgsc_ref_seq_1.1/BSseq-ref/";

#############################


my $read_path_triticum = "/jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/bismark_mapping/";
my $output_dir = "/jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/methylation_extraction/";

### lists of samples (text file containing the sample name and the .bam file, separated by tabs; Sample	Sample_folder/Sample_1.trimmed.bam):
my $input_for_bismark = "/jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/bismark_mapping/input_for_bismark.txt";


#############################

#open the input file and go through the lines one by one so go to each directories where the fastq.gz should be located
chdir("$read_path_triticum") or die "couldn't move to input directory";

open (INPUT_FILE, "$input_for_bismark") || die "couldn't open the input file $input_for_bismark!";
		    while (my $line = <INPUT_FILE>) {
			chomp $line;
my @array = split(/\t/,$line);
#print "\nmy line was: $line\n";

#print "\nmy array: @array\n";
#print "\narray element 1: @array[0]\n";

my $sample = $array[0];
my $fastq1 = $array[1];

chdir("$read_path_triticum") or die "couldn't move to specific read directory";


my $SLURM_header = <<"SLURM";
#!/bin/bash
#
# SLURM batch script to launch parallel tasks
#
#SBATCH -p jic-long
#SBATCH -t 28-00:00
#SBATCH -c 4
#SBATCH --mem=120000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=delfi.dorussen@jic.ac.uk
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/methylation_extraction/slurm_output/%x.%N.%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/methylation_extraction/slurm_output/%x.%N.%j.err


SLURM

 my $tmp_file = "$output_dir/tmp/bismark.$sample";


  open (SLURM, ">$tmp_file") or die "Couldn't open temp file\n";
  $SLURM_header = $SLURM_header;
  print SLURM "$SLURM_header\n\n";
  print SLURM "\ncd $read_path_triticum\n";


  print SLURM "set -e\n";

	print SLURM "source package 33c48798-0827-4add-8153-909c1bd83e89\n";
	print SLURM "bismark_methylation_extractor -p --bedGraph -o $output_dir --CX_context --cytosine_report --genome_folder $index $fastq1\n";

	close SLURM;
  system("sbatch $tmp_file");
 # unlink $tmp_file;

}

	    close(INPUT_FILE);
