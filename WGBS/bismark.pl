#!/usr/bin/perl -w

# script obtained from Philippa Borrill, modified by Delfi Dorussen for Bismark mapping
#
# Aim of script is to run bismark on WGBS data for multiple samples to a common reference

#### paths and references:
my $index = '/jic/scratch/groups/Philippa-Borrill/References/iwgsc_ref_seq_1.1/BSseq-ref';

#############################


my $read_path_triticum = "/jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/trim/";
my $output_dir = "/jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/bismark_mapping/";

### lists of samples (text file containing the sample name and each of the two trimmed read files, separate by tabs; Sample	Sample_1.trimmed.fq.gz	Sample_2.trimmed.fq.gz):
my $input_for_bismark = "/jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/trim/input_for_bismark.txt";

## Sample Sample_1.fq.gz Sample_2.fq.gz

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
my $fastq2 = $array[2];


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
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/bismark_mapping/slurm_output/%x.%N.%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_WGBS/bismark_mapping/slurm_output/%x.%N.%j.err


SLURM

 my $tmp_file = "$output_dir/tmp/bismark.$sample";


  open (SLURM, ">$tmp_file") or die "Couldn't open temp file\n";
  $SLURM_header = $SLURM_header;
  print SLURM "$SLURM_header\n\n";
  print SLURM "\ncd $read_path_triticum\n";


  print SLURM "set -e\n";

	print SLURM "source package 33c48798-0827-4add-8153-909c1bd83e89\n";
	print SLURM "bismark $index -1 $fastq1 -2 $fastq2 -o $output_dir/$sample/ -p 4\n";

	close SLURM;
  system("sbatch $tmp_file");
 # unlink $tmp_file;

}

	    close(INPUT_FILE);
