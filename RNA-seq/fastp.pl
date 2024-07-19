#!/usr/bin/perl -w

# Delfi Dorussen, adapted from Philippa Borrill
#
# Aim of script is to run fastp for mutliple samples

#### paths and references:
my $read_path_triticum = "/jic/scratch/groups/Philippa-Borrill/Delfi/MET1_RNA/";
my $output_dir = "/jic/scratch/groups/Philippa-Borrill/Delfi/MET1_RNA/trim/";

### lists of samples (text file containing the sample name and each of the two read files, separate by tabs; Sample	Sample_1.fq.gz	Sample_2.fq.gz):
my $input_list_dir = "/jic/scratch/groups/Philippa-Borrill/Delfi/MET1_RNA/";
my $input_for_trim = "/jic/scratch/groups/Philippa-Borrill/Delfi/MET1_RNA/list_of_samples.txt";

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
# SLURM batch script to launch parallel tasks
#
#SBATCH -p jic-medium
#SBATCH -t 0-04:00
#SBATCH -c 4
#SBATCH --mem=120000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=delfi.dorussen@jic.ac.uk
#SBATCH -J trim
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_RNA/trim/slurm_output/%x.%N.%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_RNA/trim/slurm_output/%x.%N.%j.err
SLURM

 my $tmp_file = "$output_dir/tmp/trim.$sample";


  open (SLURM, ">$tmp_file") or die "Couldn't open temp file\n";
  $SLURM_header = $SLURM_header;
  print SLURM "$SLURM_header\n\n";
  print SLURM "\ncd $read_path_triticum\n";

  print SLURM "set -e\n";

  print SLURM "source package fastp-0.20.0\n";

	### this part all needs editing!!!! ###
	print SLURM "fastp --in1 $f1 --in2 $f2 --out1 $output_dir/$sample"."_1.trimmed.fq.gz --out2 $output_dir/$sample"."_2.trimmed.fq.gz -l 80 -g --detect_adapter_for_pe\n";

    close SLURM;
  	system("sbatch $tmp_file");
 # unlink $tmp_file;

}

            close(INPUT_FILE);

