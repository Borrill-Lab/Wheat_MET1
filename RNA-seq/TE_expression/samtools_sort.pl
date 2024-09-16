#!/usr/bin/perl -w

# Delfi Dorussen
#
# Aim of script is to run samtools for multiple samples for met1 RNA-seq data

#### paths and references:
my $read_path_triticum = "/jic/scratch/groups/Philippa-Borrill/Delfi/MET1_RNA/trim/";
my $output_dir = "/jic/scratch/groups/Philippa-Borrill/Delfi/MET1_TE/hisat_mapping/";

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
# SLURM batch script to launch parallel samtools tasks
#
#SBATCH -p jic-medium
#SBATCH -t 1-00:00
#SBATCH -c 4
#SBATCH --mem=200000
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=delfi.dorussen@jic.ac.uk
#SBATCH -o /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_TE/hisat_mapping/slurm_output/sort-slurm-%j.out
#SBATCH -e /jic/scratch/groups/Philippa-Borrill/Delfi/MET1_TE/hisat_mapping/slurm_output/sort-slurm-%j.err
SLURM

my $tmp_file = "$output_dir/tmp/trim.$sample";


  open (SLURM, ">$tmp_file") or die "Couldn't open temp file\n";
  $SLURM_header = $SLURM_header;
  print SLURM "$SLURM_header\n\n";
  print SLURM "\ncd $output_dir\n";

  print SLURM "set -e\n";

  print SLURM "source package c92263ec-95e5-43eb-a527-8f1496d56f1a\n";

	### this part all needs editing!!!! ###
	print SLURM "samtools sort -@ 4 -o $sample.sorted.bam $sample.sam\n";

	close SLURM;
  	system("sbatch $tmp_file");
 # unlink $tmp_file;

}

            close(INPUT_FILE);





