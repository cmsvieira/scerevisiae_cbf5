#!/usr/bin/env perl
use strict;
use FileHandle;
use Bio::Adventure;
my $out = FileHandle->new(">>outputs/log.txt");
my $d = qx'date';
print $out "###Started /cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782/scripts/81parsesnp_scerevisiae.sh at ${d}";
chdir("/cbcb/nelsayed-scratch/atb/rnaseq/scerevisiae_cbf5_2017/preprocessing/v2/hpgl0782");
my $h = Bio::Adventure->new();


## Parse the SNP data and generate a modified scerevisiae genome.
##  This should read the file:
## outputs/vcfutils_scerevisiae/hpgl0782.bcf
##  and provide 4 new files:
## outputs/vcfutils_scerevisiae/hpgl0782_scerevisiae_count.txt
## outputs/vcfutils_scerevisiae/hpgl0782_scerevisiae_all.txt
## outputs/vcfutils_scerevisiae/hpgl0782_scerevisiae_pct.txt
## and a modified genome: outputs/vcfutils_scerevisiae/hpgl0782_scerevisiae_modified.fasta


use Bio::Adventure::SNP;
Bio::Adventure::SNP::Make_SNP_Ratio(
  $h,
  input => 'outputs/vcfutils_scerevisiae/hpgl0782.bcf',
  output => 'outputs/vcfutils_scerevisiae/hpgl0782_scerevisiae',
  species => 'scerevisiae',
  vcf_cutoff => '10',
  vcf_minpct => '0.8',
);

## The following lines give status codes and some logging
my $jobid = "";
$jobid = $ENV{SLURM_JOBID} if ($ENV{SLURM_JOBID});
my $end_d = qx'date';
print $out "####Finished ${jobid} 81parsesnp_scerevisiae.sh at ${end_d}.";
close($out);

