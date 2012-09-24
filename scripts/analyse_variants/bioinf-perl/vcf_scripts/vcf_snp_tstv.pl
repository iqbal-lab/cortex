#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use VCFFile;
use UsefulModule;

## Config
my $csvsep = ",";
#

sub print_usage
{
  for my $err (@_)
  {
    print STDERR "Error: $err\n";
  }
  
  print STDERR
"Usage: ./vcf_snp_tstv.pl [in.vcf]
  Transition (Ts): A <-> G; Transversion (Tv) is all other substitutions
  Although there are twice as many Tv as Tv combinations, Ts make up
  approximately 2/3rds of all SNPs[1]
  [1] Collins DW, Jukes TH (April 1994). 'Rates of transition and
  transversion in coding sequences since the human-rodent divergence'\n";

  exit;
}

## Test for filtering
my $skip_failed_vars = 0;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  $skip_failed_vars = 1;
}
##

if(@ARGV > 1)
{
  print_usage();
}

my $vcf_file = shift;

#
# Open VCF Handle
#
my $vcf_handle;

if(defined($vcf_file) && $vcf_file ne "-")
{
  open($vcf_handle, $vcf_file)
    or print_usage("Cannot open VCF file '$vcf_file'\n");
}
elsif(-p STDIN) {
  # STDIN is connected to a pipe
  open($vcf_handle, "<&=STDIN") or print_usage("Cannot read pipe");
}
else
{
  print_usage("Must specify or pipe in a VCF file");
}

#
# Read VCF
#
my $vcf = new VCFFile($vcf_handle);

# Skip non-PASS variants if -p passed
if($skip_failed_vars) { $vcf->set_filter_failed(undef);}

my $num_of_transitions = 0;
my $num_of_transversions = 0;

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  my $ref = uc($vcf_entry->{'true_REF'});
  my $alt = uc($vcf_entry->{'true_ALT'});

  if(length($ref) == 1 && length($alt) == 1)
  {
    my $snp = $ref.$alt;

    if($snp eq "AG" || $snp eq "GA" || $snp eq "CT" || $snp eq "TC")
    {
      $num_of_transitions++;
    }
    else
    {
      $num_of_transversions++;
    }
  }
}

close($vcf_handle);

print "Transition (Ts) A<->G;  Transversion: all other SNPs\n";
print "Twice as many Tv values than Ts, but 2/3rds of SNPs are Ts\n";
print "Ts: ".num2str($num_of_transitions)."; " .
      "Tv: ".num2str($num_of_transversions)."\n";

if($num_of_transversions > 0)
{
  print sprintf("%.3f", $num_of_transitions / $num_of_transversions) .
        " Ts / Tv\n";
}

print "Total SNPs: ".num2str($num_of_transitions+$num_of_transversions)."\n";

