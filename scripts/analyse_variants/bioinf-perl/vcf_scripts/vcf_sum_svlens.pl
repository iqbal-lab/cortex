#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use VCFFile;
use UsefulModule;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_sum_svlens.pl [file.vcf]\n";
  print STDERR "  Sum polarised SVLEN (uses tag AALEN).\n";
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
if($skip_failed_vars) { $vcf->set_filter_failed(undef); }

my $vcf_entry;

my $num_of_variants = 0;
my $num_of_variants_polarised = 0;

my $num_mnps = 0;
my $sum_mnps = 0;

my $num_insertions = 0;
my $sum_insertions = 0;
my $num_deletions = 0;
my $sum_deletions = 0;

while(defined($vcf_entry = $vcf->read_entry()))
{
  $num_of_variants++;
  my $aalen = $vcf_entry->{'INFO'}->{'AALEN'};
  
  if(defined($aalen) && $aalen ne ".")
  {
    $num_of_variants_polarised++;

    if($aalen > 0)
    {
      $num_insertions++;
      $sum_insertions += $aalen;
    }
    elsif($aalen < 0)
    {
      $num_deletions++;
      $sum_deletions += $aalen;
    }
    else
    {
      $num_mnps++;
      $sum_mnps += length($vcf_entry->{'true_REF'});
    }
  }
}

my $num_of_indels = $num_insertions + $num_deletions;

print "".num2str($num_of_variants_polarised)." / ".num2str($num_of_variants)." ".
      sprintf("%.2f", (100*$num_of_variants_polarised / $num_of_variants)) . "% " .
      " variants polarised\n";

print "Mean Ins/Del per indel: " .
      sprintf("%.3f", ($sum_insertions+$sum_deletions) / $num_of_indels)."bp\n";

if($num_insertions > 0)
{
  print "Total inserted bases = ".num2str($sum_insertions) ." bp in " .
        num2str($num_insertions) . " insertions (" . 
        sprintf("%.2f", $sum_insertions/$num_insertions) . "bp each)\n";
}

if($num_deletions > 0)
{
  print "Total deleted bases = ".num2str($sum_deletions) ." bp in " .
        num2str($num_deletions) . " deletions (" . 
        sprintf("%.2f", $sum_deletions/$num_deletions) . "bp each)\n";
}

if($num_mnps > 0)
{
  print "Total bases in MNPs = ".num2str($sum_mnps) ." bp in " .
        num2str($num_mnps) . " variants (" . 
        sprintf("%.2f", $sum_mnps/$num_mnps) . "bp each)\n";
}

close($vcf_handle);
