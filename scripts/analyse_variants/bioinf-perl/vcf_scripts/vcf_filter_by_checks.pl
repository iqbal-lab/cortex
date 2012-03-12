#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use VCFFile;
use UsefulModule; # num2str

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }
  
  print STDERR "Usage: ./vcf_filter_by_checks.pl [--invert] [in.vcf]\n";
  print STDERR "  Prints variants that passed the filtering\n";
  print STDERR "  If [in.vcf] is '-', reads from stdin\n";
  print STDERR "  --invert  Print variants that failed\n";
  exit;
}

if(@ARGV > 2) {
  print_usage();
}

my $invert = 0;
my $vcf_file;

if(@ARGV == 2)
{
  if($ARGV[0] =~ /^-?-invert?$/i)
  {
    shift;
    $invert = 1;
  }
  else
  {
    print_usage();
  }
}

$vcf_file = shift;

#
# Open VCF Handle
#
my $vcf_handle;

if(defined($vcf_file) && $vcf_file ne "-") {
  open($vcf_handle, $vcf_file) or die("Cannot open VCF file '$vcf_file'\n");
}
elsif(-p STDIN) {
  # STDIN is connected to a pipe
  open($vcf_handle, "<&=STDIN") or die("Cannot read pipe");
}
else
{
  print_usage("Must specify or pipe in a VCF file");
}

#
# Read VCF
#
my $vcf = new VCFFile($vcf_handle);

$vcf->print_header();

my $num_of_filtered_entries = 0;
my $total_num_entries = 0;

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry))
{
  $total_num_entries++;

  if(($vcf_entry->{'FILTER'} =~ /^PASS$/i) != $invert)
  {
    $num_of_filtered_entries++;
    $vcf->print_entry($vcf_entry);
  }
}

# Print filtered rate
my $printed_percent = 100 * $num_of_filtered_entries / $total_num_entries;

print STDERR "vcf_filter_by_checks.pl: " . num2str($num_of_filtered_entries) .
             " / " . num2str($total_num_entries) . " " .
             "(" . sprintf("%.2f", $printed_percent) . "%) variants printed\n";

close($vcf_handle);
