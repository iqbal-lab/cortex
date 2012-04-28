#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use VCFFile;
use UsefulModule; # for pretty_fraction

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_remove_dupes.pl [file.vcf]\n" .
"  Remove entries that match position, REF and ALT alleles. Assumes sorted VCF.\n";
  exit;
}

if(@ARGV > 1) {
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

$vcf->print_header();

my $num_of_entries = 0;
my $num_of_printed = 0;

my $prev_entry = $vcf->read_entry();

if(defined($prev_entry))
{
  $vcf->print_entry($prev_entry);
  $num_of_entries++;
  $num_of_printed++;
}

my $next_entry;

while(defined($next_entry = $vcf->read_entry()))
{
  $num_of_entries++;

  if($next_entry->{'CHROM'} ne $prev_entry->{'CHROM'} ||
     $next_entry->{'POS'} != $prev_entry->{'POS'} ||
     uc($next_entry->{'REF'}) ne uc($prev_entry->{'REF'}) || 
     uc($next_entry->{'ALT'}) ne uc($prev_entry->{'ALT'}))
  {
    $vcf->print_entry($next_entry);
    $num_of_printed++;
  }

  $prev_entry = $next_entry;
}

print STDERR "vcf_remove_dupes.pl: " .
      pretty_fraction($num_of_printed, $num_of_entries) . " variants printed\n";

close($vcf_handle);
