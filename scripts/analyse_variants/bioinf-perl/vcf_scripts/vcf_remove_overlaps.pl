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

  print STDERR "" .
"Usage: ./vcf_remove_overlaps.pl [--padding <p>] [file.vcf]
  Remove entries that match position, REF and ALT alleles. Assumes sorted VCF
  with duplicates removed. Example:
    vcf-sort in.vcf | vcf_remove_dupes.pl | vcf_remove_overlaps.pl > out.vcf
    
  --padding <p>  require at least <p> bp between variants [default: 0]\n";
  exit;
}

my $padding = 0;

if($ARGV[0] =~ /^-?-padding$/i)
{
  shift;
  $padding = shift;
  
  if($padding !~ /^\d+$/)
  {
    print_usage("--padding <p> requires a positive integer");
  }
}

my $vcf_file = shift;

if(@ARGV > 0)
{
  # Too many arguments passed
  print_usage();
}

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

my ($prev_entry, $next_entry);
my $print_prev = 0;

if(defined($prev_entry = $vcf->read_entry()))
{
  $num_of_entries++;
  $print_prev = 1;
}

while(defined($next_entry = $vcf->read_entry()))
{
  $num_of_entries++;

  if($prev_entry->{'CHROM'} ne $next_entry->{'CHROM'} ||
     $prev_entry->{'true_POS'} + length($prev_entry->{'true_REF'}) + $padding
      < $next_entry->{'true_POS'})
  {
    if($print_prev)
    {
      $vcf->print_entry($prev_entry);
      $num_of_printed++;
    }

    $print_prev = 1;
  }
  else
  {
    $print_prev = 0;
  }

  $prev_entry = $next_entry;
}

if($print_prev)
{
  $vcf->print_entry($prev_entry);
  $num_of_printed++;
}

print STDERR "vcf_remove_overlaps.pl: " .
      pretty_fraction($num_of_printed, $num_of_entries) . " variants printed\n";

close($vcf_handle);
