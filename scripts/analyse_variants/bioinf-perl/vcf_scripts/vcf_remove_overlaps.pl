#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(min max sum);

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
  Remove entries that overlap on the reference. Assumes sorted VCF with
  duplicates removed.

  --padding <p>  require at least <p> bp between variants [default: 0]

  Example:
    vcf-sort in.vcf | vcf_remove_dupes.pl | vcf_remove_overlaps.pl > out.vcf\n";

  exit;
}

my $padding = 0;

if(@ARGV > 0 && $ARGV[0] =~ /^-?-padding$/i)
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

my ($prev_entry, $entry);
my $var_end;
my $prev_chrom;
my $print_prev = 1;

if(defined($prev_entry = $vcf->read_entry()))
{
  $num_of_entries++;
  $prev_chrom = $prev_entry->{'CHROM'};
  $var_end = $prev_entry->{'true_POS'} + length($prev_entry->{'REF'}) - 1;
}

while(defined($entry = $vcf->read_entry()))
{
  $num_of_entries++;

  my $same_chrom = ($entry->{'CHROM'} eq $prev_entry->{'CHROM'});

  if(!$same_chrom || $var_end + $padding < $entry->{'true_POS'})
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

  my $entry_end = $entry->{'true_POS'} + length($entry->{'REF'}) - 1;
  $var_end = $same_chrom ? max($var_end, $entry_end) : $entry_end;
  $prev_entry = $entry;
}

if($print_prev)
{
  $vcf->print_entry($prev_entry);
  $num_of_printed++;
}

print STDERR "vcf_remove_overlaps.pl: " .
      pretty_fraction($num_of_printed, $num_of_entries) . " variants printed\n";

close($vcf_handle);
