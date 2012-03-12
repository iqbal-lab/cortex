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

  print STDERR "Usage: ./vcf_filter_by_info.pl [--invert] <file.vcf> " .
               "[<INFO_FIELD> <VALUE_REGEX>] ..\n";
  print STDERR "  Isaac Turner <isaac.turner\@dtc.ox.ac.uk> 2011/03/26\n";
  exit;
}

if(@ARGV < 1)
{
  print_usage();
}

my $invert = 0;

if($ARGV[0] =~ /^-?-invert$/i)
{
  $invert = 1;
  shift;
}

my $vcf_file = shift;

if((@ARGV % 2) != 0)
{
  print_usage();
}

my %searches = ();

for(my $i = 0; $i < @ARGV; $i+=2)
{
  $searches{$ARGV[$i]} = $ARGV[$i+1];
}

#
# Open VCF Handle
#
my $vcf_handle;

if(defined($vcf_file) && $vcf_file ne "-")
{
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

while(defined($vcf_entry = $vcf->read_entry()))
{
  $total_num_entries++;

  my $info_hashref = $vcf_entry->{'INFO'};
  my $flags_hashref = $vcf_entry->{'INFO_flags'};

  my $match = 1;

  my ($key,$search);

  for my $key (keys %searches)
  {
    if(!defined($flags_hashref->{$key}) &&
       (!defined($info_hashref->{$key}) ||
        $info_hashref->{$key} !~ /$searches{$key}/i))
    {
      $match = 0;
      last;
    }
  }

  if($match != $invert)
  {
    $num_of_filtered_entries++;
    $vcf->print_entry($vcf_entry);
  }
}

# Print stats
my $printed_percent = 100 * $num_of_filtered_entries / $total_num_entries;

print STDERR "vcf_filter_by_info.pl: " . num2str($num_of_filtered_entries) .
             " / " . num2str($total_num_entries) . " " .
             "(" . sprintf("%.2f", $printed_percent) . "%) variants printed\n";

close($vcf_handle);
