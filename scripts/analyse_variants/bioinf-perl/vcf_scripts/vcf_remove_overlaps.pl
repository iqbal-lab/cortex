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
"Usage: ./vcf_remove_overlaps.pl [OPTIONS] [file.vcf]
  Remove entries that overlap on the reference. Assumes sorted VCF with
  duplicates removed.

  OPTIONS:
  --padding <p>, -p <p>        Require at least <p> bp between variants [default: 0]
  --filter_txt <txt>, -f <txt> Add to / set the filter column instead of removing
  --invert, -i                 Invert selection (print overlapping variants only)

  Example:
    vcf-sort in.vcf | vcf_remove_dupes.pl | vcf_remove_overlaps.pl > out.vcf\n";

  exit;
}

## Test for filtering
my $failed_vars_out = undef;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  open($failed_vars_out, ">-");
}
##

my $padding = 0;
my $filter_txt;
my $invert = 0;

while(@ARGV >= 1)
{
  if($ARGV[0] =~ /^-?-p(adding)?$/i)
  {
    shift;
    $padding = shift;

    if($padding !~ /^\d+$/)
    {
      print_usage("--padding <p> requires a positive integer");
    }
  }
  elsif($ARGV[0] =~ /^-?-f(ilter_txt)?$/i)
  {
    shift;

    if(!defined($filter_txt = shift))
    {
      print_usage("--filter_txt <txt>  requires a FILTER tag");
    }
  }
  elsif($ARGV[0] =~ /^-?-i(nvert)?$/i)
  {
    shift;
    $invert = 1;
  }
  else
  {
    last;
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

# Print non-PASS variants straight to stdout if -p passed
if(defined($failed_vars_out)) { $vcf->set_filter_failed($failed_vars_out);}

$vcf->print_header();

my $num_of_entries = 0;
my $num_of_printed = 0;
my $num_overlaps = 0;

my ($prev_entry, $curr_entry);
my $var_end;
my $prev_chrom;
my $print_prev = 1;
my $print_curr;
my $prev_overlapping;

if(defined($prev_entry = $vcf->read_entry()))
{
  $num_of_entries++;
  $prev_chrom = $prev_entry->{'CHROM'};
  $var_end = $prev_entry->{'true_POS'} + length($prev_entry->{'REF'}) - 1;
}

while(defined($curr_entry = $vcf->read_entry()))
{
  $num_of_entries++;

  my $same_chrom = ($curr_entry->{'CHROM'} eq $prev_entry->{'CHROM'});

  $print_curr = (!$same_chrom || $var_end + $padding < $curr_entry->{'true_POS'} );

  $prev_overlapping = !($print_curr && $print_prev);

  if($prev_overlapping)
  {
    $num_overlaps++;
  }

  if($prev_overlapping == $invert)
  {
    $vcf->print_entry($prev_entry);
    $num_of_printed++;
  }
  elsif(defined($filter_txt))
  {
    vcf_add_filter_txt($prev_entry, $filter_txt);
    $vcf->print_entry($prev_entry);
    $num_of_printed++;
  }

  my $entry_end = $curr_entry->{'true_POS'} + length($curr_entry->{'REF'}) - 1;
  $var_end = $same_chrom ? max($var_end, $entry_end) : $entry_end;

  $prev_entry = $curr_entry;
  $print_prev = $print_curr;
}

if(!$print_prev)
{
  $num_overlaps++;
}

if($print_prev != $invert)
{
  $vcf->print_entry($prev_entry);
  $num_of_printed++;
}
elsif(defined($filter_txt))
{
  vcf_add_filter_txt($prev_entry, $filter_txt);
  $vcf->print_entry($prev_entry);
  $num_of_printed++;
}

print STDERR "vcf_remove_overlaps.pl: " .
             pretty_fraction($num_overlaps, $num_of_entries) .
             " variants overlapping\n";

print STDERR "vcf_remove_overlaps.pl: " .
             pretty_fraction($num_of_entries-$num_overlaps, $num_of_entries) .
             " variants non-overlapping\n";

print STDERR "vcf_remove_overlaps.pl: " .
             pretty_fraction($num_of_printed, $num_of_entries) .
             " variants printed\n";

close($vcf_handle);
