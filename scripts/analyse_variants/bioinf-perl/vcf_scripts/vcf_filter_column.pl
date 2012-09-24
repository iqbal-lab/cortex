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
  
  print STDERR "" .
"Usage: ./vcf_filter_column.pl [OPTIONS] [in.vcf]
 Prints variants that passed the filtering. If [in.vcf] is '-', reads from stdin.
  
 OPTIONS:
  --require <filter>  Require that a variant fits a certain filter.  Can use
                      mulitple times.  
  --ignore <filter>   Ignore a certain filter.  Multiple use allowed
  --ignore_others     Ignore all other filters
  --invert            Invert results

  --show_stats        Show stats for filter types\n";

  exit;
}

## Test for filtering
my $failed_vars_out = undef;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  open($failed_vars_out, ">-");
}
##

my %required = ();
my %ignored = ();
my $ignore_others = 0;
my $invert = 0;
my $show_stats = 0;
my $vcf_file;

while(@ARGV > 0)
{
  if($ARGV[0] =~ /^-?-invert?$/i)
  {
    shift;
    $invert = 1;
  }
  elsif($ARGV[0] =~ /^-?-show_stats?$/i)
  {
    shift;
    $show_stats = 1;
  }
  
  elsif($ARGV[0] =~ /^-?-ignore_others?$/i)
  {
    shift;
    $ignore_others = 1;
  }
  elsif($ARGV[0] =~ /^-?-require$/i)
  {
    shift;
    my $arg = shift;
    if(!defined($arg))
    {
      print_usage("--required <filter> needs a filter name");
    }
    $required{uc($arg)} = 1;
  }
  elsif($ARGV[0] =~ /^-?-ignore$/i)
  {
    shift;
    my $arg = shift;
    if(!defined($arg))
    {
      print_usage("--ignore <filter> needs a filter name");
    }
    $ignored{$arg} = 1;
  }
  elsif($ARGV[0] =~ /^-/)
  {
    print_usage("Unknown argument: '$ARGV[0]'");
  }
  else
  {
    last;
  }
}

$vcf_file = shift;

if(@ARGV > 0)
{
  print_usage("Excess arguments");
}

my $num_of_required_filters = scalar(keys %required);
my $num_of_ignored_filters = scalar(keys %ignored);

if($ignore_others && $num_of_ignored_filters)
{
  print STDERR "vcf_filter_column.pl: Note --ignore_others makes other " .
               "--ignore arguments redundant";
}

if($num_of_required_filters == 0)
{
  $ignored{'PASS'} = 1;
}

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

# Print non-PASS variants straight to stdout if -p passed
if(defined($failed_vars_out)) { $vcf->set_filter_failed($failed_vars_out);}

my %columns = $vcf->get_columns_hash();
if(!defined($columns{'FILTER'}))
{
  print_usage("VCF file is missing 'FILTER' column");
}

$vcf->print_header();

my %filter_counts = ();
my %filter_combos_counts = ();

my $num_of_filtered_entries = 0;
my $num_of_vars = 0;

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry))
{
  $num_of_vars++;

  # Convert to uppercase
  my @filters = map {uc($_)} split(",", $vcf_entry->{'FILTER'});

  if($show_stats)
  {
    for my $filter (@filters)
    {
      $filter_counts{$filter}++;
    }

    my $combo = join(",", sort @filters);
    $filter_combos_counts{$combo}++;
  }

  my $print = print_var(@filters);

  if($print)
  {
    $num_of_filtered_entries++;
    $vcf->print_entry($vcf_entry);
  }
}

print STDERR "vcf_filter_column.pl: " .
             pretty_fraction($num_of_filtered_entries, $num_of_vars) .
             " variants printed\n";

if($show_stats)
{
  for my $filter (sort keys %filter_counts)
  {
    print STDERR "vcf_filter_column.pl: " .
                 pretty_fraction($filter_counts{$filter}, $num_of_vars) .
                 " caught by filter '$filter'\n";
  }

  print STDERR "vcf_filter_column.pl: Filter combinations:\n";

  for my $filter (sort keys %filter_combos_counts)
  {
    print STDERR "vcf_filter_column.pl: " .
                 pretty_fraction($filter_combos_counts{$filter}, $num_of_vars) .
                 " '$filter' " .
                 "[".(print_var(split(",", $filter)) ? "" : "not ")."printed]\n";
  }
}

close($vcf_handle);

sub print_var
{
  my @filters = @_;

  # must have the same number of required filters as were specified on cmdline
  my $num_of_required = 0;
  my $fail = 0;

  for my $filter (@filters)
  {
    if($required{$filter})
    {
      $num_of_required++;
    }
    elsif(!$ignore_others && !defined($ignored{$filter}))
    {
      $fail = 1;
      last;
    }
  }

  return ($num_of_required == $num_of_required_filters && !$fail) != $invert;
}
