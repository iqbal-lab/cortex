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
  
  print STDERR "Usage: ./vcf_filter_by_homopolymer.pl [--invert] <max_run> [in.vcf]\n";
  print STDERR "  Filter out homopolymer runs LONGER than <max_run> (min = 1)\n";
  print STDERR "  Variant is homopolymer run if one allele is explained by a\n";
  print STDERR "  run of > max_run bases.  \n";
  exit;
}

## Test for filtering
my $failed_vars_out = undef;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  open($failed_vars_out, ">-");
}
##

if(@ARGV < 1 || @ARGV > 3)
{
  print_usage();
}

my $invert = 0;

if(@ARGV > 1 && $ARGV[0] =~ /^-?-inv(ert)?$/i)
{
  $invert = 1;
  shift;
}

if(@ARGV > 2)
{
  print_usage();
}

my $max_run = shift;
my $vcf_file = shift; # optional

if($max_run !~ /^\d+$/ || $max_run < 1)
{
  print_usage("<max_run> must be an integer >= 1");
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

# Print non-PASS variants straight to stdout if -p passed
if(defined($failed_vars_out)) { $vcf->set_filter_failed($failed_vars_out);}

$vcf->print_header();

my $num_of_entries_printed = 0;
my $num_of_entries = 0;

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  $num_of_entries++;

  my $left_flank = $vcf_entry->{'INFO'}->{'left_flank'};
  my $right_flank = $vcf_entry->{'INFO'}->{'right_flank'};
  my $ref_allele = $vcf_entry->{'true_REF'};
  my $alt_allele = $vcf_entry->{'true_ALT'};

  if(!defined($left_flank) || !defined($right_flank))
  {
    print_usage("'left_flank' or 'right_flank' INFO field missing: " .
                "use vcf_add_flanks.pl");
  }

  my $homopolymer = is_homopolymer_run($left_flank, $ref_allele, $right_flank) ||
                    is_homopolymer_run($left_flank, $alt_allele, $right_flank);

  if($homopolymer != $invert)
  {
    $vcf->print_entry($vcf_entry);
    $num_of_entries_printed++;
  }
}

# Print stats
my $printed_percent = 100 * $num_of_entries_printed / $num_of_entries;

print STDERR "vcf_filter_by_homopolymer.pl: " .
             num2str($num_of_entries_printed) . " / " .
             num2str($num_of_entries) . " " .
             "(" . sprintf("%.2f", $printed_percent) . "%) variants printed\n";

close($vcf_handle);


#
# Methods
#

sub is_homopolymer_run
{
  my ($left_flank, $allele, $right_flank) = @_;

  if(length($allele) == 0)
  {
    return 0;
  }

  my $rev_right_flank = reverse($right_flank);
  my $rev_allele = reverse($allele);

  return (is_homopolymer_run_flank($left_flank, $allele) ||
          is_homopolymer_run_flank($rev_right_flank, $rev_allele));
}

sub is_homopolymer_run_flank
{
  my ($left_flank, $allele) = @_;

  my $char = substr($allele, 0, 1);
  
  my ($flank_match, $allele_match);

  if($left_flank =~ /((?:$char)+)$/i)
  {
    $flank_match = $1;
  }

  if($allele =~ /^((?:$char)+)$/i)
  {
    $allele_match = $1;
  }

  return (defined($flank_match) && defined($allele_match) &&
          length($flank_match) + length($allele_match) > $max_run);
}
