#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use VCFFile;
use UsefulModule; # for num2str & pretty_fraction

sub print_usage
{
  for my $err (@_)
  {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./vcf_filter_variants.pl [--invert] <INDEL|CLEAN_INDEL|SNP|MNP|NP> [in.vcf]
  Prints variants of a particular type
  --invert      Everything but option

  Options:
   INDEL        SVLEN != 0
   CLEAN_INDEL  (SVLEN != 0) & (one allele 0bp)
   SNP          both alleles 1bp
   MNP          (SVLEN == 0) & (allele length > 1)
   NP           SVLEN == 0
   COMPLEX      !snp && (length(allele1) > 0) & (length(allele2) > 0)\n";

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

my $filter_invert = 0;

if(@ARGV > 1 && $ARGV[0] =~ /^-?-c(lean)?$/i)
{
  $filter_invert = 1;
  shift;
}

my $filter = shift;
$filter = uc($filter);

my @filters = qw(INDEL CLEAN_INDEL SNP MNP NP COMPLEX);

my @plural = map {$_."S"} @filters;

if(grep {$filter eq $_} @plural)
{
  # Chop 's' off
  $filter = substr($filter, 0, -1);
}
elsif(!grep {$_ eq $filter} @filters)
{
  print_usage("Unknown option: '$filter' " .
              "(expected one of: ".join(",",@filters).")");
}

my $vcf_file = shift;

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

$vcf->print_header();

my $num_of_variants = 0;
my $num_of_printed = 0;

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  $num_of_variants++;

  my $print;

  if($filter eq "INDEL")
  {
    $print = ($vcf_entry->{'INFO'}->{'SVLEN'} != 0);
  }
  elsif($filter eq "CLEAN_INDEL")
  {
    $print = ($vcf_entry->{'INFO'}->{'SVLEN'} != 0 &&
              (length($vcf_entry->{'true_REF'}) == 0 ||
               length($vcf_entry->{'true_ALT'}) == 0));
  }
  elsif($filter eq "SNP")
  {
    $print = (length($vcf_entry->{'true_REF'}) == 1 &&
              length($vcf_entry->{'true_ALT'}) == 1);
  }
  elsif($filter eq "MNP")
  {
    $print = ($vcf_entry->{'INFO'}->{'SVLEN'} == 0 &&
              length($vcf_entry->{'true_REF'}) > 1);
  }
  elsif($filter eq "NP")
  {
    $print = ($vcf_entry->{'INFO'}->{'SVLEN'} == 0);
  }
  elsif($filter eq "COMPLEX")
  {
    my $ref_len = length($vcf_entry->{'true_REF'});
    my $alt_len = length($vcf_entry->{'true_ALT'});

    $print = ($ref_len > 0 && $alt_len > 0 && ($ref_len != 1 || $alt_len != 1));
  }

  if($print != $filter_invert)
  {
    $vcf->print_entry($vcf_entry);
    $num_of_printed++;
  }
}

# Print filtered rate
print STDERR "vcf_filter_variants.pl: " .
             "(" . ($filter_invert ? "not " : "") . $filter . ") " .
             pretty_fraction($num_of_printed, $num_of_variants) . " " .
             "variants printed\n";

close($vcf_handle);
