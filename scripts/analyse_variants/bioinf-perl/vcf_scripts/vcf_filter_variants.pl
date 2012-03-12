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
    chomp($err);
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./vcf_filter_variants.pl [--invert] <INDEL|CLEAN_INDEL|SNP|MNP|NP> " .
  "[in.vcf]\n" .
"  Prints variants were SVLEN != 0\n" .
"  --invert      Everything but option\n" .
"  \n" .
"  Options:\n" .
"   INDEL        SVLEN != 0\n" .
"   CLEAN_INDEL  (SVLEN != 0) & (one allele 0bp)\n" .
"   SNP          both alleles 1bp\n" .
"   MNP          (SVLEN == 0) & (allele length > 1)\n" .
"   NP           SVLEN == 0\n";

  exit;
}

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

my @filters = qw(INDEL INDELS CLEAN_INDEL CLEAN_INDELS SNP SNPS MNP MNPS NP NPS);

if(!grep {$_ eq $filter} @filters)
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

$vcf->print_header();

my $num_of_variants = 0;
my $num_of_printed = 0;

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  $num_of_variants++;

  my $print;

  if($filter eq "INDEL" || $filter eq "INDELS")
  {
    $print = ($vcf_entry->{'INFO'}->{'SVLEN'} != 0);
  }
  elsif($filter eq "CLEAN_INDEL" || $filter eq "CLEAN_INDELS")
  {
    $print = ($vcf_entry->{'INFO'}->{'SVLEN'} != 0 &&
              (length($vcf_entry->{'true_REF'}) == 0 ||
               length($vcf_entry->{'true_ALT'}) == 0));
  }
  elsif($filter eq "SNP" || $filter eq "SNPS")
  {
    $print = (length($vcf_entry->{'true_REF'}) == 1 &&
              length($vcf_entry->{'true_ALT'}) == 1);
  }
  elsif($filter eq "MNP" || $filter eq "MNPS")
  {
    $print = ($vcf_entry->{'INFO'}->{'SVLEN'} == 0 &&
              length($vcf_entry->{'true_REF'}) > 1);
  }
  elsif($filter eq "NP" || $filter eq "NPS")
  {
    $print = ($vcf_entry->{'INFO'}->{'SVLEN'} == 0);
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
