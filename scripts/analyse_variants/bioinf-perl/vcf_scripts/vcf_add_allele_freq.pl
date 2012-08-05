#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use VCFFile;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_add_allele_freq.pl [file.vcf]\n";
  print STDERR "  Add (or correct) allele frequency INFO tag in a VCF file\n";
  print STDERR "  Adds AF (Alt allele Frequency) INFO tag: 0..1\n";
  exit;
}

## Test for filtering
my $failed_vars_out = undef;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  open($failed_vars_out, ">-");
}
##

my $vcf_file = shift;

if(@ARGV > 0)
{
  print_usage("Too many arguments!");
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
elsif(-p STDIN)
{
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

my @sample_names = $vcf->get_list_of_sample_names();

if(@sample_names == 0)
{
  print_usage("VCF doesn't have any genotyped samples in it");
}

my $tag = "Allele Frequency, for each ALT allele, in the same order as listed";

$vcf->add_header_tag("INFO", "AF", 1, "Float", $tag);
$vcf->print_header();

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  # Works for multiple ALT frequencies separated by commas
  my @alleles = split(/,/, $vcf_entry->{'ALT'});

  my %allele_counts = ();
  my $ploidy;

  for my $sample (@sample_names)
  {
    if(!defined($vcf_entry->{$sample}->{'GT'}))
    {
      print_usage("vcf_add_allele_freq.pl: Missing genotype " .
                  "[sample:$sample; var:" . $vcf_entry->{'ID'} . " " .
                  "format:" . join(":", @{$vcf_entry->{'FORMAT'}}) . "]");
    }
  
    my @sample_alleles = split(/[\/\|]/, $vcf_entry->{$sample}->{'GT'});

    if(defined($ploidy) && @sample_alleles != $ploidy)
    {
      print_usage("vcf_add_allele_freq.pl: Ploidy mismatch " .
                  "[ploidy:$ploidy sample:$sample; " .
                  "var:" . $vcf_entry->{'ID'} . " " .
                  "GT:" . $vcf_entry->{$sample}->{'GT'} . "]");
    }

    $ploidy = scalar(@sample_alleles);

    for my $sample_allele (@sample_alleles)
    {
      $allele_counts{$sample_allele}++;
    }
  }

  my @alt_freqs = ();
  my $total_ploidy = $ploidy * @sample_names;

  for(my $i = 1; $i <= @alleles; $i++)
  {
    my $allele_count = defined($allele_counts{$i}) ? $allele_counts{$i} : 0;
    push(@alt_freqs, $allele_count / $total_ploidy);
  }

  $vcf_entry->{'INFO'}->{'AF'} = join(",", @alt_freqs);

  $vcf->print_entry($vcf_entry);
}

close($vcf_handle);
