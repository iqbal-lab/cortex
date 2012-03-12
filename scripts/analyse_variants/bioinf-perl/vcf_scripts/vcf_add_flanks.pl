#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use GeneticsModule;
use UsefulModule; # num2str
use VCFFile;
use FASTNFile;

use List::Util qw(min max);

#   Isaac Turner <isaac.turner\@dtc.ox.ac.uk> 2011/09/21

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_add_flanks.pl [FILTERS] <flank_size> <in.vcf> " .
               "<ref_name> <ref1.fa ..>\n";
  print STDERR "  Adds the flanking sequence to VCF files\n";
  print STDERR "  --filter_ref_match: remove variants where ref allele " .
               "doesn't match\n";
  print STDERR "  --filter_Ns: remove variants where flank+[ref/alt]+flank " .
               "contains Ns or flanks aren't long enough\n";
  print STDERR "  - <ref_name> is used in the VCF header\n";
  print STDERR "  - If <in.vcf> is '-', reads from stdin\n";
  exit;
}

if(@ARGV < 4)
{
  print_usage();
}

my $filter_by_ref_match = 0;
my $filter_Ns = 0;

for my $i (1,2)
{
  if($ARGV[0] =~ /^-?-filter_ref(?:_match)?$/i)
  {
    shift(@ARGV);
    $filter_by_ref_match = 1;
  }
  elsif($ARGV[0] =~ /^-?-filter_ns$/i)
  {
    shift(@ARGV);
    $filter_Ns = 1;
  }
}

my $flank_size = shift;
my $vcf_file = shift;
my $ref_name = shift;
my @ref_files = @ARGV;

if($flank_size !~ /^\d+$/ || $flank_size < 0)
{
  print_usage("Flank size must be a +ve integers ('$flank_size' invalid)");
}

#
# Open VCF Handle
#
my $vcf_handle;

if($vcf_file ne "-")
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
# Check reference genome files exist
#
for my $ref_file (@ref_files)
{
  if(!(-e $ref_file)) {
    print_usage("Reference file '$ref_file' doesn't exist");
  }
  elsif(!(-r $ref_file)) {
    print_usage("Reference file '$ref_file' isn't readable");
  }
}

#
# Load reference
#
my ($ref_genomes_hashref, $qual) = read_all_from_files(@ref_files);

my %ref_genomes = %$ref_genomes_hashref;

# Correct '1..' to 'chr1...' etc
while(my ($key,$value) = each(%ref_genomes))
{
  my $new_key = get_clean_chr_name($key);

  if($new_key ne $key)
  {
    $ref_genomes{$new_key} = uc($ref_genomes{$key});
    delete($ref_genomes{$key});
  }
  else {
    $ref_genomes{$key} = uc($ref_genomes{$key});
  }
}

#
# Read VCF
#
my $vcf = new VCFFile($vcf_handle);

# Add tags to header and print
my $tag_description = "$flank_size bp adjacent in ref genome '$ref_name'";
$vcf->add_header_tag("INFO", "left_flank", 1, "String", $tag_description);
$vcf->add_header_tag("INFO", "right_flank", 1, "String", $tag_description);
$vcf->print_header();

my $num_ref_mismatch = 0;
my $num_ns_in_alleles_or_flanks = 0;
my $num_printed = 0;
my $total_num_entries = 0;

my $vcf_entry;

my %missing_chrs = ();

while(defined($vcf_entry = $vcf->read_entry()))
{
  $total_num_entries++;
  my $chr = $vcf_entry->{'CHROM'};
  
  if(!defined($ref_genomes{$chr}))
  {
    $missing_chrs{$chr} = 1;
    next;
  }
  
  # Get coordinates, convert to zero based
  my $var_start = $vcf_entry->{'true_POS'}-1;
  my $ref_allele = $vcf_entry->{'true_REF'};
  my $alt_allele = $vcf_entry->{'true_ALT'};
  
  my $var_length = length($ref_allele);

  my $left_flank_start = max(0, $var_start-$flank_size);
  my $left_flank_length = $var_start - $left_flank_start;

  my $right_flank_start = $var_start + $var_length;
  my $right_flank_length = min($flank_size,
                               length($ref_genomes{$chr}) - $right_flank_start);

  my $vcf_info = $vcf_entry->{'INFO'};

  if($left_flank_start < 0 ||
     $left_flank_start + $left_flank_length > length($ref_genomes{$chr}))
  {
    print STDERR "Out of bounds: $left_flank_start $left_flank_length " .
                 "($chr ".length($ref_genomes{$chr}).")\n";
    
    open(my $stderr, ">&", STDERR) or die("Cannot open STDERR");
    $vcf->print_entry($vcf_entry, $stderr);
    #die();
  }
  elsif($right_flank_start < 0 ||
        $right_flank_start + $right_flank_length > length($ref_genomes{$chr}))
  {
    print STDERR "Out of bounds: $right_flank_start $right_flank_length " .
                 "($chr ".length($ref_genomes{$chr}).")\n";
    
    open(my $stderr, ">&", STDERR) or die("Cannot open STDERR");
    $vcf->print_entry($vcf_entry, $stderr);
    #die();
  }

  $vcf_info->{'left_flank'} = substr($ref_genomes{$chr},
                                     $left_flank_start, $left_flank_length);

  $vcf_info->{'right_flank'} = substr($ref_genomes{$chr},
                                      $right_flank_start, $right_flank_length);
  
  my $ref_seq = substr($ref_genomes{$chr}, $var_start, $var_length);
  
  if($filter_by_ref_match && $ref_seq ne uc($ref_allele))
  {
    $num_ref_mismatch++;
  }
  elsif($filter_Ns &&
        ($vcf_entry->{'REF'} =~ /n/i ||
         $vcf_entry->{'ALT'} =~ /n/i ||
         $vcf_info->{'left_flank'} =~ /n/i ||
         $vcf_info->{'right_flank'} =~ /n/i ||
         $left_flank_length < $flank_size ||
         $right_flank_length < $flank_size))
  {
    $num_ns_in_alleles_or_flanks++;
  }
  else
  {
    # No filtering or passed filtering
    $num_printed++;
    $vcf->print_entry($vcf_entry);
  }
}

my @missing_chr_names = sort keys %missing_chrs;

if(@missing_chr_names > 0)
{
  print STDERR "vcf_add_flanks.pl: Missing chromosomes: " .
               join(", ", @missing_chr_names) . "\n";
}

if($num_ref_mismatch > 0)
{
  my $mismatch_percent = 100 * $num_ref_mismatch / $total_num_entries;

  print STDERR "vcf_add_flanks.pl: " .
               num2str($num_ref_mismatch) . " / " .
               num2str($total_num_entries) . " " .
               "(" . sprintf("%.2f", $mismatch_percent) . "%) " .
               "variants removed for not matching the reference\n";
}

if($num_ns_in_alleles_or_flanks > 0)
{
  my $ns_percent = 100 * $num_ns_in_alleles_or_flanks / $total_num_entries;

  print STDERR "vcf_add_flanks.pl: " .
               num2str($num_ns_in_alleles_or_flanks) . " / " .
               num2str($total_num_entries) . " " .
               "(" . sprintf("%.2f", $ns_percent) . "%) " .
               " variants removed for containing Ns in ref/alt allele or " .
               "flanks too short\n";
}

if($filter_by_ref_match || $filter_Ns)
{
  my $printed_percent = 100 * $num_printed / $total_num_entries;

  print STDERR "vcf_filter_by_checks.pl: " .
               num2str($num_printed) . " / " .
               num2str($total_num_entries) . " " .
               "(" . sprintf("%.2f", $printed_percent) . "%) " .
               "variants printed\n";
}

close($vcf_handle);
