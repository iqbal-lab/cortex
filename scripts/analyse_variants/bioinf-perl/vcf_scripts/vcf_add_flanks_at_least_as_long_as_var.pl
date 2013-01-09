#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use VCFFile;
use RefGenome;
use UsefulModule; # num2str

use List::Util qw(min max);

#   Isaac Turner <isaac.turner\@dtc.ox.ac.uk> 2011/09/21

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./vcf_add_flanks.pl [FILTERS] <flank_size> <in.vcf> <ref_name> <ref1.fa ..>
  Adds the flanking sequence to VCF files
  --filter_ref_match: remove variants where ref allele doesn't match
  --filter_Ns: remove variants where flank+[ref/alt]+flank contains Ns or flanks
               aren't long enough
  - <ref_name> is used in the VCF header
  - If <in.vcf> is '-', reads from stdin\n";

  exit;
}

## Test for filtering
my $failed_vars_out = undef;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  open($failed_vars_out, ">-");
}
##

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

my $vcf = new VCFFile($vcf_handle);

# Print non-PASS variants straight to stdout if -p passed
if(defined($failed_vars_out)) { $vcf->set_filter_failed($failed_vars_out);}

#
# Load reference files
#
my $genome = new RefGenome(uppercase => 1);
$genome->load_from_files(@ref_files);

#
# Read VCF
#

# Add tags to header and print
my $tag_description = "$flank_size bp adjacent in ref genome '$ref_name'";
$vcf->add_header_tag("INFO", "left_flank", 1, "String", $tag_description);
$vcf->add_header_tag("INFO", "right_flank", 1, "String", $tag_description);
$vcf->print_header();

my $num_printed = 0;
my $num_flanks_added = 0;
my $num_ref_mismatch = 0;
my $num_ns_in_alleles_or_flanks = 0;
my $total_num_entries = 0;

my $vcf_entry;

my %missing_chrs = ();

while(defined($vcf_entry = $vcf->read_entry()))
{
  $total_num_entries++;
  my $chr = $vcf_entry->{'CHROM'};

  my $print = 1;

  if(!$genome->chr_exists($chr))
  {
    $missing_chrs{$chr} = 1;
  }
  else
  {
    # Get coordinates, convert to zero based
    my $var_start = $vcf_entry->{'true_POS'}-1;

    my $ref_allele = $vcf_entry->{'true_REF'};
    my $alt_allele = $vcf_entry->{'true_ALT'};
  
    my $var_length = length($ref_allele);
    my $alt_length = length($alt_allele);
    my $max_allele_len =$var_length;
    if ($alt_length>$var_length)
    {
	$max_allele_len =$alt_length;
    }
    my $ref_chrom_length = $genome->get_chr_length($chr);
    
    my $vcf_info = $vcf_entry->{'INFO'};

    if($var_start < 0 || $var_start + $var_length >= $ref_chrom_length)
    {
      print STDERR "vcf_add_flanks.pl: variant outside of reference genome " .
                   "bounds: " . $chr . ":" . $vcf_entry->{'POS'} . " " .
                   "(ref '".$genome->guess_chrom_fasta_name($chr)."' has " .
                   "length: " . num2str($ref_chrom_length) . ")\n";

      open(my $stderr, ">&", STDERR) or die("Cannot open STDERR");
      $vcf->print_entry($vcf_entry, $stderr);
    }
    else
    {
      ($vcf_info->{'left_flank'}, $vcf_info->{'right_flank'})
        = vcf_get_flanks($vcf_entry, $genome, get_max($flank_size, $max_allele_len));
  
      my $ref_seq = $genome->get_chr_substr($chr, $var_start, $var_length);

      if($filter_by_ref_match && $ref_seq ne uc($ref_allele))
      {
        $num_ref_mismatch++;
        $print = 0;
      }
      elsif($filter_Ns &&
            ($vcf_entry->{'REF'} =~ /n/i ||
             $vcf_entry->{'ALT'} =~ /n/i ||
             $vcf_info->{'left_flank'} =~ /n/i ||
             $vcf_info->{'right_flank'} =~ /n/i ||
             length($vcf_info->{'left_flank'}) < $flank_size ||
             length($vcf_info->{'right_flank'}) < $flank_size))
      {
        $num_ns_in_alleles_or_flanks++;
        $print = 0;
      }
      else
      {
        $num_flanks_added++;
        # No filtering or passed filtering
      }
    }
  }
  
  if($print)
  {
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

if($filter_by_ref_match)
{
  print STDERR "vcf_add_flanks.pl: " .
               pretty_fraction($num_ref_mismatch, $total_num_entries) . " " .
               "variants removed for not matching the reference\n";
}

if($filter_Ns)
{
  print STDERR "vcf_add_flanks.pl: " .
               pretty_fraction($num_ref_mismatch, $total_num_entries) . " " .
               "variants removed for not matching the reference\n";
}

print STDERR "vcf_add_flanks.pl: " .
             pretty_fraction($num_flanks_added, $total_num_entries) . " " .
             "variants had flanks added\n";

print STDERR "vcf_add_flanks.pl: " .
             pretty_fraction($num_printed, $total_num_entries) . " " .
             "variants printed\n";

close($vcf_handle);


sub get_max
{
    my ($a, $b) = @_;
    if ($a>$b)
    {
	return $a;
    }
    else
    {
	return $b;
    }
}
