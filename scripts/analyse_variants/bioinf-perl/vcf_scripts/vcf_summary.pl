#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use List::Util qw(min);

use VCFFile;
use UsefulModule;

## Config
my $csvsep = ",";
#

sub print_usage
{
  for my $err (@_)
  {
    print STDERR "Error: $err\n";
  }
  
  print STDERR
"Usage: ./vcf_summary.pl [--max_indel <abs_svlen>] [in.vcf]
  Print summary stats\n";

  exit;
}

my $max_indel;

if(@ARGV > 1 && $ARGV[0] =~ /^-?-max_indel$/i)
{
  shift;
  $max_indel = shift;

  if(!defined($max_indel))
  {
    print_usage("No max_indel value given");
  }
  elsif(!$max_indel =~ /^\d+$/)
  {
    print_usage("Max indel must be a positive integer greater than zero");
  }
}
 
my $vcf_file = shift;

if(@ARGV > 0)
{
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

my $num_of_variants = 0;
my $num_of_complex = 0;
my $num_of_clean_indels = 0;

my $num_of_transitions = 0;
my $num_of_transversions = 0;

my $num_of_small_indels = 0;
my $num_of_small_tandem = 0;
my $num_of_small_non_tandem = 0;

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  $num_of_variants++;

  my $ref = uc($vcf_entry->{'true_REF'});
  my $alt = uc($vcf_entry->{'true_ALT'});
  my $svlen = $vcf_entry->{'INFO'}->{'SVLEN'};
  my $indel = get_clean_indel($vcf_entry);

  if(length($ref) == 1 && length($alt) == 1)
  {
    my $snp = $ref.$alt;

    if($snp eq "AG" || $snp eq "GA" || $snp eq "CT" || $snp eq "TC")
    {
      $num_of_transitions++;
    }
    else
    {
      $num_of_transversions++;
    }
  }
  elsif(defined($indel))
  {
    $num_of_clean_indels++;

    if((!defined($max_indel) || abs($svlen) <= $max_indel))
    {
      $num_of_small_indels++;

      # tandem
      my $left_flank = $vcf_entry->{'INFO'}->{'left_flank'};
      my $right_flank = $vcf_entry->{'INFO'}->{'right_flank'};

      if(defined($left_flank) && defined($right_flank))
      {
        if(is_slippage($left_flank, $indel, $right_flank))
        {
          $num_of_small_tandem++;
        }
        else
        {
          $num_of_small_non_tandem++;
        }
      }
    }
  }
  else
  {
    # Not SNP or clean indel -> COMPLEX
    $num_of_complex++;
  }
}

close($vcf_handle);

# Print results
my $num_of_snps = $num_of_transitions+$num_of_transversions;
print "SNPs:         ".pretty_fraction($num_of_snps, $num_of_variants)." of variants\n";
print "clean indels: ".pretty_fraction($num_of_clean_indels, $num_of_variants)." of variants\n";
print "complex:      ".pretty_fraction($num_of_complex, $num_of_variants)." of variants\n";

# Print Ts/Tv results
print "Ts:    ".pretty_fraction($num_of_transitions, $num_of_snps)." of snps\n" .
      "Tv:    ".pretty_fraction($num_of_transversions, $num_of_snps)." of snps\n";

if($num_of_transversions > 0)
{
  print "Ts/Tv: " .
        sprintf("%.3f", $num_of_transitions / $num_of_transversions) . "\n";
}

# small indel results
print "small clean indels: ".pretty_fraction($num_of_small_indels, $num_of_variants)." of variants\n";

my $total_tandem_nontandem = $num_of_small_tandem + $num_of_small_non_tandem;
print "small tandem:       ".pretty_fraction($num_of_small_tandem, $total_tandem_nontandem)." of small indels\n";
print "small non-tandem:   ".pretty_fraction($num_of_small_non_tandem, $total_tandem_nontandem)." of small indels\n";

sub is_slippage
{
  my ($left_flank, $indel, $right_flank) = @_;

  my $left_flank_rev = reverse($left_flank);
  my $indel_rev = reverse($indel);

  my $slippage_dist = compare_slippage_indel_flank($left_flank_rev, $indel_rev) +
                      compare_slippage_indel_flank($right_flank, $indel);

  my $is_slippage = ($slippage_dist >= length($indel));

  return $is_slippage;
}

sub compare_slippage_indel_flank
{
  my ($str1, $str2) = @_;

  $str1 = uc($str1);
  $str2 = uc($str2);

  my $len = min(length($str1), length($str2));

  for(my $i = 0; $i < $len; $i++)
  {
    if(substr($str1,$i,1) ne substr($str2,$i,1))
    {
      return $i;
    }
  }

  return $len;
}
