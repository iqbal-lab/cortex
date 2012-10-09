#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use List::Util qw(sum min max);

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
  Print summary stats.  Does not do any filtering.  
  To only use 'PASS' variants, filter first with vcf_filter_column.pl\n";

  exit;
}

## Test for filtering
my $skip_failed_vars = 0;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  $skip_failed_vars = 1;
}
##

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

# Skip non-PASS variants if -p passed
if($skip_failed_vars) { $vcf->set_filter_failed(undef);}

my $num_of_variants = 0;
my $num_of_complex = 0;
my $num_of_clean_indels = 0;

my $num_of_transitions = 0;
my $num_of_transversions = 0;

my $num_of_indels = 0;
my $num_of_tandem = 0;
my $num_of_non_tandem = 0;

# Ins/del ratios
my $num_of_del = 0;
my $num_of_ins = 0;
my @ins_sizes = ();
my @del_sizes = ();

my $num_of_tandem_del = 0;
my $num_of_tandem_ins = 0;
my @ins_sizes_tandem = ();
my @del_sizes_tandem = ();

my $num_of_non_tandem_del = 0;
my $num_of_non_tandem_ins = 0;
my @ins_sizes_non_tandem = ();
my @del_sizes_non_tandem = ();


my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  $num_of_variants++;

  my $ref = uc($vcf_entry->{'true_REF'});
  my $alt = uc($vcf_entry->{'true_ALT'});
  my $svlen = $vcf_entry->{'INFO'}->{'SVLEN'};
  my $indel = vcf_get_clean_indel($vcf_entry);

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
      $num_of_indels++;

      # tandem
      my $left_flank = $vcf_entry->{'INFO'}->{'left_flank'};
      my $right_flank = $vcf_entry->{'INFO'}->{'right_flank'};

      if(defined($left_flank) && defined($right_flank))
      {
        my $is_tandem = is_slippage($left_flank, $indel, $right_flank);

        if($is_tandem)
        {
          $num_of_tandem++;
        }
        else
        {
          $num_of_non_tandem++;
        }

        my $polarised_len;

        if(defined($polarised_len = $vcf_entry->{'INFO'}->{'AALEN'}) &&
           $polarised_len ne ".")
        {
          if($polarised_len > 0)
          {
            $num_of_ins++;
            push(@ins_sizes, $polarised_len);
          }
          else
          {
            $num_of_del++;
            push(@del_sizes, -$polarised_len);
          }

          if($is_tandem)
          {
            if($polarised_len > 0)
            {
              $num_of_tandem_ins++;
              push(@ins_sizes_tandem, $polarised_len);
            }
            else
            {
              $num_of_tandem_del++;
              push(@del_sizes_tandem, -$polarised_len);
            }
          }
          else
          {
            if($polarised_len > 0)
            {
              $num_of_non_tandem_ins++;
              push(@ins_sizes_non_tandem, $polarised_len);
            }
            else
            {
              $num_of_non_tandem_del++;
              push(@del_sizes_non_tandem, -$polarised_len);
            }
          }
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

@ins_sizes = sort {$a <=> $b} @ins_sizes;
@del_sizes = sort {$a <=> $b} @del_sizes;

@ins_sizes_tandem = sort {$a <=> $b} @ins_sizes_tandem;
@del_sizes_tandem = sort {$a <=> $b} @del_sizes_tandem;

@ins_sizes_non_tandem = sort {$a <=> $b} @ins_sizes_non_tandem;
@del_sizes_non_tandem = sort {$a <=> $b} @del_sizes_non_tandem;

# Print results
print "Note: includes ALL variants (PASS and those that fail any filters)\n";
print "-- Overview --\n";

my $num_of_snps = $num_of_transitions+$num_of_transversions;

print "SNPs:         " .
      pretty_fraction($num_of_snps, $num_of_variants) . " of variants\n";
print "clean indels: " . pretty_fraction($num_of_clean_indels, $num_of_variants) .
      " of variants\n";
print "complex:      " . pretty_fraction($num_of_complex, $num_of_variants) .
      " of variants\n";

# Print Ts/Tv results
print "-- SNPs --\n";
print "Ts:    ".pretty_fraction($num_of_transitions, $num_of_snps)." of snps\n" .
      "Tv:    ".pretty_fraction($num_of_transversions, $num_of_snps)." of snps\n";

if($num_of_transversions > 0)
{
  print "Ts/Tv: " .
        sprintf("%.3f", $num_of_transitions / $num_of_transversions) . "\n";
}

# small indel results
print "-- Indels --\n";
print "clean indels: " . pretty_fraction($num_of_indels, $num_of_variants) .
      " of variants\n";
print "  ins: ".pretty_fraction($num_of_ins, $num_of_ins+$num_of_del)."\n";
print "  ins distrib: ".print_arr_stats(\@ins_sizes)."\n";
print "  del distrib: ".print_arr_stats(\@del_sizes)."\n";

my $total_tandem_nontandem = $num_of_tandem + $num_of_non_tandem;

print "  tandem:       " .
      pretty_fraction($num_of_tandem, $total_tandem_nontandem) .
      " of small indels\n";

if($num_of_tandem > 0)
{
  my $polarised_tandem = $num_of_tandem_ins + $num_of_tandem_del;

  print "    ins: ".pretty_fraction($num_of_tandem_ins, $polarised_tandem)."\n";
  print "    ins distrib: ".print_arr_stats(\@ins_sizes_tandem)."\n";
  print "    del distrib: ".print_arr_stats(\@del_sizes_tandem)."\n";
}

print "  non-tandem:   " .
      pretty_fraction($num_of_non_tandem, $total_tandem_nontandem) . " " .
      "of small indels\n";

if($num_of_non_tandem > 0)
{
  my $polarised_non_tandem = $num_of_non_tandem_ins +
                             $num_of_non_tandem_del;
  print "    ins: " . pretty_fraction($num_of_non_tandem_ins,
                                      $polarised_non_tandem) . "\n";

  print "    ins distrib: " . print_arr_stats(\@ins_sizes_non_tandem) . "\n";
  print "    del distrib: " . print_arr_stats(\@del_sizes_non_tandem) . "\n";
}

sub print_arr_stats
{
  my ($arr) = @_;
  
  my $elements = scalar(@$arr);

  if($elements == 0)
  {
    return "";
  }

  my $mean = sum(@$arr) / $elements;
  my $median;

  if($elements == 1)
  {
    $median = $arr->[0];
  }
  elsif($elements % 2 == 0)
  {
    $median = ($arr->[$elements/2] + $arr->[$elements/2+1]) / 2;
  }
  else
  {
    $median = $arr->[$elements/2];
  }

  return "min: $arr->[0]; mean: " . sprintf("%.3f", $mean) . "; " .
         "median: " . sprintf("%.1f", $median) . "; max: $arr->[$elements-1];";
}

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
