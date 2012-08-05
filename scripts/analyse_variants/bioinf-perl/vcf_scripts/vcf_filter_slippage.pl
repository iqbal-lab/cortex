#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(min max);

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

  print STDERR "" .
"Usage: ./vcf_filter_slippage.pl [OPTIONS] [in.vcf]
  Prints and flags slippage indels.  Slippage can only be on clean indels.
  Requires INFO tags 'left_flank' and 'right_flank'

  --only <tandem|non>  Print only tandem/non-tandem events (+ SNPs etc)
  --non_flag <tag>     INFO flag non-tandem (slippage < SVLEN)
  --tandem_flag <tag>  INFO flag tandem (slippage >= SVLEN)
  --dist_tag <tag>     INFO tag with amount of slippage (bp shared with flank)\n";
  exit;
}

## Test for filtering
my $failed_vars_out = undef;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  open($failed_vars_out, ">-");
}
##

my ($filter, $tandem_flag, $non_tandem_flag, $dist_tag);

while(@ARGV > 0)
{
  if($ARGV[0] =~ /^-?-only$/i)
  {
    shift;
    $filter = lc(shift);

    if(defined($filter) && $filter =~ /^-*non[_ \-]?tandem$/i)
    {
      $filter = "non";
    }

    if(!defined($filter) || ($filter ne "tandem" && $filter ne "non"))
    {
      print_usage("--only needs to be given 'tandem' or 'non'");
    }
  }
  elsif($ARGV[0] =~ /^-?-tandem_flag$/i)
  {
    shift;
    $tandem_flag = shift;
    
    if(!defined($tandem_flag))
    {
      print_usage("--tandem_flag needs a name to be specified");
    }
  }
  elsif($ARGV[0] =~ /^-?-non_flag$/i)
  {
    shift;
    $non_tandem_flag = shift;
    
    if(!defined($non_tandem_flag))
    {
      print_usage("--non_flag needs a name to be specified");
    }
  }
  elsif($ARGV[0] =~ /^-?-dist_tag$/i)
  {
    shift;
    $dist_tag = shift;
    
    if(!defined($dist_tag))
    {
      print_usage("--dist_flag needs a name to be specified");
    }
  }
  elsif(@ARGV > 1)
  {
    print_usage("Unknown argument '$ARGV[0]'");
  }
  else
  {
    last;
  }
}

my $vcf_file = shift;

if(@ARGV > 0)
{
  print_usage();
}
elsif(defined($filter))
{
  if(defined($tandem_flag) && $filter eq "tandem")
  {
    print_usage("Did you mean to flag tandems and filter them out??");
  }
  if(defined($non_tandem_flag) && $filter eq "non")
  {
    print_usage("Did you mean to flag non-tandems and filter them out??");
  }
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

# Add tags to VCF header
if(defined($dist_tag))
{
  $vcf->add_header_tag("INFO", $dist_tag, 1, "Integer", "Distance of slippage");
}

if(defined($tandem_flag))
{
  $vcf->add_header_tag("INFO", $tandem_flag, 0, "Flag",
                       "Tandem indels (distance of slippage >= length)");
}

if(defined($non_tandem_flag))
{
  $vcf->add_header_tag("INFO", $non_tandem_flag, 0, "Flag",
                       "Non-tandem indels (distance of slippage < length)");
}

$vcf->print_header();

my $num_of_variants = 0;
my $num_of_indels = 0;
my $num_of_tandem = 0;
my $num_of_printed = 0;

my $total_slippage = 0;

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  $num_of_variants++;

  my $ref = $vcf_entry->{'true_REF'};
  my $alt = $vcf_entry->{'true_ALT'};

  my $long_allele = $ref;
  my $short_allele = $alt;

  if(length($long_allele) < length($short_allele))
  {
    # swap
    ($long_allele, $short_allele) = ($short_allele, $long_allele);
  }

  my $is_tandem = 0;
  my $slippage_dist = 0;

  if(length($short_allele) == 0 && length($long_allele) > 0)
  {
    $num_of_indels++;

    if(!defined($vcf_entry->{'INFO'}->{'left_flank'}) ||
       !defined($vcf_entry->{'INFO'}->{'right_flank'}))
    {
      warn($vcf_entry->{'ID'}.": left_flank / right_flank INFO tags not set");
    }

    my $left_flank_rev = reverse($vcf_entry->{'INFO'}->{'left_flank'});
    my $right_flank = $vcf_entry->{'INFO'}->{'right_flank'};

    my $indel = $long_allele;
    my $indel_rev = reverse($indel);

    $slippage_dist = get_match($left_flank_rev, $indel_rev) +
                     get_match($right_flank, $indel);

    $total_slippage += $slippage_dist;

    if($slippage_dist >= length($indel))
    {
      $is_tandem = 1;
      $num_of_tandem++;
    }
  }

  if(defined($tandem_flag) && $is_tandem)
  {
    $vcf_entry->{'INFO_flags'}->{$tandem_flag} = 1;
  }
  elsif(defined($non_tandem_flag) && !$is_tandem)
  {
    $vcf_entry->{'INFO_flags'}->{$non_tandem_flag} = 1;
  }

  if(defined($dist_tag))
  {
    $vcf_entry->{'INFO'}->{$dist_tag} = $slippage_dist;
  }

  if(!defined($filter) || ($filter eq "tandem") == $is_tandem)
  {
    $vcf->print_entry($vcf_entry);
    $num_of_printed++;
  }
}

# Print
print STDERR "vcf_filter_slippage.pl: " .
              pretty_fraction($num_of_printed, $num_of_variants) . " " .
              "variants printed";

if(defined($filter))
{
  print STDERR " (".($filter eq "non" ? "" : "non-")."tandem indels removed)";
}

print STDERR "\n";

print STDERR "vcf_filter_slippage.pl: " .
              pretty_fraction($num_of_indels, $num_of_variants) . " " .
              "variants were clean indels\n";

if($num_of_indels > 0)
{
  print STDERR "vcf_filter_slippage.pl: " .
                pretty_fraction($num_of_tandem, $num_of_indels) . " " .
                "indels are tandem\n";

  print STDERR "vcf_filter_slippage.pl: " .
                pretty_fraction($num_of_indels-$num_of_tandem, $num_of_indels) .
                " indels are non-tandem\n";

  print STDERR "vcf_filter_slippage.pl: mean slippage per indel is " .
               sprintf("%.3f", $total_slippage / $num_of_indels) . "bp\n";
}

close($vcf_handle);

sub get_match
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
