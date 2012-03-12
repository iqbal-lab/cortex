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
"Usage: ./vcf_filter_slippage.pl [--invert | --flag <flag>] [in.vcf]\n" .
"  Prints (or flags) slippage indels.  Slippage can only be on clean indels.\n".
"  Requires INFO tags 'left_flank' and 'right_flank'\n" .
"  \n" .
"  --invert       Print or flag non slippage indels\n" .
"  --tag <tag>    INFO tag to add slippage\n";
  exit;
}

if(@ARGV > 4)
{
  print_usage();
}

my $invert = 0;
my $flag;

while(@ARGV > 0)
{
  if($ARGV[0] =~ /^-?-invert$/i)
  {
    shift;
    $invert = 1;
  }
  elsif($ARGV[0] =~ /^-?-flag$/i)
  {
    shift;
    $flag = shift;
    
    if(!defined($flag))
    {
      print_usage("Flag name needs to be specified");
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

if(defined($flag))
{
  # Add tag to VCF header
  my $description = $invert ? "Not slippage indels" : "Slippage indels";
  $vcf->add_header_tag("INFO", $flag, 0, "Flag", $description);
}

$vcf->print_header();

my $num_of_variants = 0;
my $num_of_printed = 0;

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

  my $is_slippage = 0;
  my $slippage_dist = 0;

  if(length($short_allele) == 0 && length($long_allele) > 0)
  {
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
    
    $is_slippage = ($slippage_dist >= length($indel));
  }

  if(defined($flag))
  {
    # Print all, tag those that are filtered
    $vcf_entry->{'INFO'}->{$flag} = $slippage_dist;
    $vcf->print_entry($vcf_entry);
    
    if($is_slippage)
    {
      $num_of_printed++;
    }
  }
  elsif($is_slippage != $invert)
  {
    # Just print / filter
    $vcf->print_entry($vcf_entry);
    $num_of_printed++;
  }
}

# Print filtered rate
print STDERR "vcf_filter_slippage.pl: " .
              "(" . (!defined($flag) && $invert ? "not " : "") . "slippage) " .
              pretty_fraction($num_of_printed, $num_of_variants) . " " .
              "variants printed\n";

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
