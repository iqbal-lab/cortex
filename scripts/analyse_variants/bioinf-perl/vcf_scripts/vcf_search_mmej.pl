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
"Usage: ./vcf_search_mmej.pl  [in.vcf]\n" .
"  Prints (or flags) MMEJ indels. \n".
"  Requires INFO tags 'left_flank' and 'right_flank'\n" .
"  \n" .
"  --tag <tag>    INFO tag to add MMEJ\n";
  exit;
}

if(@ARGV > 4)
{
  print_usage();
}

my $flag_lh_or_rh_indel;
my $flag_lh_rg;


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


# Add tag to VCF header
my $flag1 = "MMEJ_FLANK_INDEL";
my $description1 = "MMEJ-like indels. At least one flank shares >5bp homology (measured as length of longest word shared) with the indel";
$vcf->add_header_tag("INFO", $flag1, 0, "Flag", $description1);
my $flag2 = "MMEJ_FLANK_FLANK";
my $description2 = "How much homology do the flanks of MMEJ-like indels share. Length of longest word shared";
$vcf->add_header_tag("INFO", $flag2, 0, "Flag", $description2);


$vcf->print_header();

my $num_of_variants = 0;
my $num_of_printed = 0;

my $vcf_entry;
my $is_microhom_lh_indel;
my $is_microhom_rh_indel;
my $is_microhom_lh_rh;
my $best_microhom_lh_indel;
my $best_microhom_rh_indel;
my $best_microhom_lh_rh;

while(defined($vcf_entry = $vcf->read_entry()))
{
    my $type =  $vcf_entry->{'INFO'}->{'SVTYPE'};
    if ($type =~ /COMPLEX/)
    {
	next;
    }
    my $ref = $vcf_entry->{'true_REF'};
    my $alt = $vcf_entry->{'true_ALT'};

  my $long_allele = $ref;
  my $short_allele = $alt;

  if(length($long_allele) < length($short_allele))
  {
    # swap
    ($long_allele, $short_allele) = ($short_allele, $long_allele);
  }

  if (  (length($long_allele)==1) && (length($short_allele)==1) )
  {
      ## is a SNP
      next;
  }
    


  $num_of_variants++;


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

    $best_microhom_lh_indel = get_micro($vcf_entry->{'INFO'}->{'left_flank'}, $indel);
    $best_microhom_rh_indel = get_micro($vcf_entry->{'INFO'}->{'right_flank'}, $indel);
    $best_microhom_lh_rh    = get_micro($vcf_entry->{'INFO'}->{'left_flank'}, 
					$vcf_entry->{'INFO'}->{'right_flank'});

    #$slippage_dist = get_match($left_flank_rev, $indel_rev) +
    #                 get_match($right_flank, $indel);
    
    $is_microhom_lh_indel=0;
    $is_microhom_rh_indel=0;
    $is_microhom_lh_rh=0;
    if ($best_microhom_lh_indel>=5)
    {
	$is_microhom_lh_indel=1;
    }
    if ($best_microhom_rh_indel>=5)
    {
	$is_microhom_rh_indel=1;
    }
    if ($best_microhom_lh_rh>=5)
    {
	$is_microhom_lh_rh=1;
    }

  }


  if (($is_microhom_lh_indel==1) || ($is_microhom_rh_indel==1))
  {
      $vcf_entry->{'INFO'}->{$flag1} = get_max($best_microhom_lh_indel, $best_microhom_rh_indel);
  }
  if ($is_microhom_lh_rh==1)
  {
      $vcf_entry->{'INFO'}->{$flag2} = $best_microhom_lh_rh;
  }
  $vcf->print_entry($vcf_entry);
  if (($is_microhom_lh_indel==1) || ($best_microhom_rh_indel==1) || ($is_microhom_lh_rh==1) )
  {
      $num_of_printed++;
  }

}

# Print filtered rate
print STDERR "vcf_search_mmej.pl: " .
              pretty_fraction($num_of_printed, $num_of_variants) . " " .
              "variants printed\n";

close($vcf_handle);

sub get_max
{
    my ($a, $b) = @_;
    if ($a>=$b)
    {
	return $a;
    }
    return $b;
}
sub get_micro
{
    my ($str1, $str2) = @_;
    my $cmd = "/home/zam/dev/hg/CORTEX_release/scripts/analyse_variants/needleman_wunsch-0.3.0/needleman_wunsch --zam $str1 $str2";
    my $ret = qx{$cmd};
    my @sp = split(/\n/, $ret);
    if ( ($sp[0] !~ /Br1/) || ($sp[2] !~ /Br2/) )
    {
	die("Cannot parse $ret");
    }
    my $res = $sp[1];

    my @sp2 = split(//, $res);
    my $i=0;
    my $best=0;
    my $current=0;
    for ($i=0; $i<scalar(@sp2); $i++)
    {
	if ($sp2[$i] eq "|")
	{
	    $current++;
	}
	else
	{
	    if ($current>$best)
	    {
		$best = $current;
	    }
	    $current=0;
	}
    }
    if ($current>$best)
    {
	$best = $current;
    }
    
    return $best;
}
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
