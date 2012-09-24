#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(min max);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use VCFFile;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_print_indel_fasta.pl [file.vcf]\n";
  print STDERR "  Print FASTA of 5' branch, ancestral, alt allele & 3' branch\n";
  print STDERR "  On different lines.  Needs AA tag to be 0 or 1\n";
  exit;
}

## Test for filtering
my $skip_failed_vars = 0;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  $skip_failed_vars = 1;
}
##

if(@ARGV > 1)
{
  print_usage();
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

# Skip non-PASS variants if -p passed
if($skip_failed_vars) { $vcf->set_filter_failed(undef); }

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  # Works for multiple ALT frequencies separated by commas
  my $id = $vcf_entry->{'ID'};
  my $ref = $vcf_entry->{'true_REF'};
  my $alt = $vcf_entry->{'true_ALT'};

  my $aa = $vcf_entry->{'INFO'}->{'AA'};
  my $svlen = $vcf_entry->{'INFO'}->{'SVLEN'};
  my $min_allele_len = min(length($ref), length($alt));

  my $left = $vcf_entry->{'INFO'}->{'left_flank'};

  if(defined($left))
  {
    print ">".$id."_5p\n";
    print "$left\n";
  }

  if($svlen != 0 && $min_allele_len == 0)
  {
    # Clean Indel
    
    my $indel = "";
    if(defined($aa) && ($aa eq "0" || $aa eq "1"))
    {
      $indel = (($aa eq "0") == ($svlen > 0)) ? "ins" : "del";
    }
    else
    {
      $indel = "indel";
    }
    
    print ">".$id."_".$indel."\n";
    print "".($svlen > 0 ? $alt : $ref) . "\n";
  }
  elsif(defined($aa) && ($aa eq "0" || $aa eq "1"))
  {
    # ancestral info
    print ">".$id."_anc\n";
    print "".($aa eq "0" ? $ref : $alt) . "\n";
    print ">".$id."_alt\n";
    print "".($aa eq "0" ? $alt : $ref) . "\n";
  }
  else
  {
    # no ancestral
    print ">".$id."_ref\n";
    print "$ref\n";
    print ">".$id."_alt\n";
    print "$alt\n";
  }

  my $right = $vcf_entry->{'INFO'}->{'right_flank'};

  if(defined($right))
  {
    print ">".$id."_3p\n";
    print "$right\n\n";
  }
}

close($vcf_handle);
