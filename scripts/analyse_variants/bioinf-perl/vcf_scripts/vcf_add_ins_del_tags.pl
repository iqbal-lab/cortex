#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use VCFFile;
use List::Util qw(min max);

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./vcf_add_ins_del_tags.pl [file.vcf]
  Adds INFO tags using AA and SVLEN INFO tags:
   * INDEL=INS|DEL (for clean insertions / deletions)
   * AALEN=(derived-ancestral allele length)\n";

  exit;
}

## Test for filtering
my $failed_vars_out = undef;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  open($failed_vars_out, ">-");
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

# Print non-PASS variants straight to stdout if -p passed
if(defined($failed_vars_out)) { $vcf->set_filter_failed($failed_vars_out);}

# Add tag and print header
my $tag_description = "Indel mutation: INS=insertion, DEL=deletion or .=unknown";
$vcf->add_header_tag("INFO", "INDEL", 1, "String", $tag_description);
$vcf->print_header();

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  my $ref_len = length($vcf_entry->{'true_REF'});
  my $alt_len = length($vcf_entry->{'true_ALT'});
  
  if(min($ref_len, $alt_len) == 0 && max($ref_len, $alt_len) > 0)
  {
    # SVLEN is length(alt) - length(ref)
    if(!defined($vcf_entry->{'INFO'}->{'AA'}))
    {
      $vcf_entry->{'INFO'}->{'INDEL'} = '.';
      $vcf_entry->{'INFO'}->{'AALEN'} = '.';
    }
    elsif($vcf_entry->{'INFO'}->{'AA'} eq "0")
    {
      # ref allele is ancestral
      my $ins = ($vcf_entry->{'INFO'}->{'SVLEN'} > 0);
      $vcf_entry->{'INFO'}->{'INDEL'} = $ins ? 'INS' : 'DEL';
      $vcf_entry->{'INFO'}->{'AALEN'} = $alt_len - $ref_len;
    }
    elsif($vcf_entry->{'INFO'}->{'AA'} eq "1")
    {
      # alt allele is ancestral
      my $ins = ($vcf_entry->{'INFO'}->{'SVLEN'} < 0);
      $vcf_entry->{'INFO'}->{'INDEL'} = $ins ? 'INS' : 'DEL';
      $vcf_entry->{'INFO'}->{'AALEN'} = $ref_len - $alt_len;
    }
    else
    {
      $vcf_entry->{'INFO'}->{'INDEL'} = '.';
      $vcf_entry->{'INFO'}->{'AALEN'} = '.';
    }
  }

  $vcf->print_entry($vcf_entry);
}

close($vcf_handle);
