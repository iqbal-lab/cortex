#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use VCFFile;

# Config #
my $csvsep = ",";
#

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_print_info_tags.pl <file.vcf> <infotag1 ..>\n";
  print STDERR "  Prints comma separated info tag values from VCF entries\n";
  print STDERR "  If <file.vcf> is '-' reads from STDIN\n";
  exit;
}

## Test for filtering
my $skip_failed_vars = 0;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  $skip_failed_vars = 1;
}
##

if(@ARGV < 2)
{
  print_usage();
}

my $vcf_file = shift;

my @tags = @ARGV;

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

# Print CSV header
print join($csvsep, @tags)."\n";

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  my @data = ();

  for my $tag (@tags)
  {
    my $d = defined($vcf_entry->{'INFO'}->{$tag}) ? $vcf_entry->{'INFO'}->{$tag} : "";
    push(@data, $d);
  }

  print join($csvsep, @data)."\n";
}

close($vcf_handle);
