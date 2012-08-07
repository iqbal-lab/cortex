#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use VCFFile;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_print_breakpoint.pl <flank_size> [in.vcf]\n";

  exit;
}

## Test for filtering
my $skip_failed_vars = 0;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  $skip_failed_vars = 1;
}
##

if(@ARGV < 1 || @ARGV > 2)
{
  print_usage();
}

my $flank_size = shift;
my $vcf_file = shift;

if($flank_size !~ /^\d+$/ || $flank_size <= 0)
{
  print_usage("Invalid flank size '$flank_size'");
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
if($skip_failed_vars) { $vcf->set_filter_failed(undef); }

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  my $aa = $vcf_entry->{'INFO'}->{'AA'};

  if(defined($aa) && ($aa eq "0" || $aa eq "1"))
  {
    my $left_flank = $vcf_entry->{'INFO'}->{'left_flank'};
    my $right_flank = $vcf_entry->{'INFO'}->{'right_flank'};

    if(!defined($left_flank) || !defined($right_flank))
    {
      print_usage("Missing left/right flank '".$vcf_entry->{'ID'}."'");
    }
    elsif(length($left_flank) < $flank_size || length($right_flank) < $flank_size)
    {
      print_usage("Flank is too short on variant '".$vcf_entry->{'ID'}."'");
    }

    print ">".$vcf_entry->{'ID'}."\n";
    print substr($left_flank, -$flank_size) .
          $vcf_entry->{$aa eq "0" ? 'true_REF': 'true_ALT'} .
          substr($right_flank, 0, $flank_size) . "\n";
  }
}

close($vcf_handle);
