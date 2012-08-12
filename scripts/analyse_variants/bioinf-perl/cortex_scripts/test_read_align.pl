#!/usr/bin/perl

use strict;
use warnings;

use CortexCovgFile;

#
# Get distribution of indel sizes
#
# Isaac Turner <isaac.turner@dtc.ox.ac.uk>
# 13 Nov 2011

## Config
my $csvsep = ",";
#

sub print_usage
{
  for my $err (@_)
  {
    print STDERR "Error: $err\n";
  }
  
  print STDERR "Usage: ./test_read_align.pl [.colour_covgs]\n";
  print STDERR "  Print number of indels of each size\n";
  exit;
}

if(@ARGV > 2)
{
  print_usage();
}

my $covg_file = shift;

#
# Open .colour_covgs handle
#
my $covg_handle;

if(defined($covg_file) && $covg_file ne "-")
{
  open($covg_handle, $covg_file)
    or print_usage("Cannot open .colour_covgs file '$covg_file'");
}
elsif(-p STDIN) {
  # STDIN is connected to a pipe
  open($covg_handle, "<&=STDIN") or print_usage("Cannot read pipe");
}
else
{
  print_usage("Must specify or pipe in a Cortex .colour_covgs file");
}

my $covgfile = new CortexCovgFile($covg_handle);

# Start reading bubbles
my ($read_name, $sequence, $colours_hashref) = $covgfile->read_align_entry();

while(defined($read_name))
{
  print "$read_name ($sequence)\n";
  
  for my $col (sort keys %$colours_hashref)
  {
    print "  ".join(",",@{$colours_hashref->{$col}})."\n";
  }

  ($read_name, $sequence, $colours_hashref) = $covgfile->read_align_entry();
}

close($covg_handle);

