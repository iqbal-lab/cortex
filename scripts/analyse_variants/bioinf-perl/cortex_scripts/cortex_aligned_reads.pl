#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(min max sum);

use CortexCovgFile;

#
# Isaac Turner <isaac.turner@dtc.ox.ac.uk>
# 23 July 2012
#

## Config
my $csvsep = ",";
#

sub print_usage
{
  for my $err (@_)
  {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./cortex_aligned_reads.pl [.colour_covgs]
  Print out read names that have coverage >0 in any kmer in any colour\n";

  exit;
}

if(@ARGV > 1)
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

# Start reading aligned entries
my ($read_name, $sequence, $colours_arrref);

while((($read_name, $sequence, $colours_arrref) = $covgfile->read_align_entry()) &&
      defined($read_name))
{
  for my $col_covgs (@$colours_arrref)
  {
    if(max(@$col_covgs) > 0)
    {
      print "$read_name\n";
      #print join(" ", @$col_covgs) . "\n";
      #print "$sequence\n";
    }
  }
}

close($covg_handle);


