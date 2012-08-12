#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use CortexCovgFile;

## Config
my $csvsep = ",";
#

sub print_usage
{
  for my $err (@_)
  {
    print STDERR "Error: $err\n";
  }
  
  print STDERR "Usage: ./cortex_snp_count.pl <kmer_size> [.colour_covgs]\n";
  print STDERR "  Print number of SNPs as well as MNPs of each size\n";
  exit;
}

if(@ARGV < 1 || @ARGV > 2)
{
  print_usage();
}

my $kmer_size = shift;
my $covg_file = shift;

if($kmer_size !~ /^\d+$/ || $kmer_size <= 0)
{
  print_usage("<kmer_size> value invalid '$kmer_size'");
}

#
# Open .colour_covgs handle
#
my $covg_handle;

if(defined($covg_file) && $covg_file ne "-")
{
  open($covg_handle, $covg_file)
    or print_usage("Cannot open .colour_covgs file '$covg_file'\n");
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

my @mnps_sizes = ();

# Start reading bubbles
my ($flank_5p, $flank_3p, $branches) = $covgfile->read_bubble_entry();

while(defined($flank_5p))
{
  my $branch1_len = length($branches->[0]->{'seq'});
  my $branch2_len = length($branches->[1]->{'seq'});

  if($branch1_len == $branch2_len)
  {
    my $size = $branch2_len - $kmer_size;
    $mnps_sizes[$size]++;
  }

  ($flank_5p, $flank_3p, $branches) = $covgfile->read_bubble_entry();
}

close($covg_handle);

print "NPsize".$csvsep."count\n";
for(my $i = 1; $i < @mnps_sizes; $i++)
{
  print $i . $csvsep . (defined($mnps_sizes[$i]) ? $mnps_sizes[$i] : 0) . "\n";
}
