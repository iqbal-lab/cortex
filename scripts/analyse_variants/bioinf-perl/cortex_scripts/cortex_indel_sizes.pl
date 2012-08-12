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
  
  print STDERR "Usage: ./cortex_indel_sizes.pl [--kmer <k>] [.colour_covgs]\n";
  print STDERR "  Print number of indels of each size\n";
  print STDERR "  --kmer <k>  only count 'clean' indels when called with k\n";
  exit;
}

my $k_limit;

if(@ARGV > 3)
{
  print_usage();
}
elsif(@ARGV >= 2)
{
  my $arg = shift;
  $k_limit = shift;
  
  if($arg !~ /^-?-kmer$/i)
  {
    print_usage("Unexpected argument '$arg'");
  }
  elsif($k_limit !~ /^\d+$/ || $k_limit <= 0)
  {
    print_usage("--kmer value invalid '$k_limit'");
  }
}

my $covg_file = shift;

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

my @indel_sizes = ();

# Start reading bubbles
my ($flank_5p, $flank_3p, $branches) = $covgfile->read_bubble_entry();

while(defined($flank_5p))
{
  my $branch1_len = length($branches->[0]->{'seq'});
  my $branch2_len = length($branches->[1]->{'seq'});

  if(!defined($k_limit) || $branch1_len <= $k_limit || $branch2_len <= $k_limit)
  {
    my $indel_size = abs($branch1_len - $branch2_len);
    $indel_sizes[$indel_size]++;
  }

  ($flank_5p, $flank_3p, $branches) = $covgfile->read_bubble_entry();
}

close($covg_handle);

print "indel_size".$csvsep."count\n";

my $i = 0;

while(!defined($indel_sizes[$i]) || $indel_sizes[$i] == 0)
{
  $i++;
}

for(; $i < @indel_sizes; $i++)
{
  print $i . $csvsep . (defined($indel_sizes[$i]) ? $indel_sizes[$i] : 0) . "\n";
}
