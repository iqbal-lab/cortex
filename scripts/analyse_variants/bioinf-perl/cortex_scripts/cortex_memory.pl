#!/usr/bin/perl

use strict;
use warnings;

use POSIX qw(ceil);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use UsefulModule; # num2str

#
# Calculate amount of memory required by cortex
#
# Isaac Turner <isaac.turner@dtc.ox.ac.uk>
# 29 Jan 2011

sub usage
{
  for my $err (@_) {
    print "Error: $err\n";
  }
  
  print "Usage: ./cortex_memory.pl [options] <mem_height> <mem_width>\n";
  print "  Options: --kmer_size -k <k>   K-mer size of De Bruijn graph [default: 31]\n";
  print "           --colours -c <num>   Number of colours in graph    [default: 1]\n";
  exit;
}

if(@ARGV < 2 || @ARGV > 6 || @ARGV % 2 != 0) {
  usage();
}

my $memHeight = $ARGV[@ARGV-2];
my $memWidth = $ARGV[@ARGV-1];

my $kmerSize = undef;
my $numOfColours = undef;

# Loop through options
for(my $i = 0; $i < @ARGV-2; $i+=2)
{
  if($ARGV[$i] =~ m/^\-\-kmer\_size$/i || $ARGV[$i] =~ m/^\-k$/i)
  {
    $kmerSize = $ARGV[$i+1];

    if($kmerSize !~ /^[0-9]+$/ || $kmerSize < 1 || $kmerSize > 63) {
      usage("invalid kmer");
    }
  }
  elsif($ARGV[$i] =~ m/^\-\-colours$/i || $ARGV[$i] =~ m/^\-c$/i)
  {
    $numOfColours = $ARGV[$i+1];
    
    if($numOfColours !~ /^[0-9]+$/ || $numOfColours == 0) {
      usage("invalid number of colours");
    }
  }
  else
  {
    usage();
  }
}

if(!defined($kmerSize)) {
  $kmerSize = 31;
}

if(!defined($numOfColours)) {
  $numOfColours = 1;
}

if($memHeight !~ /^[0-9]+$/ || $memHeight < 1) {
  usage("invalid memHeight");
}

if($memWidth !~ /^[0-9]+$/ || $memWidth < 1) {
  usage("invalid memHeight");
}

print "width: $memWidth; height: $memHeight; colours: " .
      "$numOfColours; kmer_size: $kmerSize\n";

my $num_of_hash_entries = 2**$memHeight * $memWidth;

print num2str($num_of_hash_entries) . " hash table entries\n";

# Round entry size to nearest 8 bytes and multiply by the number of entries
my $bytes = $num_of_hash_entries *
            ceil((8*ceil($kmerSize/31) + 5*$numOfColours + 1)/8)*8;

print "Memory: " . num2str($bytes) . " bytes (" . mem2str($bytes) . ")\n";

print "Command: cortex_var_" . ($kmerSize > 31 ? "63" : "31") .
      "_c" . $numOfColours . " --kmer_size $kmerSize " .
      "--mem_width $memWidth --mem_height $memHeight\n";
