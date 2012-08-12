#!/usr/bin/perl

# Redundant
# Now moved into ~/perl/lib/FastaModule.pm

use strict;
use warnings;

use UsefulModule;
use FASTNFile;

if(@ARGV != 1)
{
  print "Usage: ./fastq_size.pl <file>\n";
  print " Estimate the number of reads in a file\n";
  exit;
}

my $file = shift;

my ($read_length, $file_size, $num_of_reads, $num_of_bases)
  = estimate_fastq_size($file);

print "first read: " . num2str($read_length) . " bases long\n";
print "file size: " . num2str($file_size) . " bytes\n";
print "estimated number of reads: " . num2str($num_of_reads) . "\n";
print "estimated number of bases: " . num2str($num_of_bases) . "\n";
