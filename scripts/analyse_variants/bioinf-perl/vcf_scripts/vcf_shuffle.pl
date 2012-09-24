#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(shuffle);
use Fcntl qw(SEEK_SET);

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_shuffle.pl <in.vcf>\n";

  exit;
}

if(@ARGV != 1)
{
  print_usage();
}

my $vcf_file = shift;

my $vcf_handle;

open($vcf_handle, $vcf_file)
  or print_usage("Cannot open VCF input file '$vcf_file'\n");

my $line;
my $pos = 0;
my @line_starts = ();

while(defined($line = <$vcf_handle>) && substr($line, 0, 1) eq '#')
{
  print $line;
  $pos += length($line);
}

while(defined($line))
{
  if($line !~ /^\s*$/)
  {
    push(@line_starts, $pos);
  }

  $pos += length($line);
  $line = <$vcf_handle>;
}

@line_starts = shuffle(@line_starts);

for my $line_start (@line_starts)
{
  seek($vcf_handle, $line_start, SEEK_SET);
  $line = <$vcf_handle>;
  print $line;
}

close($vcf_handle);
