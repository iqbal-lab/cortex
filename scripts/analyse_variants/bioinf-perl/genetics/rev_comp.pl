#!/usr/bin/perl

use warnings;
use strict;

use GeneticsModule;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./rev_comp.pl [seq]\n";
  print STDERR "  Reverse complements a DNA sequence\n";
  exit;
}

my $seq;

if(@ARGV == 0 && (-p STDIN))
{
  # STDIN is connected to a pipe
  open(IN, "<&=STDIN") or print_usage("Cannot read pipe");
  $seq = <IN>;
  close(IN);
}
elsif(@ARGV == 1)
{
  $seq = shift;
}
else {
  print_usage();
}

chomp($seq);

print rev_comp($seq)."\n";
