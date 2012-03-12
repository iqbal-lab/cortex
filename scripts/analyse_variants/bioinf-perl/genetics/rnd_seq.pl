#!/usr/bin/perl

use warnings;
use strict;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./rnd_seq.pl <len>\n";
  print STDERR "  Generates random DNA sequence\n";
  exit;
}

my $len;

if(@ARGV == 0 && (-p STDIN))
{
  # STDIN is connected to a pipe
  open(IN, "<&=STDIN") or print_usage("Cannot read pipe");
  $len = <IN>;
  close(IN);
}
elsif(@ARGV == 1)
{
  $len = shift;
}
else {
  print_usage();
}

chomp($len);

if($len !~ /^\d+$/)
{
  print_usage("Invalid length value '$len'");
}

my @bases = qw{a c g t};

for(my $i = 0; $i < $len; $i++)
{
  print $bases[int(rand() * 4)];
}

print "\n";
