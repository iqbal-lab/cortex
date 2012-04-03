#!/usr/bin/perl

use strict;
use warnings;

if(@ARGV < 2 || @ARGV > 4)
{
  print "usage: ./clean_fasta.pl [-p|-n] <ref.fa> <out.fa>\n";
  print "  Copies sequence from <ref.fa> to <out.fa> only if name is not chrUn* or *random*\n";
  print "  -p  Print sequence names that are copied\n";
  print "  -n  Print sequence names that are not copied\n";
  exit;
}

my $print_clean = 0;
my $print_unclean = 0;

while(@ARGV > 2)
{
  my $arg = shift;

  if($arg =~ /^-n$/i) {
    $print_unclean = 1;
  }
  elsif($arg =~ /^-p$/i) {
    $print_clean = 1;
  }
  else {
    print "  Invalid argument '$arg'\n";
    exit;
  }
}

my $file_in = shift;
my $file_out = shift;

open(IN, $file_in) or die("Cannot open input fasta '$file_in'");
open(OUT, '>'.$file_out) or die("Cannot write fasta file '$file_out'");

my $name = "";
my $seq = "";
my $printing = 0;

my $line;

while(defined($line = <IN>))
{
  if($line =~ /^>/)
  {
    if(length($line) > 1 && $line !~ /^>chrUn/ && $line !~ /random/)
    {
      # Save this sequence to the output FASTA file
      print OUT $line;
      $printing = 1;

      if($print_clean) {
        print "+".substr($line,1);
      }
    }
    else
    {
      $printing = 0;

      if($print_unclean) {
        print "-".substr($line,1);
      }
    }
  }
  elsif($printing)
  {
    print OUT $line;
  }
}

close(OUT);
close(IN);
