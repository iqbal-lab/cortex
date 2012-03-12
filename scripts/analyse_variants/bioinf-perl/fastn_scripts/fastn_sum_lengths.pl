#!/usr/bin/perl

use strict;
use warnings;

use UsefulModule;

use List::Util qw(max);

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./fastn_lengths.pl [file]\n";
  print STDERR " Print the length of each sequence\n";
  print STDERR " If file not given (or '-') reads from STDIN\n";
  exit;
}

if(@ARGV == 1 && $ARGV[0] =~ /^--?h(elp)?$/i)
{
  print_usage();
}

#
# Open FASTA/Q handle
#
my $bases_read = 0;

if(@ARGV == 0)
{
  my $fastn_handle = open_stdin();
  
  print_read_lengths($fastn_handle);
  
  close($fastn_handle);
}
else
{
  for my $file (@ARGV)
  {
    my $fastn_handle;

    if($file eq "-")
    {
      $fastn_handle = open_stdin();
    }
    else
    {
      open($fastn_handle, $file)
        or print_usage("Cannot open FASTA/Q file '$file'");
    }
    
    print_read_lengths($fastn_handle);
    
    close($fastn_handle);
  }
}

print "$bases_read\n";

sub open_stdin
{
  my $stdin_handle;

  if(-p STDIN) {
    # STDIN is connected to a pipe
    open($stdin_handle, "<&=STDIN") or print_usage("Cannot read STDIN pipe");
  }
  else
  {
    print_usage("Must specify or pipe in a FASTA/Q file");
  }
  
  return $stdin_handle;
}

sub print_read_lengths
{
  my ($handle) = @_;

  my $line = <$handle>;
  my $is_fastq = 0;

  if($line =~ /^>/) {
    $is_fastq = 0;
  }
  elsif($line =~ /^@/) {
    $is_fastq = 1;
  }
  else {
    print_usage("Cannot identify file as fasta or fastq");
  }

  chomp($line);

# FASTQ example:
#@WTCHG_1025:6:1:1020:19385#0
#CAACTTGAAATCTGAAGAACAGATGAGTGATCCTTAATAGTTTTCTTTTCN
#+
####################################################

  if($is_fastq)
  {
    my $name = $line;

    while(defined($name))
    {
      $line = <$handle>; # read sequene line

      chomp($line);

      $bases_read += length($line);
  
      $line = <$handle>; # read '+'
      $line = <$handle>; # read quality scores
      $name = <$handle>; # read next entry @...
    }
  }
  else
  {
    while(defined($line = <$handle>))
    {
      chomp($line);

      if($line !~ /^>/)
      {
        $bases_read += length($line);
      }
    }
  }
}
