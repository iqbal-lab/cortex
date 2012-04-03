#!/usr/bin/perl

use strict;
use warnings;

use FASTNFile;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./fastq2fastq.pl [files..]\n";
  print STDERR "  Convert FASTQ files to FASTA\n";
  exit;
}

my @files = @ARGV;

my @read_stdin = grep {$_ eq "-"} @files;

if(@read_stdin > 1)
{
  print STDERR "Warning: reading from stdin more than once (multiple '-'s)\n";
}

if(@files == 0)
{
  # Read from STDIN
  read_stdin();
}

for my $file (@files)
{
  if($file eq "-")
  {
    read_stdin();
  }
  else
  {
    my $handle;
    open($handle, $file) or print_usage("Cannot open fasta/q file '$file'");
    parse_file($handle);
    close($handle);
  }
}

sub read_stdin
{
  if(!(-p STDIN))
  {
    # STDIN is not connected to a pipe
    print_usage("Must specify or pipe in a FASTA/Q file");
  }

  my $handle;
  open($handle, "<&=STDIN") or die("Cannot read pipe");

  parse_file($handle, "stdin");

  close($handle);
}

sub parse_file
{
  my ($handle, $descriptor) = @_;

  my $fastn_file = new FASTNFile($handle, $descriptor);

  my ($title, $seq);

  while((($title, $seq) = $fastn_file->read_next()) && defined($title))
  {
    print_FASTA($title, $seq);
  }
}
