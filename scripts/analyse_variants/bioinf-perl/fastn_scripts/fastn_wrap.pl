#!/usr/bin/perl

use strict;
use warnings;

use FASTNFile;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./fastn_wrap.pl <line_wrap> <FASTA|FASTQ> [files..]
  Wrap lines.  line_width of 0 means no wrap.  Takes FASTA/Q input files.  You
  need to specify the format of the output you require.\n";

  exit;
}

if(@ARGV < 2)
{
  print_usage();
}

my $line_wrap = shift;
my $file_output = uc(shift);
my @files = @ARGV;

if($line_wrap !~ /^\d+$/)
{
  print_usage("Line wrap is not a positive integer (0 for no wrap)");
}
elsif($file_output ne "FASTA" && $file_output ne "FASTQ")
{
  print_usage("File output format is neither FASTA nor FASTQ");
}

my @read_stdin = grep {$_ eq "-"} @files;

if(@read_stdin > 1)
{
  print STDERR "Warning: reading from stdin more than once (multiple '-'s)\n";
}

if(@files == 0)
{
  # Read from STDIN
  push(@files, '-');
}

for my $file (@files)
{
  my $fastn_handle;

  if($file eq "-")
  {
    $fastn_handle = open_stdin("Cannot read ref -- need to pipe in fasta/fastq");
  }
  else
  {
    open($fastn_handle, $file)
      or print_usage("Cannot open FASTA/Q file '$file'");
  }

  parse_file($fastn_handle);

  close($fastn_handle);
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

  my ($title, $seq, $qual);

  if($file_output eq "FASTA")
  {
    while((($title, $seq) = $fastn_file->read_next()) && defined($title))
    {
      print_FASTA($title, $seq, $line_wrap);
    }
  }
  elsif($file_output eq "FASTQ")
  {
    while((($title, $seq, $qual) = $fastn_file->read_next()) && defined($title))
    {
      print_FASTQ($title, $seq, $qual, $line_wrap);
    }
  }
}
