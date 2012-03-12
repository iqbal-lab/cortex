#!/usr/bin/perl

use strict;
use warnings;

## Config
my $csvsep = ",";
##

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./fasta_acgtn.pl [files..]\n";
  print STDERR "  Replace non-ACGTN bases with N\n";
  exit;
}

my @files = @ARGV;

my @read_stdin = grep {$_ eq "-"} @files;

if(@read_stdin > 1)
{
  print_usage("Cannot read from stdin more than once (too many '-'s)");
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

  parse_file($handle);

  close($handle);
}

sub parse_file
{
  my ($handle) = @_;
  my $line;
  while(defined($line = <$handle>))
  {
    chomp($line);

    if($line !~ /^>/) {
      $line =~ s/[^acgtn]/N/gi;
    }

    print "$line\n";
  }
}
