#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(min max sum);

use FASTNFile;
use UsefulModule;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./fastn_qualities.pl [files..]
  Summarise FASTQ quality data.\n";

  exit;
}

my @files = @ARGV;

my @read_stdin = grep {$_ eq "-"} @files;

if(@read_stdin > 1)
{
  print STDERR "Warning: reading from stdin more than once (multiple '-'s)\n";
}

my $min_quality = 0;
my $max_quality = 0;
my $sum_qualities = 0;
my $num_of_bases = 0;

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

# Print results
print num2str($num_of_bases) . " bases\n";
print num2str($min_quality) . " min quality\n";
print num2str($max_quality) . " max quality\n";
print num2str($sum_qualities / $num_of_bases,undef,2) . " mean quality\n";


#
# Subroutines
#

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

  my ($title, $seq, $qualities);

  while((($title, $seq, undef, $qualities) = $fastn_file->read_next()) && defined($title))
  {
    my $seq_len = length($seq);
    $num_of_bases += $seq_len;

    for(my $i = 0; $i < $seq_len; $i++)
    {
      my $qual = ord(substr($qualities, $i, 1)) - 33;
      $min_quality = min($qual, $min_quality);
      $max_quality = max($qual, $max_quality);
      $sum_qualities += $qual;
    }
  }
}
