#!/usr/bin/perl

use strict;
use warnings;

use FASTNFile;
use UsefulModule;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./fastn_gc.pl [files..]
  Prints GC content.\n";

  exit;
}

my @files = @ARGV;

my @read_stdin = grep {$_ eq "-"} @files;

if(@read_stdin > 1)
{
  print STDERR "Warning: reading from stdin more than once (multiple '-'s)\n";
}

my $gc_bases = 0;
my $n_bases = 0;
my $all_bases = 0;

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
my $non_n_bases = $all_bases - $n_bases;

print pretty_fraction($gc_bases, $all_bases) . " of all bases are GC.\n";
print pretty_fraction($gc_bases, $non_n_bases) . " of non-N bases are GC.\n";


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

  my ($title, $seq);

  while((($title, $seq) = $fastn_file->read_next()) && defined($title))
  {
    my $seq_len = length($seq);
    $all_bases += $seq_len;

    for(my $i = 0; $i < $seq_len; $i++)
    {
      my $base = uc(substr($seq, $i, 1));

      if($base eq "N")
      {
        $n_bases++;
      }
      elsif($base eq "G" || $base eq "C")
      {
        $gc_bases++;
      }
    }
  }
}
