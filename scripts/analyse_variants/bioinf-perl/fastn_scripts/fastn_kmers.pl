#!/usr/bin/perl

use strict;
use warnings;

use UsefulModule;

## Config
my $csvsep = ",";
##

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./fastn_kmers.pl <K> [files..]\n";
  print STDERR "  Prints the occurances of each possible kmer of length k\n";
  print STDERR "  If file not given, or equals '-', reads from stdin)\n";
  exit;
}

if(@ARGV < 1)
{
  print_usage();
}

my $kmer_size = shift;
my @files = @ARGV;

if(@files == 1 && $files[0] eq "-") {
  @files = ();
}

if($kmer_size !~ /^\d+$/ || $kmer_size <= 0)
{
  print_usage("kmer_size not a valid integer greater than zero ('$kmer_size')");
}

my %kmers = ();

if(@files == 0)
{
  # Add stdin
  push(@files, "-");
}

# Read from files
for my $file (@files)
{
  my $handle;

  if($file eq "-")
  {
    $handle = open_stdin("Cannot read file -- need to pipe in fasta/fastq");
  }
  else
  {
    open($handle, $file)
      or print_usage("Cannot open FASTA/Q file '$file'");
  }

  load_file($handle);
  
  close($handle);
}

#
# Print kmers
#
print "kmer".$csvsep."count\n";

for my $kmer (sort {$a cmp $b} keys %kmers)
{
  print $kmer.$csvsep.$kmers{$kmer}."\n";
}


## FASTQ example:
#@WTCHG_1025:6:1:1020:19385#0
#CAACTTGAAATCTGAAGAACAGATGAGTGATCCTTAATAGTTTTCTTTTCN
#+
####################################################

## FASTA example:
#>chr1
#acgtnnnn
#tcga
#>chr2

sub load_file
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
    print_usage("Cannot identify file type (line: '$line')");
  }

  chomp($line);

  if($is_fastq)
  {
    while(defined($line))
    {
      $line = <$handle>; # read sequene line
      chomp($line); # remove new line characters
      $line = uc($line); # uppercase

      # Count all kmers
      for(my $i = 0; $i <= length($line)-$kmer_size; $i++)
      {
        $kmers{substr($line,$i,$kmer_size)}++;
      }
  
      $line = <$handle>; # read '+'
      $line = <$handle>; # read quality scores
      $line = <$handle>; # read next entry @...
    }

  }
  else
  {
    my $last_kmer = '';

    while(defined($line = <$handle>))
    {
      chomp($line);

      if($line =~ /^>/)
      {
        $last_kmer = '';
      }
      else
      {
        $line = uc($line); # uppercase

        if(length($last_kmer) > 1)
        {
          # Count kmers overhanging prev read line
          my $prev_kmer = $last_kmer;

          for(my $i = 1; $i < $kmer_size && $i <= length($line); $i++)
          {
            $last_kmer = substr($prev_kmer,$i) . substr($line,0,$i);
            $kmers{$last_kmer}++;
          }
        }
      
        # Count kmers on this line
        for(my $i = 0; $i <= length($line)-$kmer_size; $i++)
        {
          $last_kmer = substr($line,$i,$kmer_size);
          $kmers{$last_kmer}++;
        }
      }
    }
  }
}
