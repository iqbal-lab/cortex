#!/usr/bin/perl

use strict;
use warnings;

use UsefulModule;
use FASTNFile;

use List::Util qw(max);

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "usage: ./fastn_total_length.pl [file1 ..]\n";
  print STDERR " Print the length of the longest sequence\n";
  print STDERR " If [file] is omitted '-', reads from STDIN\n";
  exit;
}

if(@ARGV > 0 && $ARGV[0] =~ /^--?h(elp)?$/i)
{
  print_usage();
}

#
# Open FASTA/Q handle
#

my $max_read_length = 0;
my $total_read_length = 0;
my $num_of_reads = 0;

my @files = @ARGV;

if(@files == 0)
{
  push(@files, '-');
}

for my $file (@files)
{
  my $fastn_handle;
  my $descriptor;

  if($file eq "-")
  {
    $fastn_handle = open_stdin("Cannot read file -- need to pipe in fasta/fastq");
    $descriptor = "stdin";
  }
  else
  {
    open($fastn_handle, $file) or die("Cannot open FASTA/Q file '$file'");
    $descriptor = $file;
  }
  
  get_max_read_in_fastn($fastn_handle, $descriptor);
  
  close($fastn_handle);
}

my $mean_read_length = $total_read_length / $num_of_reads;

print "Total length (bp): ".num2str($total_read_length)."\n";
print "Longest read (bp): ".num2str($max_read_length)."\n";
print "Number of reads:   ".num2str($num_of_reads)."\n";
print "Mean length (bp):  ".num2str(sprintf("%.2f", $mean_read_length))."\n";


sub get_max_read_in_fastn
{
  my ($handle, $descriptor) = @_;

  my $fastn_file = new FASTNFile($handle, $descriptor);

  if($fastn_file->is_fastq())
  {
    # FASTQ
    my ($title,$seq) = $fastn_file->read_next();

    while(defined($title))
    {
      my $seq_length = length($seq);

      $max_read_length = max($seq_length, $max_read_length);
      $total_read_length += $seq_length;
      $num_of_reads++;

      ($title, $seq) = $fastn_file->read_next();
    }
  }
  else
  {
    # FASTA - memory effecient method:
    my ($line, $peek, $curr_length);
    
    while(defined($line = $fastn_file->read_line()))
    {
      my $seq_length = 0;

      while(defined($peek = $fastn_file->peek_line()) && $peek !~ /^[>]/)
      {
        $line = $fastn_file->read_line();
        chomp($line);
        $seq_length += length($line);
      }
  
      $max_read_length = max($seq_length, $max_read_length);
      $total_read_length += $seq_length;
      $num_of_reads++;
    }
  }
}
