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

if(@ARGV == 0)
{
  my $fastn_handle = open_stdin();
  
  get_max_read_in_fastn($fastn_handle, "stdin");
  
  close($fastn_handle);
}
else
{
  for my $file (@ARGV)
  {
    my $fastn_handle;
    my $descriptor;
  
    if($file eq "-")
    {
      $fastn_handle = open_stdin();
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
}

print "Total length (bp): ".num2str($total_read_length)."\n";
print "Longest read (bp): ".num2str($max_read_length)."\n";

sub open_stdin
{
  my $stdin_handle;

  if(-p STDIN) {
    # STDIN is connected to a pipe
    open($stdin_handle, "<&=STDIN") or die("Cannot read pipe");
  }
  else
  {
    print_usage("Must specify or pipe in a FASTA/Q file");
  }
  
  return $stdin_handle;
}

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

      ($title, $seq) = $fastn_file->read_next();
    }
  }
  else
  {
    # FASTA - memory effecient method:
    my ($line, $peak, $curr_length);
    
    while(defined($line = $fastn_file->read_line()))
    {
      my $seq_length = 0;

      while(defined($peak = $fastn_file->peak_line()) && $peak !~ /^[>]/)
      {
        $line = $fastn_file->read_line();
        chomp($line);
        $seq_length += length($line);
      }
  
      $max_read_length = max($seq_length, $max_read_length);
      $total_read_length += $seq_length;
    }
  }
}
