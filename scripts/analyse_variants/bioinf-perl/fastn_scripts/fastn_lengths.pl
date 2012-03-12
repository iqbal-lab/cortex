#!/usr/bin/perl

use strict;
use warnings;

use UsefulModule;

use List::Util qw(max);

my $csvsep = ",";

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./fastn_lengths.pl [--bases] [file]\n";
  print STDERR " Print the length of each sequence.  " .
               "If file not given (or '-') reads from STDIN\n";
  print STDERR " --bases only count [ACGT] bases, not `N's\n";
  exit;
}

my $bases_only = 0;

if(@ARGV == 1 && $ARGV[0] =~ /^-?-h(elp)?$/i)
{
  print_usage();
}
elsif(@ARGV > 0 && $ARGV[0] =~ /^-?-b(ases?)?$/i)
{
  shift;
  $bases_only = 1;
}
elsif(@ARGV > 1)
{
  print_usage();
}

print "chr,length\n";

#
# Open FASTA/Q handle
#
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
    my $name = substr($line,1);

    while(defined($name))
    {
      $line = <$handle>; # read sequene line

      chomp($name);
      chomp($line);

      print "$name" . $csvsep . count_length($line) . "\n";
  
      $line = <$handle>; # read '+'
      $line = <$handle>; # read quality scores
      $name = <$handle>; # read next entry @...
    }
  }
  else
  {
    my $name = substr($line,1);

    my $curr_length = 0;
  
    while(defined($line = <$handle>))
    {
      chomp($line);

      if($line =~ /^>/)
      {
        print $name.$csvsep.$curr_length."\n";

        $name = substr($line,1);
        $curr_length = 0;
      }
      else
      {
        $curr_length += count_length($line);
      }
    }

    print $name.$csvsep.$curr_length."\n";
  }
}

sub count_length
{
  my ($str) = @_;

  if($bases_only)
  {
    $str =~ s/[^acgt]//gi;
    return length($str);
  }
  else
  {
    return length($str);
  }
}
