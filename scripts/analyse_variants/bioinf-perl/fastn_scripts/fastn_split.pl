#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(min max);
use File::Basename;

use FASTNFile;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./fastn_split.pl <-b <bases>|-e <entries>> [file]
  Split fasta/q file into files of <bases> bp or <entries> seqeuences. If file
  not given (or '-') reads from STDIN\n";

  exit;
}

if(@ARGV < 2 || @ARGV > 3)
{
  print_usage();
}

my $arg = shift;
my $limit = shift;

my $count_bases_not_entries;

if($arg =~ /^-?-b(ases)?$/i)
{
  $count_bases_not_entries = 1;
}
elsif($arg =~ /^-?-e(ntries?)?$/i)
{
  $count_bases_not_entries = 0;
}
else
{
  print_usage("Unknown argument '$arg' - expected -b <limit> or -e <limit>");
}

if($limit !~ /^\d+$/ || $limit == 0)
{
  print_usage("Invalid limit '$limit' -- positive integer expected");
}

my $fastn_file = shift;

#
# Open FASTA/Q Handle
#
my $fastn_fh;

if(defined($fastn_file) && $fastn_file ne "-")
{
  open($fastn_fh, $fastn_file)
    or print_usage("Cannot open FASTA/Q file '$fastn_file'\n");
}
elsif(-p STDIN) {
  # STDIN is connected to a pipe
  open($fastn_fh, "<&=STDIN") or print_usage("Cannot read pipe");
}
else
{
  print_usage("Must specify or pipe in a FASTA/Q file");
}

# Create FASTNFile
my $fastn = new FASTNFile($fastn_fh);

# Get filename
if(defined($fastn_file))
{
  $fastn_file = basename($fastn_file);
}
else
{
  my $ext = $fastn->is_fastq() ? ".fq" : ".fa";
  $fastn_file = "split".$ext;
  
  for(my $i = 1; -e $fastn_file; $i++)
  {
    $fastn_file = "split-".$i.$ext;
  }
}

my $out_counter = -1;
my $out_fh;

my ($title, $seq, $qual);

my $file_counter = 0;

while((($title, $seq, undef, $qual) = $fastn->read_next()) && defined($title))
{
  if($file_counter >= $limit || !defined($out_fh))
  {
    # Next file
    if(defined($out_fh))
    {
      close($out_fh);
    }

    $file_counter = 0;
    $out_counter++;
    my $out_file = $fastn_file.".".$out_counter;

    open($out_fh, ">$out_file")
      or print_usage("Cannot open output '$out_file'");
  }

  $fastn->print_entry($title, $seq, $qual, 0, $out_fh);

  if($count_bases_not_entries)
  {
    $file_counter += length($seq);
  }
  else
  {
    $file_counter++;
  }
}

if(defined($out_fh))
{
  close($out_fh);
}
