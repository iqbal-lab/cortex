#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use FASTNFile;
use RefGenome;
use UsefulModule;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./fastn_select.pl [-d|-f <names.txt>|-r <readname>] [in1.fa in2.fq ..]
  Prints reads with given names, from FASTA/Q files (or STDIN if '-').

   --allow_dupes|-d     if multiple entries have the same name
   --file|-f <f.pos>    is a file with a readname on each line
   --read|-r <readname> specify a read name to print\n";

  exit;
}

if(@ARGV < 1)
{
  print_usage();
}

my @search_files = ();
my %read_names = ();
my $allow_dupes = 0;

while(@ARGV >= 1)
{
  if(@ARGV >= 2 && $ARGV[0] =~ /^-?-f(ile)?$/i)
  {
    shift;
    push(@search_files, shift);
  }
  elsif(@ARGV >= 2 && $ARGV[0] =~ /^-?-r(ead)?$/i)
  {
    shift;
    my $readname = shift;
    $read_names{$readname} = 1;
  }
  elsif($ARGV[0] =~ /^-?-allow_dupes?$/i || lc($ARGV[0]) eq "-d")
  {
    shift;
    $allow_dupes = 1;
  }
  else
  {
    last;
  }
}

my @ref_files = @ARGV;

# Check we can read ref files

for my $ref_file (@ref_files)
{
  if($ref_file ne '-' && !(-r $ref_file))
  {
    print_usage("Cannot read fasta/q file '$ref_file'");
  }
}

# Load searches from files
for my $search_file (@search_files)
{
  open(SEARCH, $search_file)
    or print_usage("Couldn't open file of coords '$search_file'");

  my $search_line;

  while(defined($search_line = <SEARCH>))
  {
    chomp($search_line);
    $read_names{$search_line} = 1;
  }

  close(SEARCH);
}


# Check we have some read names to look for
if(scalar(keys %read_names) == 0)
{
  print_usage("No read names given");
}

#
# Open FASTA/Q handles
#
if(@ref_files == 0)
{
  push(@ref_files, '-');
}

for my $file (@ref_files)
{
  my $fastn_handle;

  if($file eq "-")
  {
    $fastn_handle = open_stdin("Cannot read file -- need to pipe in fasta/fastq");
  }
  else
  {
    open($fastn_handle, $file)
      or print_usage("Cannot open FASTA/Q file '$file'");
  }
  
  search_file($fastn_handle);
  
  close($fastn_handle);
  
  if(keys(%read_names) == 0)
  {
    # All entries have been found - no need to open any more files
    last;
  }
}

# Done

if(keys(%read_names) > 0)
{
  # Print warnings for search results with no match
  for my $read_name (keys %read_names)
  {
    if(!$allow_dupes || $read_names{$read_name} == 1)
    {
      print STDERR "fastn_substr.pl: Couldn't find '$read_name'\n";
    }
  }
}


#
# Functions
#

sub search_file
{
  my ($handle) = @_;

  my $fastn = new FASTNFile($handle);

  # Read in
  my ($name, $seq, $qual);

  while((($name, $seq, undef, $qual) = $fastn->read_next()) && defined($name))
  {
    if(defined($read_names{$name}))
    {
      if($allow_dupes)
      {
        $read_names{$name} = 0;
      }
      else
      {
        delete($read_names{$name});

        if(scalar(keys(%read_names)) == 0)
        {
          last;
        }
      }

      $fastn->print_entry($name, $seq, $qual);
    }
  }

  close($handle);
}
