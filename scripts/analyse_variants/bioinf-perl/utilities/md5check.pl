#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(reduce);

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"usage: ./md5check.pl <file1 ..>
  Compares multiple md5sum output lists.

  Example output:
    06ba2ae0897750fe14ac0b8b51b8cd4e  file1.txt.gz [2]
    Error: file2.txt.gz [2]
      60b725f10c9c85c70d97880dfe8191b3 [md5sums1.txt]
      b026324c6904b2a9cb4b88d6d61c81d1 [md5sums2.txt]
    26ab0db90d72e28ad0ba1e22ee510510  file3.txt.gz [2]

  Note: The number of hash entries for each file is shown in [x].
        If there is a mismatch between hashes, all hashes and their sources
        are shown.\n";
}

if(@ARGV < 1)
{
  print_usage();
}

my %file_hashes = ();

for my $file (@ARGV)
{
  open(FILE, $file) or die("Cannot open file '$file'");

  my $line;
  while(defined($line = <FILE>))
  {
    chomp($line);

    if($line !~ /^\s*$/)
    {
      my ($hash, $name) = split("  ", $line);
  
      if(!defined($file_hashes{$name}))
      {
        $file_hashes{$name} = {};
        $file_hashes{$name}->{$hash} = [$file];
      }
      elsif(!defined($file_hashes{$name}->{$hash}))
      {
        $file_hashes{$name}->{$hash} = [$file];
      }
      else
      {
        push(@{$file_hashes{$name}->{$hash}}, $file);
      }
    }
  }

  close(FILE);
}

my @files = sort keys %file_hashes;

for my $file (@files)
{
  my @hashes_for_file = keys %{$file_hashes{$file}};

  my $num_of_entries = 0;
  map {$num_of_entries += scalar(@{$file_hashes{$file}->{$_}})} @hashes_for_file;

  if(@hashes_for_file > 1)
  {
    print "Error: " . $file . " [$num_of_entries]\n";
    for my $hash (@hashes_for_file)
    {
      print "  $hash [" . join(" ", @{$file_hashes{$file}->{$hash}}) . "]\n";
    }
  }
  else
  {
    # File only has one hash - show the number of md5sum files that list it
    print "$hashes_for_file[0]  $file [$num_of_entries]\n";
  }
}
