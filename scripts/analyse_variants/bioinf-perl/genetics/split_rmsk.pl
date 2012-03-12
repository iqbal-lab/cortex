#!/usr/bin/perl

use warnings;
use strict;

use File::Path qw(mkpath);

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./split_rmsk.pl [rmsk.txt]\n";
  print STDERR "  Separates rmsk.txt file by repeat class\n";
  exit;
}

if(@ARGV > 1)
{
  print_usage();
}

my $rmsk_file = shift;

#
# Open rmsk.txt Handle
#
my $rmsk_handle;

if(defined($rmsk_file) && $rmsk_file ne "-")
{
  open($rmsk_handle, $rmsk_file)
    or print_usage("Cannot open RMSK file '$rmsk_file'\n");
}
elsif(-p STDIN)
{
  # STDIN is connected to a pipe
  open($rmsk_handle, "<&=STDIN") or print_usage("Cannot read pipe");
}
else
{
  print_usage("Must specify or pipe in a RMSK file");
}

mkpath('rmsk');

my %handles = ();

my $line; 
while(defined($line = <$rmsk_handle>))
{
  my @cols = split(/\t/, $line);
  my $class = $cols[11];

  if(!defined($handles{$class}))
  {
    my $file = "rmsk/rmsk.".$class.".txt";
    open($handles{$class}, ">$file") or die("Cannot open rmsk out '$file'");
  }

  my $handle = $handles{$class};

  print $handle $line;
}

close($rmsk_handle);

for my $handle (values %handles)
{
  close($handle);
}
