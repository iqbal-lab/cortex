#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use CortexCovgFile;

sub print_usage
{
  for my $err (@_) {
    print "Error: $err\n";
  }

  print STDERR "usage: ./print_5p_flank.pl [--trim <bp>] [.colour_covgs]\n";
  print STDERR "  Print a FASTA of 5p_flank from each bubble\n";
  print STDERR "  If [.colour_covgs] omitted (or '-') read from STDIN\n";
  print STDERR "   --trim <t> trim sequences longer than <t> (takes right most bases)\n";
  
  exit;
}

if(@ARGV > 3)
{
  # Too many arguments
  print_usage();
}

my $trim;

if(@ARGV >= 2)
{
  if($ARGV[0] =~ /^--?t(rim)?$/i)
  {
    shift;
    $trim = shift;
    
    if($trim !~ /^\d+$/ || $trim == 0)
    {
      print_usage("--trim must be a positive integer >= 1");
    }
  }
  else
  {
    print_usage("Unknown argument '$ARGV[0]'");
  }
}

my $covgs_file = shift;

#
# Open .colour_covgs Handle
#
my $covgs_handle;

if(defined($covgs_file) && $covgs_file ne "-") {
  open($covgs_handle, $covgs_file)
    or die("Cannot open .colour_covgs file '$covgs_file'\n");
}
elsif(-p STDIN) {
  # STDIN is connected to a pipe
  open($covgs_handle, "<&=STDIN") or die("Cannot read pipe");
}
else
{
  print_usage("Must specify or pipe in a .colour_covgs file");
}

# Read outgroup covgs ref entry
my $cortex = new CortexCovgFile($covgs_handle);

my ($flank_5p, $flank_3p, $branches, $col_llk) = $cortex->read_bubble_entry();

while(defined($flank_5p))
{
  my $seq = $flank_5p->{'seq'};

  if(defined($trim) && length($seq) > $trim)
  {
    $seq = substr($seq, -$trim);
  }

  print ">" . $flank_5p->{'name'} . "\n";
  print $seq . "\n";

  ($flank_5p, $flank_3p, $branches, $col_llk) = $cortex->read_bubble_entry();
}

close($covgs_handle);

