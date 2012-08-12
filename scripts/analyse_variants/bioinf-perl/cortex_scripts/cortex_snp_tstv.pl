#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use CortexCovgFile;
use UsefulModule;

## Config
my $csvsep = ",";
#

sub print_usage
{
  for my $err (@_)
  {
    print STDERR "Error: $err\n";
  }
  
  print STDERR "" .
"Usage: ./cortex_snp_tstv.pl <kmer_size> [.colour_covgs]\n" .
"  Transition (Ts): A <-> G; Transversion (Tv) is all other substitutions\n" .
"  Although there are twice as many Tv as Tv combinations, Ts make up \n".
"  approximately 2/3rds of all SNPs[1].\n".
"  [1] Collins DW, Jukes TH (April 1994). 'Rates of transition and\n".
"  transversion in coding sequences since the human-rodent divergence'\n";

  exit;
}

if(@ARGV < 1 || @ARGV > 2)
{
  print_usage();
}

my $kmer_size = shift;
my $covg_file = shift;

if($kmer_size !~ /^\d+$/ || $kmer_size <= 0)
{
  print_usage("<kmer_size> value invalid '$kmer_size'");
}

#
# Open .colour_covgs handle
#
my $covg_handle;

if(defined($covg_file) && $covg_file ne "-")
{
  open($covg_handle, $covg_file)
    or print_usage("Cannot open .colour_covgs file '$covg_file'\n");
}
elsif(-p STDIN) {
  # STDIN is connected to a pipe
  open($covg_handle, "<&=STDIN") or print_usage("Cannot read pipe");
}
else
{
  print_usage("Must specify or pipe in a Cortex .colour_covgs file");
}

my $covgfile = new CortexCovgFile($covg_handle);

my $num_of_transitions = 0;
my $num_of_transversions = 0;

# Start reading bubbles
my ($flank_5p, $flank_3p, $branches) = $covgfile->read_bubble_entry();

while(defined($flank_5p))
{
  my $branch1_len = length($branches->[0]->{'seq'});
  my $branch2_len = length($branches->[1]->{'seq'});

  if($branch1_len == $kmer_size+1 && $branch2_len == $kmer_size+1)
  {
    my $base1 = uc(substr($branches->[0]->{'seq'}, 0, 1));
    my $base2 = uc(substr($branches->[1]->{'seq'}, 0, 1));

    my $snp = $base1.$base2;

    if($snp eq "AG" || $snp eq "GA" || $snp eq "CT" || $snp eq "TC")
    {
      $num_of_transitions++;
    }
    else
    {
      $num_of_transversions++;
    }
  }

  ($flank_5p, $flank_3p, $branches) = $covgfile->read_bubble_entry();
}

print "Transition (Ts) A<->G;  Transversion: all other SNPs\n";
print "Twice as many Tv values than Ts, but 2/3rds of SNPs are Ts\n";
print "Ts: ".num2str($num_of_transitions)."; " .
      "Tv: ".num2str($num_of_transversions)."\n";

if($num_of_transversions > 0)
{
  print sprintf("%.3f", $num_of_transitions / $num_of_transversions) .
        " Ts / Tv\n";
}

print "Total SNPs: ".num2str($num_of_transitions+$num_of_transversions)."\n";

close($covg_handle);
