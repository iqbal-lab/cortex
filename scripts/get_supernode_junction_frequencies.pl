#!/usr/bin/perl -w
use strict;

my $file = shift; ## file of dumped/printed out supernodes
my $total_num_supernodes = shift;

open(FILE,$file)||die("Cannot open $file");

my %junc_type_freq=();

$junc_type_freq{"00"}=0; #00 means no edges out at either end
$junc_type_freq{"01"}=0; #01 means one edge out at one end, and none at other end
$junc_type_freq{"02"}=0; # ..similarly
$junc_type_freq{"03"}=0;
$junc_type_freq{"04"}=0;
$junc_type_freq{"11"}=0;# These next 4 should be zero at the end. If there is one exit at one end, the supernode should be extended
$junc_type_freq{"12"}=0;#
$junc_type_freq{"13"}=0;
$junc_type_freq{"14"}=0;
$junc_type_freq{"22"}=0;
$junc_type_freq{"23"}=0;
$junc_type_freq{"24"}=0;
$junc_type_freq{"33"}=0;
$junc_type_freq{"34"}=0;
$junc_type_freq{"44"}=0;


while(<FILE>)
{
    my $line = $_;
    chomp $line;
    
    if ($line =~ /^\>/)
    {
	##for a supernode in a graph
	if ($line =~ /min_coverage:(\d+) max_coverage:(\d+)/)
	{
	    my $min=$1;
	    my $max=$2;

	    if ($max>$max_expected_occurrences)
	    {
		die("This node $line is has max covg $max, > than max expected occurrences you entered as arg 3: $max_expected_occurrences");
	    }
	    
	    ($occ_to_min_freq{$min}) ++;
	    ($occ_to_max_freq{$max}) ++;

	}
	else
	{
	    die("Bad format line: $line");
	}
    }
    else
    {
	##ignore sequence
    }

}
close(FILE);


print "Supernode occurrence frequencies:\n";

for ($i=0; $i< $max_expected_occurrences; $i++)
{
    print "$i\t";
    print $occ_to_min_freq{$i};
    print "\t";
    print $occ_to_max_freq{$i};    
    print "\n";
}
