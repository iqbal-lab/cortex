#!/usr/bin/perl -w
use strict;

my $file = shift; ## file of dumped/printed out supernodes
my $max_expected_occurrences=shift; ## max number of times you expect any supernode to be seen in genome path. is "coverage" of the supernode in the graph

open(FILE,$file)||die("Cannot open $file");

##hash takes covg (ie number of times supernode exists in graph) to frequency (number of supernodes that have this coverage)
## Since graph is built of reads in general, not from a single fasta, a supernode need not have constant coverage/number of occurrences
my %occ_to_min_freq=();
my %occ_to_max_freq=();


my $i;
for ($i=0; $i< $max_expected_occurrences; $i++)
{
    $occ_to_min_freq{$i}=0;
    $occ_to_max_freq{$i}=0;
}

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
