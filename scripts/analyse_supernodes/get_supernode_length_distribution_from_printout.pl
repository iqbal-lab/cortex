#!/usr/bin/perl -w
use strict;

my $file = shift; ## file of dumped/printed out supernodes
my $bin_size = shift;
my $max_expected_sup_length=shift;


open(FILE,$file)||die("Cannot open $file");

my %length_bin_to_freq=();
my $i;
for ($i=0; $i< $max_expected_sup_length; $i=$i+$bin_size)
{
    $length_bin_to_freq{$i}=0;
}

while(<FILE>)
{
    my $line = $_;
    chomp $line;
    
    if ($line =~ /^\>/)
    {
	if ($line =~ /length:(\d+)/)
	{
	    my $len = $1;
	    if ($len>$max_expected_sup_length)
	    {
		die("This node $line is longer than max expected length you entered as arg 3: $max_expected_sup_length");
	    }
	    
	    ($length_bin_to_freq{ $bin_size* int($len/$bin_size)}) ++;

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


print "Supernode length frequencies, binsize $bin_size:\n";

for ($i=0; $i< $max_expected_sup_length; $i=$i+$bin_size)
{
    print "$i\t";
    print $length_bin_to_freq{$i};
    print "\n";
}
