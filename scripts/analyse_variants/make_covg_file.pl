#!/usr/bin/perl -w
use strict;


my $callfile = shift; ## callfile output by Cortex
my $number_of_colours = shift;
my $ref_colour = shift;## if one of the colours is the reference, specify here. If none if them is a reference, set this to -1.


my $outfile = $callfile.".covg_for_classifier";
if (-e $outfile)
{
    die("Halting make_covg_file.pl as the output file $outfile already exists. Do you really mean to overwrite it? Delete it first and rerun, or save it somewhere else\n");
}

open(CALLS, $callfile) or die("Cannot open file '$callfile'");
open(OUT, ">$outfile") or die("Cannot open file '$outfile'");


#print OUT "VAR\tREF_FILTER\tBR1_LEN\tBR2_LEN\tREF_COV_BR1\tREF_COV_BR2\tSAMPLES...\n";

my $name = "";
my $reads_on_branch=0;

while (<CALLS>)
{
    my $line = $_;
    
    if ($line =~ /^\>(\S*var_\d+)_5p_flank/)
    {
	$name = $1;
	my $lenbr1;
	my $lenbr2;

	print OUT "$name\t";
	my @array_counts_for_each_colour=(); ## will be an array of strings of the form "1/14", "2,3", etc where "a,b" at index i means colour i has count a on branch1  and count b on branch 2

	my $all_nodes_br1_in_ref=0;
	my $all_nodes_br2_in_ref=0;

	$line = <CALLS>;
	$line = <CALLS>;
	$line = <CALLS>;## br1 seq
	chomp $line;
	$lenbr1 = length($line);
	$line = <CALLS>;
	$line = <CALLS>;
	chomp $line;
	$lenbr2 = length($line);
	$line = <CALLS>;
	$line = <CALLS>;


	# what do we expect next:
	#Colour  br1_median_covg br2_median_covg br1_min_covg    br2_min_covg
	#0       0       0       0       0
	#1       91      86      67      70
	#2       73      71      51      62

	while ($line !~ /^Colour/)
	{
	    	$line = <CALLS>;
	}
	my $n=0;
	while ($n<=$number_of_colours-1)
	{
	    $line = <CALLS>;
	    my @sp = split(/\s+/, $line);
	    my $col = $sp[0];
	    if ($n!=$col)
	    {
		die("Unexpected ordering of colours in $line\n");
	    }
	    my $br1_median = $sp[1];
	    my $br2_median = $sp[2];
	    my $br1_min = $sp[3];
	    my $br2_min = $sp[4];
	    $array_counts_for_each_colour[$col] = $br1_median."\t".$br2_median;
	    if ($col == $ref_colour)
	    {
		if ($br1_min>0)
		{
		    $all_nodes_br1_in_ref=1;
		}
		if ($br2_min>0)
		{
		    $all_nodes_br2_in_ref=1;
		}
	    }
	    $n++;
	}

	if ( ($all_nodes_br1_in_ref==1) && ($all_nodes_br2_in_ref==1) )
	{
	    print OUT "REF_BUBBLE\t";
	}
	else
	{
	    print OUT "NOT_IN_REF\t";
	}
	print OUT $lenbr1."\t".$lenbr2."\t";
	if ($ref_colour>=0)
	{
	    print OUT $array_counts_for_each_colour[$ref_colour];
	    print OUT "\t";
	}
	elsif ($ref_colour==-1)
	{
	    #Zam - just modified this. If there is no ref in the graph, printnothing. Either way, num of columns we print =  2*number of colours, but if 
	    ## there is a ref colour, that is moved to the front
	    #print OUT "0";
	}
	else
	{
	    die("Ref colour is not -1 and yet is not positive\n");
	}
	my $i;

	for ($i=0; $i<$number_of_colours; $i++)
	{
	    if ($i==$ref_colour)
	    {
		if ($i==$number_of_colours-1)
		{
		    print OUT "\n";
		}
		next;
	    }
	    print OUT $array_counts_for_each_colour[$i];
	    if ($i<$number_of_colours-1)
	    {
		print OUT "\t";
	    }
	    else
	    {
		print OUT "\n";
	    }
	}

    }


}
close(CALLS);
close(OUT);
print "Completed making covg file!\n";
