#!/usr/bin/perl -w
use strict;


## Pass in filename of output from the trusted path/supernode SV caller.

my $svcalls=shift;
my $chrom = shift; ##which chromosome are these calls on?
my $min_flank_length_to_be_worth_mapping = shift;  ##Only print out a flank if it is longer than some min length


if (-e $svcalls)
{
}
else
{
    die("This file $svcalls does not exist");
}



open(IN, $svcalls)||die("Cannot open $svcalls");




my $current_output_fname = "chrom_".$chrom."_variant_flanks.fasta";

open(OUT, " >". $current_output_fname)||die("Cannot open $current_output_fname");

while (<IN>)
{
    my $line = $_;
    chomp $line;
    
    my $var_no;
    
    if ($line =~ /VARIATION: (\d+)/)
    {
	$var_no = $1;

	$line = <IN>;
	chomp $line;
	if ($line=~/5p_flank length:(\d+)/)
	{
	    my $flank5p_len = $1;
	    if ($flank5p_len >= $min_flank_length_to_be_worth_mapping)
	    {
		print OUT ">chromosome_".$chrom."_VARIATION_".$var_no."_5p_flank_length_".$flank5p_len."\n";
		$line = <IN>;
		print OUT $line; ##seq
	    }
	}
	else
	{
	    die("Problem with this line $line which should be 5p flank of variant $var_no");
	}
	##ignore variant and trusted branches
	<IN>;
	<IN>;
	<IN>;
	<IN>;

	$line = <IN>;
	chomp $line;
	if ($line=~/3p_flank length:(\d+)/)
	{
	    my $flank3p_len = $1;
	    if ($flank3p_len >=$min_flank_length_to_be_worth_mapping)
	    {
		print OUT ">chromosome_".$chrom."_VARIATION_".$var_no."_3p_flank_length_".$flank3p_len."\n";
		$line = <IN>;
		print OUT $line; ##seq
	    }
	}
	else
	{
	    die("Problem with this line $line which should be 3p flank of variant $var_no");
	}
	    

    }

}


close(IN);
close(OUT);

