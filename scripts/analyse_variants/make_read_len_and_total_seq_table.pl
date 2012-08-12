#!/usr/bin/perl -w
use strict;

my $file = shift; #log file with read length and covg info

get_data($file);


sub get_data
{
    my ($log) = @_;
    open(LOG, $log)||die();

    my $outer_done=0;
    while ($outer_done==0)
    {
	my $ln = <LOG>;
	chomp $ln;
	if ($ln =~ /SUMMARY/)
	{
	    $ln = <LOG>; #header, of form 
	    # Colour  SampleID        MeanReadLen     TotalSeq        ErrorCleaning   LowCovSupsThresh        LowCovNodesThresh       PoolagainstWhichCleaned
	    my $readlen_column=-1;
	    my $seq_column=-1;
	    my @sp_head = split(/\t/, $ln);
	    my $zam;
	    for ($zam=0; $zam<scalar(@sp_head); $zam++)
	    {
		if ($sp_head[$zam] eq "MeanReadLen")
		{
		    $readlen_column=$zam;
		}
		if ($sp_head[$zam] eq "TotalSeq")
		{
		    $seq_column=$zam;
		}
	    }
	    if ( ($readlen_column==-1) || ($seq_column ==-1) )
	    {
		die("Tell Zam, make_read_len_and_total_seq_table.pl has been broken by his recent changes.\n");
	    }
	    my $done = 0;
	    while ($done==0)
	    {
		$ln = <LOG>;
		chomp $ln;
		if ($ln !~ /\*/)
		{
		    my @sp = split(/\t/, $ln);
		    my $colour = $sp[0];
		    my $readlen = $sp[$readlen_column];
		    my $seq = $sp[$seq_column];
		    print "$colour\t$readlen\t$seq\n";
		}
		else
		{
		    $done=1;
		    $outer_done=1;
		}
	    }

	}
	
    }
    close(LOG);
}
