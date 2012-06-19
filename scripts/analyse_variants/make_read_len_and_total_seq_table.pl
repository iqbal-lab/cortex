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
	    <LOG>;
	    my $done = 0;
	    while ($done==0)
	    {
		$ln = <LOG>;
		chomp $ln;
		if ($ln !~ /\*/)
		{
		    my @sp = split(/\t/, $ln);
		    my $colour = $sp[0];
		    my $readlen = $sp[1];
		    my $seq = $sp[2];
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
