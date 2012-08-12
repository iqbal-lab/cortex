#!/usr/bin/perl -w
use strict;

my $file = shift;

## change from flank, branch, branch, flank on 4 lines to
## just the flanks

open(FILE, $file)||die();
while(<FILE>)
{
    my $line=$_;
    chomp $line;
    if ($line =~ /^>var_(\d+)_5p_flank/)
    {
	
	my $readid_5p = $line;
	chomp $readid_5p;
	if ($readid_5p =~ /^(\S+)/)
	{
	    $readid_5p=$1;
	}


	my $f5 = <FILE>;

	chomp $f5;
	my $len_f5 = length($f5);

	if ($len_f5>1000)
	{
	    $readid_5p = $readid_5p."_cut_at_1000";
	    $f5 = substr($f5, -1000);
	    $len_f5 = 1000;
	}
	
	<FILE>;<FILE>;<FILE>;<FILE>;

	my $readid_3p = <FILE>;
	chomp $readid_3p;
	if ($readid_3p =~ /^(\S+)/)
	{
	    $readid_3p=$1;
	}


	my $f3 = <FILE>;
	chomp $f3;
	my $len_f3 = length($f3);

	if ($len_f3>1000)
	{
	    $readid_3p = $readid_3p."_cut_at_1000";
	    $f3 = substr($f3, 0,1000);
	    $len_f3 = 1000;
	}

	
	print "$readid_5p\n";
	print "$f5\n";
    }
}
close(FILE);
