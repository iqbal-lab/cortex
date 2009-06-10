#!/usr/bin/perl -w
use strict;

my $list=shift;
open(LIST,$list)||die("cannot open $list");

while (<LIST>)
{
    my $filename=$_;
    chomp $filename;
    open(FILE,$filename)||die("Cannot open $filename");

    my $out="clean/".$filename;
    open(OUT, "> ".$out)||die("Cannot open $out");

    my $line_count=0;
    while(<FILE>)
    {
	my $line=$_;
	$line_count++;
	$line =~ s/\s+$/\n/;

	my @chars=split(//, $line);
	
	if ($chars[0] eq '>')
	{
	    print OUT $line;
	    next;
	}
	

	foreach my $c (@chars)
	{
	    if (($c ne "A") &&($c ne "a") &&($c ne "C") &&($c ne "c") &&($c ne "G") &&($c ne "g") &&($c ne "T") &&($c ne "t") && ($c ne "N") && ($c ne "n") && ($c ne "\n") )
	    {
		print "Bad character $c on line $line_count in file $filename. will replace with empty character\n";
		print "Line is\n$line\n";
		$line =~ s/$c//;
	    }
	}

       
	print OUT $line;

	
    }

    close(FILE);
    close(OUT);
}

close(LIST);

