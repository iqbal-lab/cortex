#!/usr/bin/perl -w
use strict;

my $file = shift;

open(FILE, $file)||die("Cannot open $file");
open(OUT, ">".$file.".clean")||die("Cannot open $file.clean");

my $count_good = 0;
my $count_bad = 0;
while (<FILE>)
{
    my $line = $_;
    my $flag = 0;

    my $readid_5p;
    my $seq5p;
    my $readid_br1;
    my $seqbr1;
    my $readid_br2;
    my $seqbr2;
    my $readid_3p;
    my $seq3p;


    if ($line =~ /^>var_5p/)
    {
	$readid_5p=$line;
	
	$seq5p =<FILE>;
	if ($seq5p =~ /N/)
	{
	    $flag=1;
	}	
	$readid_br1 = <FILE>;
	$seqbr1 =<FILE>;
	if ($seqbr1 =~ /N/)
	{
	    $flag=1;
	}	
	$readid_br2 = <FILE>;
	$seqbr2 =<FILE>;
	if ($seqbr2 =~ /N/)
	{
	    $flag=1;
	}	
	$readid_3p = <FILE>;
	$seq3p =<FILE>;
	if ($seq3p =~ /N/)
	{
	    $flag=1;
	}	

    }
    else
    {
	die("probem parsing $line");
    }
    if ($flag==0)
    {
	$count_good++;
	print OUT "$readid_5p$seq5p$readid_br1$seqbr1$readid_br2$seqbr2$readid_3p$seq3p";
    }
    else
    {
	$count_bad++;
    }
}
close(FILE);
close(OUT);

print "Good Count is $count_good\n";
print "Bad Count is $count_bad\n";
