#!/usr/bin/perl

use strict;
use warnings;

if(@ARGV != 1)
{
    print "usage: ./split_fasta.pl <ref.fa>\n";
    print "  Splits a single fasta file into a file for each sequence\n";
    exit;
}

my $file = shift;

open(FILE, $file) or die("Cannot open file '$file'");

my $name = "";
my $seq = "";

my $line;

while(defined($line = <FILE>))
{
    if($line =~ /^>/)
    {
	if(length($name) > 0)
	{
	    # Save FASTA file
	    #print "Saving $name\n";
	    open(OUT, '>'.$name.'.fa') or die("Cannot write fasta file '$name'");

	    print OUT ">$name\n";
	    print OUT "$seq";

	    close(OUT);
	}

	chomp($line);

	# Remove '>' from name
	$name = substr($line,1);
	$seq = "";
    }
    else
    {
	$seq .= $line;
    }
}

if(length($name) > 0)
{   
    # Save FASTA file
    #print "Saving name\n";
    open(OUT, '>'.$name.'.fa') or die("Cannot write fasta file '$name'");

    print OUT ">$name\n";
    print OUT "$seq";

    close(OUT);
}


close(FILE);
