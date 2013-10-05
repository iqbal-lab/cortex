#!/usr/bin/perl -w
use strict;

use File::Basename;

my $cortex_dir = "/home/zam/dev/hg/bitbucket/CORTEX_mainline/";
my $mem_height = 24;
my $mem_width = 73;
my $pop="";
my $kmer= 31;
my $vcf=shift;
my $outfile = shift;



open(OUT, ">".$outfile)||die("Cannot open the specified output file: $outfile\n");   

open(VCF, $vcf)||die();
while (<VCF>)
{
    my $line = $_;
    if ($line !~ /^\#/)
    {
	my @sp = split(/\t/, $line);
	my $name = $sp[2];
	my $ref = $sp[3];
	my $alt = $sp[4];
	my $info = $sp[7];
	my $left;
	my $right;

	## left_flank=ATTAGATTTGACCTTCAGCAAGGTCAAAGGGAGTC;right_flank=GAACTAGTCTCAGGCTTCAACATCGAAT
	if ($info =~ /left_flank=(\S+);right_flank=(\S+)/)
	{
	    $left  = $1;
	    $right = $2;
	    if ( (length($left)<$kmer) || (length($right)<$kmer) )
	    {
		die("Your vcf is annotatedw with flanks, but they are shorter than the kmer\n");
	    }
	}

	my $lflank = substr($left, -$kmer);
	my $branch1= $ref.substr($right,0, $kmer);
	my $branch2= $alt.substr($right,0, $kmer);
	my $rflank = substr($right, -(length($right)-$kmer));
	print OUT ">".$name."_5p_flank\n$lflank\n".">".$name."_branch_1\n$branch1\n".">".$name."_branch_2\n".">".$name."_3p_flank$rflank\n";
    }
}
close(VCF);
