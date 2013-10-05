#!/usr/bin/perl -w
use strict;

my $cortex_root = "/home/zam/dev/CORTEX_release_v1.0.5.21/";  ## set this!
my $list = shift; ## list of RAW vcfs, one per population
my $ref = shift; ## reference fasta

open(LIST, $list)||die();
while (<LIST>)
{
    my $lyne = $_;
    chomp $lyne;
    my $newfile = parse_file($lyne);
    my $newfile2 = remove_ref_mismatches($newfile);
}
close(LIST);

sub remove_ref_mismatches
{
    my ($file) = @_;
    my $out = $file."remv_mismatch";
    my $cmd = $cortex_root."scripts/analyse_variants/vcf-hack/bin/vcfref -s $file tests/ref.fa > $out";
    my $ret =qx{$cmd};
    print "$cmd\n$ret\n";
    return $out;
}

sub parse_file
{
    my ($file) = @_;
    open(FILE, $file)||die();
    my $out = $file.".reduced";
    open(OUT, ">".$out)||die();
    while (<FILE>)
    {
	my $line = $_;
	chomp $line;
	if ($line =~ /^\#/)
	{
	    if ($line !~ /CHROM/)
	    {
		print OUT "$line\n";
	    }
	    else
	    {
		my @head = split(/\s+/, $line);
		my @new_head = @head[0..7];
		print OUT join("\t", @new_head);
	    }
	}
	else
	{
	    my @sp = split(/\t/, $line);
	    my $filter = $sp[6];
	    if ($filter eq "PASS")
	    {
		my @newsp = @sp[0..7];
		print OUT join("\t", @newsp);
		print OUT "\n";
	    }
	}
    }
    close(FILE);
    close(OUT);
    return $out;
}

