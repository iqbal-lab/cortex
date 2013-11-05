#!/usr/bin/perl -w
use strict;

my $cortex_root = "/home/zam/dev/hg/bitbucket/CORTEX_mainline/";  ## set this!
my $list = shift; ## list of RAW vcfs, one per population
my $ref = shift; ## reference fasta
my $refname = shift;;
my $kmer_size = 31;
open(LIST, $list)||die();
while (<LIST>)
{
    my $lyne = $_;
    chomp $lyne;
    my $newfile = parse_file($lyne, $ref);
    my $newfile2 = remove_ref_mismatches($newfile);

}
close(LIST);

sub remove_ref_mismatches
{
    my ($file) = @_;
    my $out = $file.".remv_mismatch";
    my $cmd = $cortex_root."scripts/analyse_variants/vcf-hack/bin/vcfref -s $file $ref > $out";
    my $ret =qx{$cmd};
    print "$cmd\n$ret\n";
    return $out;
}

sub parse_file
{
    my ($file, $ref) = @_;
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
		print OUT "\tFORMAT\tDUMMY_SAMPLE\n";
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
		print OUT "\tGT\t0/0";#dummy
		print OUT "\n";
	    }
	}
    }
    close(FILE);
    close(OUT);

    my $flank_len = $kmer_size *2 +1;
    ##now add flanks
    my $out2 = $out.".annot_flanks";
    my $fcmd = $cortex_root."scripts/analyse_variants/bioinf-perl/vcf_scripts/vcf_add_flanks.pl $flank_len $out $refname $ref > $out2";
    my $fret = qx{$fcmd};
    print "$fcmd\n$fret\n";
    return $out2;
}

