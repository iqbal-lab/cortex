#!/usr/bin/perl -w
use strict;


my $dir=shift; ## dir containing all the sv_called_in_chrom_* files
my $threeprime_anchor=shift;

if ($dir !~ /\/$/)
{
    $dir=$dir."/";
}


## Will generate fasta files, one for each of those files, all in current dir, of flanking regions

my $f=$dir."sv_called_in_chrom_";
my @chroms=("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y");

my $i;
for ($i=0; $i<24; $i++)
{
    my $cmd = "perl  ~/dev/hg/marzam/CORTEX/scripts/parsing_and_mapping_calls/generate_flank_fasta_from_trusted_sv_calls.pl ".$f.$chroms[$i]." ".$chroms[$i]." ".$threeprime_anchor;
    print "$cmd\n";

    my $ret=qx{$cmd};
    print $ret;

}

for ($i=0; $i<24; $i++)
{
    my $cmd = "perl ~/dev/hg/scripts/qnd/split_fasta_into_chunks_with_N_reads.pl chrom_".$chroms[$i]."_variant_flanks.fasta 500 ".$chroms[$i];
    print "$cmd\n";
    
    my $ret = qx{$cmd};
    print $ret;
}
