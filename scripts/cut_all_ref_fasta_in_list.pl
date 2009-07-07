#!/usr/bin/perl -w
use strict;

my $filelist=shift;
my $desired_read_length=500;

open(FILE,$filelist)||die("Cannot open $filelist");

while(<FILE>)
{
    my $line=$_;
    chomp $line;

    my $cmd = "perl /homes/zi/dev/hg/marzam/cortex/src/scripts/cut_reference_fasta_into_manageable_reads.pl $line $desired_read_length";
    print "$cmd\n";
    my $ret=qx{$cmd};
    print $ret;
}
