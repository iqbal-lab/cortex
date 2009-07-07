#!/usr/bin/perl -w
use strict;

my $ref_fasta = shift;
my $read_len_desired =shift;

open(REF, $ref_fasta) || die("Cannot open $ref_fasta");
open(OUT, "> ".$ref_fasta.".short_reads")||die("Cannot open output file");
my $read_id = <REF>;
my $read_id_root;

##must start with >
if ($read_id !~ /^>/)
{
    die("Read id line $read_id does not start with >");
}
elsif ($read_id =~ /^(\S+)/)
{
    $read_id_root=$1;
    chomp $read_id_root;
}

my $counter=0;
my $next_read;
my $length;

while (($length = read REF, $next_read, $read_len_desired) != 0) {
    chomp $next_read;
    print OUT "$read_id_root $counter $length bases\n$next_read\n";
    $counter++;
}

close(REF);

open(REF, $ref_fasta) || die("Cannot open $ref_fasta");
<REF>;

##shift along by half a read from the start. ignore this.
$length=read REF, $next_read, $read_len_desired/2;

while (($length = read REF, $next_read, $read_len_desired) != 0) {
    chomp $next_read;
    print OUT "$read_id_root $counter $length bases\n$next_read\n";
    $counter++;
}

close(REF);


