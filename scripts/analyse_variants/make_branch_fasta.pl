#!/usr/bin/perl -w
use strict;

use Getopt::Long;


my $callfile="";
my $outfile="";
my $help='';
my $kmer = -1;

&GetOptions(
	'callfile|f:s'                         => \$callfile,
        'kmer|k:i'                             => \$kmer,
	'help'                                 => \$help,
);

if ($help)
{
    print "Script for making a fasta just of the branches of a set of calls, including the k-1 bases at the end of the 5p flank\n";
    print "Usage: perl make_branch_fasta.pl --callfile <filename of output of bubble caller or PD caller, does no matter if you used --print_colour_coverages>\n";
    print "                                 --kmer <kmer size used to call variants>\n";
    print "Will create a file with name = callfile.branches.fasta.\nWill abort if such a file already exists\n";
    exit(0);
}

my $out100 = $callfile.".branches.fasta";

check_args($callfile, $out100, $kmer);


open(FILE, $callfile)||die();
open(OUT, ">".$out100)||die();

my $k_plus_one = $kmer+1;

while (<FILE>)
{
    my $line = $_;
    if ($line =~ /var_\d+_5p_flank/)
    {
	$line =<FILE>;## 5p flank sequence
	my $last_k = substr($line, -$k_plus_one);
	chomp $last_k;

	$line =<FILE>;
	if ($line =~ /^(>\S+)/)
	{
	    $line = $1."  with previous k+1 bases prepended\n";
	}
	print OUT $line;#br1 readid
	$line =<FILE>;
	print OUT $last_k.$line;#br1 seq
	$line =<FILE>;
	if ($line =~ /^(>\S+)/)
	{
	    $line = $1."\n";
	}
	print OUT $line;#br2 readid
	$line =<FILE>;
	print OUT $last_k.$line;#br2 seq
	$line =<FILE>;
	$line =<FILE>;
    }
}
close(FILE);
close(OUT);

print "Done. Finished creating a branches file\n";

sub check_args
{
    my ($c, $f, $k) = @_;
    
    if ($c eq "")
    {
	die("You  must specify --callfile and give a file as argument\n");
    }
    if ($k==-1)
    {
	die("You must specify a kmer with --kmer\n");
    }
    
    if (!(-e $c))
    {
	die("The callfile you specify, $c, does not seem to exist\n");
    }
    if (-e($f))
    {
	die("The output filename this script wants to create, $f, already exists, and I do not want to overwrite it. Delete it and then rerun this script\n");
    }

    return;
}
