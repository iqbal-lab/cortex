#!/usr/bin/perl -w
use strict;

use File::Basename;
use File::Spec;
use Getopt::Long;
use Benchmark;
use Cwd    qw( abs_path );
use FindBin qw($Bin);


my $all_samples_index = "";
my $kmer = 31;
my $qthresh = 10;
my $bc="yes";
my $pd = "no";
my $mem_height = 0;
my $mem_width=0;
my $genome_size=3000000;
my $vcftools_dir = "";
my $stampy_bin="";
my $stampy_hash="";
my $list_ref="";
my $refbindir="";
my $num=-1;
my $cortex_dir = abs_path($0);
my $genome_size=3000000;
$cortex_dir =~ s/scripts\/calling\/par.pl//;

my $outdir="";

&GetOptions(
    ##mandatory args
    'num:i'                       =>\$num,
    'index:s'                     =>\$all_samples_index,
    'list_ref:s'                  =>\$list_ref,
    'refbindir:s'                 =>\$refbindir,
#    'index_dir:s'                 =>\$index_dir,
    'vcftools_dir:s'              =>\$vcftools_dir,
    'out_dir:s'                   =>\$outdir,
    'stampy_bin:s'                =>\$stampy_bin,
    'stampy_hash:s'               =>\$stampy_hash,
    'bc:s'                        =>\$bc,
    'pd:s'                        =>\$pd,
    'kmer:i'                      =>\$kmer,
    'mem_height:i'                =>\$mem_height,
    'mem_width:i'                 =>\$mem_width,
    'genome_size:i'                 =>\$genome_size,
    'qthresh:i'                   =>\$qthresh,
    );


if ($index_dir !~ /\/$/)
{
    $index_dir = $index_dir.'/';
}
if ($outdir !~ /\/$/)
{
    $outdir = $outdir.'/';
}

my $index_dir =$outdir."indexes/";
my $config = $outdir."config.par.txt";
my $ref_binary;
($mem_height, $mem_width, $ref_binary) = check_args($num, $index_dir, $mem_height, $mem_width, $genome_size, $refbindir, $kmer);
open(CONFIG, ">".$config)||die("Unable to create $config\n");
print CONFIG "vcftools_dir\t$vcftools_dir\n";
print CONFIG "kmer\t$kmer\n";
print CONFIG "ref_binary\t$ref_binary\n";
print CONFIG "genome_size\t$genome_size\n";
print CONFIG "mem_height\t$mem_height\n";
print CONFIG "mem_width\t$mem_width\n";
close(CONFIG);


if ($mem_height==0)

my $index = $index_dir."index_".$num;
my $c1 = "head -n $num $all_samples_index  | tail -n 1 > $index";
qx{$c1};
open(F, $index)||die();
my $line= <F>;
chomp $line;
my @sp = split(/\t/, $line);
my $sample = $sp[0];
my $odir = $outdir.$sample.'/';
my $c2 = "mkdir -p $odir";
qx{$c2};
my $log = $odir."log_bc.".$sample;
print "$log\n";

my $cmd ="perl $cortex_dir"."scripts/calling/run_calls.pl --fastaq_index $index --first_kmer $kmer --auto_cleaning yes --bc $bc --pd $pd --outdir $odir --ploidy 1 --genome_size $genome_size --mem_height $mem_height --mem_width $mem_width --qthresh $qthresh --vcftools_dir $vcftools_dir  --do_union yes --logfile $log --workflow independent --ref CoordinatesAndInCalling --list_ref_fasta $list_ref --refbindir $refbindir --stampy_bin $stampy_bin --stampy_hash $stampy_hash --outvcf $sample ";
qx{$cmd};



sub check_args
{
    my ($n, $i_dir, $mh, $mw, $g, $rbindir, $km) = @_;
    if ($km %%2==0)
    {
	die("Kmer must be an odd number\n");
    }
    elsif ($km<10) 
    {
	die("I'm just going to stop this being run for k<10, can't believe it will be useful - typo?\n");
    }
    if ($n==-1)
    {
	die("--num is a mandatory argument, specifies which sample from the INDEX file to run\n");
    }
    if (!(-d $i_dir))
    {
	my $c1 = "mkdir $i_dir";
	qx{$c1};
    }	
    if ($g==0)
    {
	die("You must enter a genome size (in bp) using --genome_size. An estimate is fine. This will be used to get memory use parameters mem ehight/width (if you do not enter them) and later on will be used for likelihood calculations\n");
    }
    if ( ($mh==0) ||  ($mw==0))
    {
	$mw=100;
	$mh = int(log(2*$g+0.5)/log(2));
    }

    if (!(-d $rbindir))
    {
	die("The specified directory $rbindir does not exist\n");
    }

    my @files = glob($rbindir."/*k".$km.".ctx");
    if (scalar @files==0)
    {
	die("The reference binary directory $rbindir does not contain any files with names ending k$kmer".".ctx; \neither you have not followed the naing convention we are asking for, \nor there is no Cortex binary graph file of the reference genome in that directory\n");
    }
    elsif (scalar @files>1)
    {
	die("There is more than one file in the reference binary directory $rbindir with name ending k$kmer".".ctx, so this script can't work out which is the reference binary\n");
    }
    return ($mh, $mw, $files[0]);
}
