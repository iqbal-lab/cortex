#!/usr/bin/perl -w
use strict;

use File::Basename;
use File::Spec;
use Getopt::Long;
use Benchmark;
use Cwd    qw( abs_path );
use FindBin qw($Bin);


my $all_samples_index = "";
my $index_dir="";
my $kmer = 31;
my $qthresh = 10;
my $bc="yes";
my $pd = "no";
my $mem_height = 16;
my $mem_width=100;
my $vcftools_dir = "";
my $stampy_bin="";
my $stampy_hash="";
my $list_ref="";
my $refbindir="";
my $num=-1;
my $cortex_dir = abs_path($0);
$cortex_dir =~ s/scripts\/calling\/par.pl//;

my $outdir="";

&GetOptions(
    ##mandatory args
    'num:i'                       =>\$num,
    'index:s'                     =>\$all_samples_index,
    'list_ref:s'                  =>\$list_ref,
    'refbindir:s'                 =>\$refbindir,
    'index_dir:s'                 =>\$index_dir,
    'vcftools_dir:s'              =>\$vcftools_dir,
    'out_dir:s'                   =>\$outdir,
    'stampy_bin:s'                =>\$stampy_bin,
    'stampy_hash:s'               =>\$stampy_hash,
    'bc:s'                        =>\$bc,
    'kmer:i'                      =>\$kmer,
    'mem_height:i'                =>\$mem_height,
    'mem_width:i'                 =>\$mem_width,
    'qthresh:i'                   =>\$qthresh,
    );

check_args($num, $index_dir);

if ($index_dir !~ /\/$/)
{
    $index_dir = $index_dir.'/';
}
if ($outdir !~ /\/$/)
{
    $outdir = $outdir.'/';
}
my $index = $index_dir."index_".$num;
my $c1 = "head -n $num $all_samples_index  | tail -n 1 > $index";
qx{$c1};
open(F, $index)||die();
my $line= <F>;
chomp $line;
my @sp = split(/\t/, $line);
my $sample = $sp[0];


my $log = $outdir."log_bc.".$sample;

my $cmd ="perl $cortex_dir"."scripts/calling/run_calls.pl --fastaq_index $index --first_kmer $kmer --auto_cleaning yes --bc $bc --pd $pd --outdir $outdir --ploidy 1  --mem_height $mem_height --mem_width $mem_width --qthresh $qthresh --vcftools_dir $vcftools_dir  --do_union no --logfile $log --workflow joint_par --ref CoordinatesAndInCalling --list_ref_fasta $list_ref --refbindir $refbindir --stampy_bin $stampy_bin --stampy_hash $stampy_hash --outvcf $sample ";
qx{$cmd};



sub check_args
{
    my ($n, $i_dir) = @_;
    if ($n==-1)
    {
	die("--num is a mandatory argument, specifies which sample from the INDEX file to run\n");
    }
if (!(-d $i_dir))
	{
	my $c1 = "mkdir $i_dir";
	qx{$c1};
	}	
}
