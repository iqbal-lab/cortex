#!/usr/bin/perl -w
use strict;

use File::Basename;
use File::Spec;
use Getopt::Long;
use Benchmark;

my $invcf = "";
my $stub = "";

##the following are set by the command-line arguments
my $cortex_dir = "";
my $mem_height = 16;
my $mem_width = 100;
my $kmer= 31;
my $list_sample_graphs="";
my $dir="";
my $max_allele=100000; ##longest allele which will be genotyped
my $g_size = 0; #genome size
my $num=0;
&GetOptions(
    ##mandatory args
    'cortex_dir:s'                        =>\$cortex_dir,
    'invcf:s'                             =>\$invcf,
    'dir:s'                               =>\$dir,
    'list_sample_graphs:s'                =>\$list_sample_graphs,# sample id \t path to binary
    'genome_size:i'                       =>\$g_size,
    'mem_height:i'                        =>\$mem_height,
    'mem_width:i'                         =>\$mem_width,
    'max_allele:i'                        =>\$max_allele,
    'num:i'                               =>\$num,
    'stub:s'                              =>\$stub,
    );


my $dircmd = "readlink -f $dir";
$dir = qx{$dircmd};
chomp $dir;

my $ctxcmd = "readlink -f $cortex_dir";
$cortex_dir = qx{$ctxcmd};
chomp $cortex_dir;

if ($cortex_dir !~ /\/$/)
{
    $cortex_dir=$cortex_dir.'/';
}
if ($dir !~ /\/$/)
{
    $dir=$dir.'/';
}

if (check_num($num, $list_sample_graphs)==-1)
{
    ## do nothing
}
else
{
    my ($sample_id, $file)  = get_sample_and_file($num, $list_sample_graphs);
    my $sample_outdir = $dir.$sample_id.'/';
    my $mk = "mkdir $sample_outdir";
    if (!(-d $sample_outdir))
    {
	qx{$mk};
    }
    my $ref_ov = "ls $dir"."\*"." | grep \"list_intersect_bubbles.ctx\" | grep -v filelist  ";
    print $ref_ov."\n";
    my $ret_ref_ov= qx{$ref_ov};
    chomp $ret_ref_ov;
    my $ov_log = $sample_outdir.$sample_id.".overlap_log";
    my $gt_log = $sample_outdir.$sample_id.".per_sample_gt_log";
    my $cmd = "perl $cortex_dir"."scripts/calling/genotype_1sample_against_sites.pl --invcf $invcf --bubble_graph $dir".$stub.".pseudo_callfile.branches.k".$kmer.".ctx --bubble_callfile $dir".$stub.".pseudo_callfile --outdir $sample_outdir --sample $sample_id --ref_overlap_bubble_graph $ret_ref_ov --genome_size $g_size --sample_graph $file  --overlap_log  $ov_log  --mem_height $mem_height --mem_width $mem_width  --cortex_dir $cortex_dir  --max_allele $max_allele > $gt_log 2>&1";
    print "$cmd\n";
    my $ret = qx{$cmd};
    print "$ret\n";
}

sub check_num
{
    my ($n, $list) = @_;
    open(LIST, $list)||die("Cannot open $list");
    my $ct=0;
    while (<LIST>)
    {
	my $line = $_;
	$ct++;
    }
    close(LIST);
    if ( ($n<0) || ($n>$ct-1) )
    {
	return -1;
    }
    else
    {
	return 0;
    }
}


sub get_sample_and_file
{
    my ($n, $list) = @_;
    open(LIST, $list)||die("Cannot open $list");
    my $ct=0;
    while (<LIST>)
    {
	my $line = $_;
	chomp $line;
	if ($ct==$n)
	{
	    my @sp = split(/\t/, $line);
	    my $id = $sp[0];
	    my $ctx = $sp[1];
	    close(LIST);
	    return ($id,  $ctx);
	}
	else
	{
	    $ct++;
	}
    }
    close(LIST);
    
    die("Coding error - should never get to this line in run_sample.pl\n");
    return ("","");
}
