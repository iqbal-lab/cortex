#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $sample = "";
my $root  = "/ddn/projects3/mcvean_res/zam/phase2_cortex/dummyrun/samples/";
my $cortex_dir = "/ddn/projects3/mcvean_res/zam/phase2_cortex/dummyrun/code/cortex_release_isaac_branch_bitbucket_sequoia/";
my $kmer= 31;
my $outdir="";

&GetOptions(
        'sample|s:s'                               => \$sample,
        'root|r:s'                               => \$root,
);



if ($sample eq "")
{
    die("You must use --sample and enter a sample name, which we assume is also a directory name inside $root");
}
else
{
    $outdir = $root.$sample.'/';
}

## Assumes fastq are downloaded (they do NOT need to be unzipped any more), 
##  and filelists are made and valid (contain only files that exists on the filelsystem)

my $se_list = $outdir.$sample."_se";
my $pe1_list = $outdir.$sample."_pe1";
my $pe2_list = $outdir.$sample."_pe2";


#### 0. Check args
if ($cortex_dir !~ /\/$/)
{
    $cortex_dir=$cortex_dir.'/';
}

my $ctx_binary = check_cortex_compiled($cortex_dir, $kmer);

if (!(-d $outdir))
{
    die("Outdir does not exist - this script expects each sample to have a directory named after it. THis dir contains the fastq and will also contain the binaries we build and everything else related to that sample\n");
}
if (!(-e $se_list))
{
    die("Abort mission. The filelist $se_list does not exist - have you reun the rpevious script??");
}

if (!(-e $pe1_list))
{
   die("Abort mission. The filelist $pe1_list does not exist - have you reun the rpevious script??");
}
if (!(-e $pe2_list))
{
    die("Abort mission. The filelist $pe2_list does not exist - have you reun the rpevious script??");
}


my $binname = $outdir.$sample.".uncleaned.q10.k31.ctx";
my $log     = $binname.".log";

my $cmd = $ctx_binary." --kmer_size $kmer --sample_id $sample --mem_height 25 --mem_width 150 --se_list $se_list --pe_list $pe1_list,$pe2_list  --quality_score_threshold 10  --dump_binary $binname --remove_pcr_duplicates > $log 2>&1";
print "$cmd\n\n";
my $ret = qx{$cmd};
print "$ret\n";
print "Finished building uncleaned graph pf $sample\n";



sub check_cortex_compiled
{
    my ($dir, $k) = @_;

    my $maxk=31;
    if (($k>32) && ($k<64) )
    {
	$maxk=63;
    }
    elsif ($k<32)
    {
	#maxk is already 31
    }
    else
    {
	die("This script expects you to use a k<64");
    }
	
    if (-e $dir."bin/cortex_var_".$maxk."_c1")
    {
	return $dir."bin/cortex_var_".$maxk."_c1";
    }
    else
    {
	die("Please go and compile Cortex in $cortex_dir for k=$kmer");
    }
}

