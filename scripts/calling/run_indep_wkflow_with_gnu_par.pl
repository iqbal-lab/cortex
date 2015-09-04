#!/usr/bin/perl -w
use strict;


use File::Basename;
use File::Spec;
use Getopt::Long;
use Benchmark;
use Cwd    qw( abs_path );
use FindBin qw($Bin);

my $cortex_dir;

BEGIN
{
    $cortex_dir = abs_path($0);
    $cortex_dir =~ s/scripts\/calling\/run_indep_wkflow_with_gnu_par.pl//;
    push ( @INC, $cortex_dir."scripts/calling/");
}

use BasicUtils qw ( check_cortex_runnable add_slash is_fastq count_bases_in_fasta create_dir_if_does_not_exist );
$cortex_dir = add_slash($cortex_dir);


## Script to prepare reference graph binary, Stampy hash of reference etc

my $vcftools_dir = "";
my $stampy_bin="";
my $all_samples_index="";
my $ref_fa="";
my $refdir = "";
my $outdir="";
my $kmer="";
my $num_procs="";
my $mem_height="";
my $prefix = "";
my $mem_width="";
&GetOptions(
    'index:s'                     =>\$all_samples_index,
    'ref_fa:s'                    =>\$ref_fa,
    'prefix:s'                    =>\$prefix,
    'dir_for_ref_objects:s'       =>\$refdir, 
    'vcftools_dir:s'              =>\$vcftools_dir,
    'outdir:s'                    =>\$outdir,
    'stampy_bin:s'               =>\$stampy_bin,
    'kmer:i'                      =>\$kmer,
    'procs:i'                     =>\$num_procs,
    'mem_height:i'                     =>\$mem_height,
    'mem_width:i'                     =>\$mem_width,
    );




### Prepare

my $prep = "perl $cortex_dir"."scripts/calling/prepare.pl --index $all_samples_index --ref_fa $ref_fa";
$prep .= " --dir_for_ref_objects $refdir --vcftools_dir $vcftools_dir --outdir $outdir ";
$prep .= " --stampy_bin $stampy_bin --kmer $kmer";
#my $ret_prep = qx{$prep};

## Build

my $num_cmd = "wc -l $all_samples_index";
my $num_samples = qx{$num_cmd};
if ($num_samples =~ /^(\d+)/)
{
    $num_samples = $1;
}
else
{
    die("index file $all_samples_index is malformed/missing. \"wc -l\" fails to count the number of lines in it\n");
}

my $build = "parallel --gnu -j $num_procs perl $cortex_dir"."scripts/calling/build_samples_parallel.pl --num {} ";
$build .= " --index $all_samples_index --outdir $outdir  --kmer $kmer ";
if ( ($mem_height ne "") && ($mem_width ne "")) 
{
    $build .= " --mem_height $mem_height --mem_width $mem_width ";
}
$build .= " ::: {1..".$num_samples."}";

#my $ret_build = qx{$build};


## Combine

my $combine = "perl $cortex_dir"."scripts/analyse_variants/combine/combine_vcfs.pl ";
$combine .= " --prefix $prefix --outdir $outdir --intersect_ref ";
#my $ret_combine = qx{$combine};
#print $ret_combine;

## Genotype (in parallel)

my $gt = "cat $outdir"."combine/list_args_for_final_step | ";
$gt .= "parallel --gnu -j $num_procs --colsep \'\\t\' ";
$gt .= $cortex_dir."scripts/calling/gt_1sample.pl ";
$gt .= "--config ".$outdir."combine/config.txt ";
$gt .= "--invcf ".$outdir."/combine/".$prefix.".sites_vcf ";
$gt .= " --sample {1} --outdir {2} --sample_graph {3} ";
print "$gt\n";
my $ret_gt = qx{$gt};
print $ret_gt;
