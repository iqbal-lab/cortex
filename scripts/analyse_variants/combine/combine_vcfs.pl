#!/usr/bin/perl -w
use strict;


### Take a set of single-sample VCFs and combine them into one "sites" VCF

my $list = shift;
my $cortex_dir = shift;
my $vcftools_dir = shift;
my $scripts_dir = shift;
my $analyse_dir = shift;
my $combine_dir = shift;
my $outdir = shift;
my $outstub = shift;
my $kmer = shift;
my $refname = shift; # eg Pf3d7_v3
my $ref_fasta = shift;
my $ref_binary = shift;
my $bubble_mem_height = shift;
my $bubble_mem_width = shift;

#make sure there is a slash on the end
if ($vcftools_dir !~ /\/$/)
{
    $vcftools_dir=$vcftools_dir.'/';
}
if ($scripts_dir !~ /\/$/)
{
    $scripts_dir=$scripts_dir.'/';
}
if ($analyse_dir !~ /\/$/)
{
    $analyse_dir=$analyse_dir.'/';
}
if ($outdir !~ /\/$/)
{
    $outdir=$outdir.'/';
}
if ($combine_dir !~ /\/$/)
{
    $combine_dir=$combine_dir.'/';
}


## 1. Cat all the VCFs, sort, remove duplicates


my $outvcf1 = $outdir.$outstub.".sites_vcf";
my $tmpvcf = $outvcf1.".tmp_delete_me";
my $header = "\#\#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDUMMY";
my $rd = $scripts_dir."vcf_remove_dupes.pl";
my $ro = $scripts_dir."vcf_remove_overlaps.pl";
print "Make sites VCF\n";

my $c1 = "head -n1 $list";
my $r1 = qx{$c1};
chomp $r1;

my $cmd1 = "head -n 300 $r1 | grep \"#\" | grep -v CHROM > $tmpvcf";
print "$cmd1\n";
my $rcmd1 = qx{$cmd1};
print "$rcmd1\n";

open(H, ">>".$tmpvcf);
print H $header."\n";
close(H);

my $cmd111="cat $list | xargs cat | grep -v \"\#\" | grep PASS | $vcftools_dir/perl/vcf-sort >> $tmpvcf";
print "$cmd111\n";
my $rcmd111 = qx{$cmd111};
print "$rcmd111\n";

my $cmd1111 = "perl $rd  --take_first --pass $tmpvcf > $outvcf1"; #| $ro --pass --filter_txt OVERLAPPING_SITE | grep -v \"\" >  $outvcf1";
print "$cmd1111\n";
my $rcmd1111=qx{$cmd1111};
print "$rcmd1111\n";

if (!(-e $outvcf1))
{
    die("ERROR - Failed to build sites vcf\n");
}

## annotate flanks
my $outvcf2 = $outvcf1.".annot_flanks";
my $af = $scripts_dir."vcf_add_flanks.pl";
my $kp4 = $kmer+4;
print "Annotate flanks\n";
my $cmd2 = "perl $af $kp4 $outvcf1 $refname $ref_fasta > $outvcf2";
print "$cmd2\n";
my $ret2 = qx{$cmd2};
print "$ret2\n";

##make pseudo callfile
my $pseudo_callfile = $outdir.$outstub.".pseudo_callfile";
my $ps = $combine_dir."build_pseudo_callfile_from_vcf.pl";
print "Make pseudo callfile\n";
my $cmd3 = "perl $ps $outvcf2 $pseudo_callfile";
print "$cmd3\n";
my $ret3 = qx{$cmd3};
print "$ret3\n";

## make branch file, and dump binaries
my $dmp = $combine_dir."dump_bubble_branch_graph.pl ";
print "Make binary graph file of the bubbles/sites\n";
my $cmd4 = "perl $dmp --mem_height $bubble_mem_height --mem_width $bubble_mem_width  --ref $ref_binary  --bubble_calls $pseudo_callfile --kmer $kmer --analyse_dir $analyse_dir  --cortex_dir $cortex_dir --out_dir $outdir ";
print "$cmd4\n";
my $ret4 =qx{$cmd4};
print "$ret4\n";


