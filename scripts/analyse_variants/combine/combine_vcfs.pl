#!/usr/bin/perl -w
use strict;


### Take a set of single-sample VCFs and combine them into one "sites" VCF and make a cortex graph of alleles and reference-intersect-bubbles

my $ref_binary="";

my $list = shift;
my $cortex_dir = shift;
my $vcftools_dir = shift;
my $outdir = shift;
my $outstub = shift;
my $kmer = shift;
my $refname = shift; # eg Pf3d7_v3 or GRC38
my $ref_fasta = shift; #one fasta file for the reference genome
$ref_binary = shift; ## cortex binary file for reference genome
my $bubble_mem_height = shift;
my $bubble_mem_width = shift;
my $run_calls_outdir = shift;# root dir below which we have sample_names and then below thatbinaries/ vcfs/ etc



&GetOptions(
    'list_vcfs:s'       => \$list,
    'cortex_dir:s'      => \$cortex_dir,
    'vcftools_dir:s'    => \$vcftools_dir,
    'outdir:s'          => \$outdir,
    'prefix:s'          => \$outstub,
    'refname:s'         => \$refname,
    'ref_fasta:s'       => \$ref_fasta,
    'rootdir_for_sample_output:s'       => \$run_calls_outdir,
    'kmer:i'            =>\$kmer,
    'mem_height:i'      =>\$mem_height,
    'mem_width:i'       =>\$mem_width,
);




#make sure there is a slash on the end
if ($cortex_dir !~ /\/$/)
{
    $cortex_dir=$cortex_dir.'/';
}
if ($vcftools_dir !~ /\/$/)
{
    $vcftools_dir=$vcftools_dir.'/';
}
if ($run_calls_outdir !~ /\/$/)
{
    $run_calls_outdir=$run_calls_outdir.'/';
}
if ($outdir !~ /\/$/)
{
    $outdir=$outdir.'/';
}
my $scripts_dir = $cortex_dir."scripts/analyse_variants/bioinf-perl/vcf_scripts/";
my $analyse_dir = $cortex_dir."scripts/analyse_variants/";
my $combine_dir = $cortex_dir."scripts/analyse_variants/combine/";

## use full paths
my $ref_fa_cmd ="readlink -f $ref_fasta";
print $ref_fa_cmd."\n";
$ref_fasta=qx{$ref_fa_cmd};
chomp $ref_fasta;
my $ref_bin_cmd ="readlink -f $ref_binary"; 
print $ref_bin_cmd."\n";
$ref_binary=qx{$ref_bin_cmd};
chomp $ref_binary;


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

my $cmd1 = "head -n 1000 $r1 | grep \"#\" | grep -v CHROM > $tmpvcf";
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

if ($ref_binary eq "")
{
    exit(0);
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

## Now just make a filelist of sample graphs

make_sample_graph_filelist($run_calls_outdir, $kmer, $outdir);



sub make_sample_graph_filelist
{
    my ($rcdir, $k, $odir) = @_;
    
    my $reg = $rcdir."\*/binaries/cleaned/k".$k.'/'."\*ctx";

    my $listing_cmd = "ls $reg";
    my $listing = qx{$listing_cmd};
    my @sp = split(/\n/, $listing);
    my $outfile = $odir."list_sample_ids_and_graphs";
    open(OUT, ">".$outfile)||die("Cannot open $outfile\n");
    foreach my $f (@sp)
    {
	if ($f=~ /\/([^\/]+)\.kmer/)
	{
	    my $id = $1;
	    print OUT "$id\t$f\n";
	}
    }
    close(OUT);
}
