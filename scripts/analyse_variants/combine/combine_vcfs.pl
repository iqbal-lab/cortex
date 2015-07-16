#!/usr/bin/perl -w
use strict;
use File::Basename;
use File::Spec;
use Getopt::Long;
use Benchmark;



### Take a set of single-sample VCFs and combine them into one "sites" VCF and make a cortex graph of alleles and reference-intersect-bubbles

my $list = "";
my $cortex_dir = "";
my $vcftools_dir = "";
my $outdir = "";
my $outstub = "";
my $kmer = 31;
my $refname = "REF"; # eg Pf3d7_v3 or GRC38
my $ref_fasta = ""; #one fasta file for the reference genome
my $ref_binary = ""; ## cortex binary file for reference genome
my $bubble_mem_height = 14;
my $bubble_mem_width = 100;
my $run_calls_outdir = "";# root dir below which we have sample_names and then below thatbinaries/ vcfs/ etc
my $mem_height =20;
my $mem_width = 100;


&GetOptions(
    'list_vcfs:s'       => \$list,
    'cortex_dir:s'      => \$cortex_dir,
    'vcftools_dir:s'    => \$vcftools_dir,
    'outdir:s'          => \$outdir,
    'prefix:s'          => \$outstub,
    'refname:s'         => \$refname,
    'ref_fasta:s'       => \$ref_fasta,
    'ref_binary:s'       => \$ref_binary,
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



push( @INC,$cortex_dir. "scripts/analyse_variants/perl_modules/Statistics-Descriptive-2.6");
push(@INC, $scripts_dir."lib/");




## use full paths
my $ref_fa_cmd ="readlink -f $ref_fasta";
print $ref_fa_cmd."\n";
$ref_fasta=qx{$ref_fa_cmd};
chomp $ref_fasta;
my $ref_bin_cmd ="readlink -f $ref_binary"; 
print $ref_bin_cmd."\n";
$ref_binary=qx{$ref_bin_cmd};
chomp $ref_binary;

my $output_config = $outdir."config.txt";
my $fh_CONFIG;
open($fh_CONFIG, ">".$output_config)||die("Cannot open output config file $output_config\n");
print $fh_CONFIG  "cortex_dir\t$cortex_dir\n";
print $fh_CONFIG "kmer\t$kmer\n";



###############################################
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
my $outvcf2 = $outvcf1.".annot_flanks";

if ($ref_binary eq "")
{
    exit(0);
}


## annotate flanks

my $af = $scripts_dir."vcf_add_flanks.pl";
my $kp4 = $kmer+4;
print "Annotate flanks\n";
my $cmd2 = "perl $af $kp4 $outvcf1 $refname $ref_fasta > $outvcf2";
print "$cmd2\n";
my $ret2 = qx{$cmd2};
print "$ret2\n";

#print $fh_CONFIG "outvcf\t$outvcf1\n";




###############################################
##make pseudo callfile
my $pseudo_callfile = $outdir.$outstub.".pseudo_callfile";
my $ps = $combine_dir."build_pseudo_callfile_from_vcf.pl";
print "Make pseudo callfile\n";
my $cmd3 = "perl $ps $outvcf2 $pseudo_callfile";
print "$cmd3\n";
my $ret3 = qx{$cmd3};
print "$ret3\n";
print $fh_CONFIG "bubble_callfile\t$pseudo_callfile\n";


###############################################    
## make branch file, and dump binaries
print "Make binary graph file of the bubbles/sites\n";

my $branches = $pseudo_callfile.".branches.fasta";
my $branches_log = $branches.".log";
my $cmd_b1 = "perl $analyse_dir"."make_branch_fasta.pl --callfile $pseudo_callfile --kmer $kmer > $branches_log 2>&1 ";
print "$cmd_b1\n";
my $ret_b1 = qx{$cmd_b1};
print "$ret_b1\n";


##calculate max allele length
my $max_allele_len = calculate_max_allele_len($branches);
if ($max_allele_len==0)
{
    die("Error creating pseudocallfile - seems to be empty\n");
}
print $fh_CONFIG "max_allele\t$max_allele_len\n";



if (!(-e $branches))
{
    die("Failed to create $branches. See logfile $branches_log\n");
}

my $branches_list = $branches.".list";
my $bn = basename($branches);
my $cmd_b2 = "echo $bn > $branches_list";
print "$cmd_b2\n";
my $ret_b2 = qx{$cmd_b2};
print "$ret_b2\n";

my $ctx_binary1 = check_cortex_compiled_1colour($cortex_dir, $kmer);
my $bubble_graph = $branches;
my $k = "k".$kmer;
$bubble_graph =~ s/fasta/$k.ctx/;
my $bubble_graph_log = $bubble_graph.".log";
my $cmd_b3 = $ctx_binary1." --se_list $branches_list --mem_height $mem_height --mem_width $mem_width --kmer_size $kmer  --dump_binary $bubble_graph  > $bubble_graph_log 2>&1";
print "$cmd_b3\n";
my $ret_b3 = qx{$cmd_b3};
print "$ret_b3\n";
print $fh_CONFIG "bubble_graph\t$bubble_graph\n";

print "Finished building a graph just of the bubble branches/alleles. Now intersect the ref binary\n";


my $ctx_binary2 = check_cortex_compiled_2colours($cortex_dir, $kmer);
my $suffix = "intersect_bubbles";
my $ref_intersect_log = basename($ref_binary.".intersect_with_bubbles.log");
$ref_intersect_log=$outdir.$ref_intersect_log;
my $ref_col_list=get_ref_col_list($ref_binary, $outdir);
my $new_ref_binary = $ref_intersect_log;
$new_ref_binary =~ s/log/ctx/;

my $cmd_b4 = $ctx_binary2." --kmer_size $kmer --multicolour_bin $bubble_graph --mem_height $mem_height --mem_width $mem_width --colour_list $ref_col_list  --load_colours_only_where_overlap_clean_colour 0 --successively_dump_cleaned_colours $suffix  > $ref_intersect_log 2>&1";
print "$cmd_b4\n";
my $ret_b4 = qx{$cmd_b4};
print "$ret_b4\n";

print $fh_CONFIG "ref_overlap_bubble_graph\t$new_ref_binary\n";

close($fh_CONFIG);


## Now just make a filelist of sample graphs

make_sample_graph_filelist($run_calls_outdir, $kmer, $outdir);



printf("DONE. Please note this has output a config file: $output_config which is needed as an argument by scripts/calling/genotype_1sample_against_sites.pl\n");


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


sub get_ref_col_list
{
    my ($ref, $odir) = @_;
    my $ref_col_list = $odir.basename($ref.".colour_list");
    my $bn_ref_col_list = basename($ref.".colour_list");
    my $ref_list = $odir.basename($ref.".list");
    my $bn_ref_list = basename($ref.".list");
    my $c1 = "ls $ref > $ref_list";
    print "$c1\n";
    my $r1=qx{$c1};
    my $c2 = "echo $bn_ref_list  > $ref_col_list";
    print "$c2\n";
    my $r2 =qx{$c2};
    print "$r2\n";
    return $ref_col_list;
}



sub check_cortex_compiled_1colour
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
	die("Please go and compile Cortex in $cortex_dir for k=$kmer with 1 colour");
    }
}


sub check_cortex_compiled_2colours
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
	
    if (-e $dir."bin/cortex_var_".$maxk."_c2")
    {
	return $dir."bin/cortex_var_".$maxk."_c2";
    }
    else
    {
	die("Please go and compile Cortex in $cortex_dir for k=$kmer with 1 colour");
    }
}

sub calculate_max_allele_len
{
    my ($fasta) = @_;
    
    my $max = 0;
    open(FA, $fasta)||die();
    while (<FA>)
    {
	my $readline = $_;
	if ($readline=~ /^>/)
	{
	    $readline = <FA>;
	    chomp $readline;
	    my $len = length($readline);
	    if ($len>$max)
	    {
		$max=$len;
	    }

	}
    }
    close(FA);
    return $max;
    
}
