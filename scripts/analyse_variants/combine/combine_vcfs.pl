#!/usr/bin/perl -w
use strict;
use File::Basename;
use File::Spec;
use Getopt::Long;
use Benchmark;
use Cwd    qw( abs_path );
use FindBin qw($Bin);



### Take a set of single-sample VCFs and combine them into one "sites" VCF and make a cortex graph of alleles and reference-intersect-bubbles

my $list = "";

my $cortex_dir = abs_path($0);
$cortex_dir =~ s/scripts\/analyse_variants\/combine\/combine_vcfs.pl//;

my $outdir = "";
my $outstub = "";

my $refname = "REF"; # eg Pf3d7_v3 or GRC38
my $ref_fasta = ""; #one fasta file for the reference genome
my $ref_binary = ""; ## cortex binary file for reference genome
my $bubble_mem_height = 14;
my $bubble_mem_width = 100;
my $mem_height =20;
my $mem_width = 100;
my $vcftools_dir = "";
my $kmer = 31;

&GetOptions(
    'list_vcfs:s'       => \$list,
    #'cortex_dir:s'      => \$cortex_dir,
    'vcftools_dir:s'    => \$vcftools_dir,
    'outdir:s'          => \$outdir,
    'prefix:s'          => \$outstub,
    'ref_fasta:s'       => \$ref_fasta,
    'ref_binary:s'       => \$ref_binary,
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

if ($outdir !~ /\/$/)
{
    $outdir = $outdir.'/';
}

my $comb_dir = $outdir."combine/";
my $c = "mkdir $comb_dir";
if (!(-d $comb_dir))
{
    qx{$c};
}

my $genome_size=0;
if (-e $outdir."config.par.txt")
{
    ($kmer, $vcftools_dir, $ref_binary, $genome_size, $mem_height, $mem_width) =
	get_args_from_par_config($outdir."config.par.txt");
}
if (-e  $outdir."config.prep.txt")
{
    $ref_fasta =get_args_from_prep_config($outdir."config.prep.txt"); 
}

my $scripts_dir = $cortex_dir."scripts/analyse_variants/bioinf-perl/vcf_scripts/";
my $analyse_dir = $cortex_dir."scripts/analyse_variants/";
my $combine_dir = $cortex_dir."scripts/analyse_variants/combine/";


if ($list eq "")
{
    if (-e $outdir."list_all_raw_vcfs")
    {
	$list = $outdir."list_all_raw_vcfs";
    }
    else
    {
	## go to the $outdir and make the list
	my $list = $outdir."list_all_raw_vcfs";
	my $cmd = "ls $outdir"."*/vcfs/*wk*raw* > $list";
	qx{$cmd};
	if (!-e ($list))
	{
	    die("Unable to create $list\n");
	}
    }
}

sub get_args_from_prep_config
{
    my ($config_file) = @_;
    open(CONFIG, $config_file)||die("Cannot open pre $config_file\n");
    my ($rf)="";
    while (<CONFIG>)
    {
	my $line = $_;
	chomp $line;
	my @sp = split(/\t/, $line);
	if ($sp[0] eq "ref_fa")
	{
	    $rf = $sp[1];
	}
    }
    close(CONFIG);

    if ($rf eq "")
    {
	die("Missing tag in pre config file\n");
    }
    return $rf;
    
}


sub get_args_from_par_config
{
    my ($config_file) = @_;
    open(CONFIG, $config_file)||die("Cannot open $config_file\n");
    my ($mh, $mw, $g, $rb, $k, $vcft)=(-1,-1,-1,-1,-1,-1);
    while (<CONFIG>)
    {
	my $line = $_;
	chomp $line;
	my @sp = split(/\t/, $line);
	if ($sp[0] eq "vcftools_dir")
	{
	    $vcft = $sp[1];
	}
	elsif ($sp[0] eq "kmer")
	{
	    $k = $sp[1];
	}
	elsif ($sp[0] eq "mem_height")
	{
	    $mh = $sp[1];
	}
	elsif ($sp[0] eq "mem_width")
	{
	    $mw = $sp[1];
	}
	elsif ($sp[0] eq "genome_size")
	{
	    $g = $sp[1];
	}
	elsif ($sp[0] eq "ref_binary")
	{
	    $rb = $sp[1];
	}
	else
	{
	    print ("unexpected tag in $config_file\n");
	    die($sp[0]);
	}
    }
    close(CONFIG);
    if ( ($k eq "-1") || ($mh eq "-1") || ($mw eq "-1") 
	 || ($rb eq "-1") || ($vcft eq "-1") || ($g eq "-1") )
    {
	die("Missing tag in $config_file\n");
    }
    return ($k, $vcft, $rb, $g, $mh, $mw);
	
}
sub get_stat
{
    return abs_path("$Bin/../perl_modules/Statistics-Descriptive-2.6/");
}

sub get_vcflib
{
    return abs_path("$Bin/../bioinf-perl/lib/");
}


BEGIN
{
    push @INC, get_stat();
    push @INC, get_vcflib();
}

use Descriptive;
use VCFFile;
my $libdir      = get_vcflib();
print "ZAMZAM libdir is $libdir\n";

## use full paths
#my $ref_fa_cmd ="readlink -f $ref_fasta";
#print $ref_fa_cmd."\n";
#$ref_fasta=qx{$ref_fa_cmd};
#chomp $ref_fasta;
#my $ref_bin_cmd ="readlink -f $ref_binary"; 
#print $ref_bin_cmd."\n";
#$ref_binary=qx{$ref_bin_cmd};
#chomp $ref_binary;

if (!(-d $comb_dir))
{
    my $c = "mkdir -p $comb_dir";
    qx{$c};
}
my $output_config = $comb_dir."config.txt";
my $fh_CONFIG;
open($fh_CONFIG, ">".$output_config)||die("Cannot open output config file $output_config\n");
print $fh_CONFIG "kmer\t$kmer\n";



###############################################
## 1. Cat all the VCFs, sort, remove duplicates


my $outvcf1 = $comb_dir.$outstub.".sites_vcf";
my $tmpvcf = $outvcf1.".tmp_delete_me";
my $header = "\#\#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDUMMY";
my $rd = $scripts_dir."vcf_remove_dupes.pl";
my $ro = $scripts_dir."vcf_remove_overlaps.pl";
print "Make sites VCF using list $list\n";

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

my $cmd111="cat $list | xargs cat | grep -v \"\#\" | grep PASS | $vcftools_dir"."/perl/vcf-sort >> $tmpvcf";
print "$cmd111\n";
my $rcmd111 = qx{$cmd111};
print "$rcmd111\n";

my $cmd1111 = "perl  -I $libdir $rd  --take_first --pass $tmpvcf > $outvcf1"; #| ≈ß$ro --pass --filter_txt OVERLAPPING_SITE | grep -v \"\" >  $outvcf1";
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
my $cmd2 = "perl  -I $libdir $af $kp4 $outvcf1 $refname $ref_fasta > $outvcf2";
print "$cmd2\n";
my $ret2 = qx{$cmd2};
print "$ret2\n";



###############################################
##make pseudo callfile
my $pseudo_callfile = $comb_dir.$outstub.".pseudo_callfile";
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
my $ref_intersect_log = basename($ref_binary.".intersect_bubbles.log");
$ref_intersect_log=$comb_dir.$ref_intersect_log;
my ($reflist, $ref_col_list) =get_ref_col_list($ref_binary, $comb_dir);
my $new_ref_binary = $reflist."_intersect_bubbles.ctx";


my $cmd_b4 = $ctx_binary2." --kmer_size $kmer --multicolour_bin $bubble_graph --mem_height $mem_height --mem_width $mem_width --colour_list $ref_col_list  --load_colours_only_where_overlap_clean_colour 0 --successively_dump_cleaned_colours $suffix  > $ref_intersect_log 2>&1";
print "$cmd_b4\n";
my $ret_b4 = qx{$cmd_b4};
print "$ret_b4\n";

print $fh_CONFIG "ref_overlap_bubble_graph\t$new_ref_binary\n";

close($fh_CONFIG);


## Now just make a filelist of sample graphs
make_sample_graph_filelist($outdir, $kmer, $comb_dir);


printf("\n\n****\nDONE.\nPlease note this has output a config file: $output_config which is needed as an argument by scripts/calling/genotype_1sample_against_sites.pl\n");



sub make_sample_graph_filelist
{
    my ($rcdir, $k, $combine_dir) = @_;
     
    my $reg = $rcdir."\*/binaries/cleaned/k".$k.'/'."\*ctx";

    my $listing_cmd = "ls $reg";
    my $listing = qx{$listing_cmd};
    my @sp = split(/\n/, $listing);
    my $outfile = $combine_dir."list_args_for_final_step";
    open(OUT, ">".$outfile)||die("Cannot open $outfile\n");
    foreach my $f (@sp)
    {
	if ($f=~ /\/([^\/]+)\.kmer/)
	{
	    my $id = $1;
	    my $sample_odir;
	    if ($f=~ /($rcdir.+)binaries/)
	    {
		$sample_odir = $1;
	    }
	    
	    my $sample_gtdir= $sample_odir."union_calls";
	    if (!(-d $sample_gtdir))
	    {
		my $c = "mkdir $sample_gtdir";
		qx{$c};
	    }
	    print OUT "$id\t$sample_gtdir\t$f\n";
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
    return ($ref_list, $ref_col_list);
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
