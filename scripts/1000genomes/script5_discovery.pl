#!/usr/bin/perl -w                                                                                                                                                         
use strict;

use Getopt::Long;

## Pass in the reference and cleaned pool                                                                                                                                  

my $cortex_dir = "/ddn/projects3/mcvean_res/zam/phase2_cortex/dummyrun/code/CORTEX_release_v1.0.5.8";
my $mem_height = 25;
my $mem_width = 150;
my $pop="";
my $kmer= 31;
my $cleaned_pool = "";
my $ref = "";
my $thresh=-1;
my $outdir = "";

&GetOptions(
         'cleaned_pool:s'                         =>\$cleaned_pool,### binary file of the cleaned population pool
         'pop:s'                                  => \$pop, #population identifier, eg LWK                                                                                 
         'thresh:i'                              => \$thresh, #what threshold did you use when cleaning the merged pool?                                                   
         'ref:s'                                 => \$ref,## full path to reference genome human_g1k_v37.proper_chroms.k31.ctx
         'outdir:s'                              =>\$outdir,

    ##the following two options have a default set that should be fine                                                                                                     
         'mem_height:i'                          =>\$mem_height,
         'mem_width:i'                           =>\$mem_width,

#        'kmer|k:i'                               => \$kmer,                                                                                                               
    );

check_args($pop, $cleaned_pool, $ref, $mem_height, $mem_width, $thresh, $outdir);
if ($outdir !~ /\/$/)
{
    $outdir = $outdir.'/';
}

if ($cortex_dir !~ /\/$/)
{
    $cortex_dir=$cortex_dir.'/';
}

my $ctx_binary = check_cortex_compiled($cortex_dir, $kmer);
my $colour_list = make_colour_list($cleaned_pool, $ref, $pop, $outdir);

my $bubbles = $outdir.$pop."_bubbles_thresh".$thresh;
my $log = $bubbles.".log";
my $cmd = $ctx_binary." --kmer_size $kmer --mem_height $mem_height --mem_width $mem_width --colour_list $colour_list --detect_bubbles1 1/1 --output_bubbles1 $bubbles  --exclude_ref_bubbles  --ref_colour 0 > $log 2>&1";
print "$cmd\n";
my $ret = qx{$cmd};
print "$ret\n";

sub make_colour_list
{
    my ($pool, $refr, $popn, $odir) = @_;
    if ($odir !~ /\/$/)
    {
	$odir = $odir.'/';
    }
    
    my $colour_list = "discovery_colour_list";
    my $ref_list = "ref_binary_list";
    my $pool_list = "pool_binary_list";
    open(COL, ">".$colour_list)||die("Cannot open $colour_list");
    print COL "$ref_list\n";
    print COL "$pool_list\n";
    close(COL);
    open(REF, ">".$ref_list)||die("Cannot open $ref_list");
    print REF "$refr\n";
    close(REF);
    open(POOL, ">".$pool_list)||die("Cannot open $pool_list");
    print POOL "$pool\n";
    close(POOL);
    return $colour_list;
}


sub check_args
{
    my ($p, $cl_pl, $REF, $h, $w,  $threshold, $odir) = @_;

    if ($odir eq "")
    {
	die("You must specify the output dir where you want the variant calls to go, with --outdir");
    }

    if (!(-d $odir))
    {
	print ("WARNING - You have specified a non-existent directory - $odir - creating it and carrying on\n");
	my $c = "mkdir -p $odir";
	check_arqx{$c};
    }
    
    if ($REF eq "")
    {
	die("You must specify the reference binary with --ref, and its filename must be human_g1k_v37.proper_chroms.k31.ctx - obviously thre path depends on where you have put it\n");
    }
    if ($REF !~ /human_g1k_v37.proper_chroms.k31.ctx/)
    {
	die("You have specified a non-standard reference $REF - we expect you to be using the same reference as everyone else - human_g1k_v37.proper_chroms.k31.ctx\n");
    }
    if ($p eq "")
    {
	die("You must specify --pop");
    }
    if ($cl_pl eq "")
    {
	die("You must specify --cleaned_pool and give the cleaned population pool binary\n");
    }


    my $num_kmers = (2**$h) * $w;
    if ($num_kmers>=4000000000)
    {
	print "You have specified mem height $h and mem_width $w, which will support $num_kmers kmers\n";
	print "This will require ";
	print 24 * $num_kmers/1000000000;
	print "  Gb\n";
    }
    else
    {
	die( "Ideally I would want to support over 4 billion to be safe\n");
    }

    if (!(-e $cl_pl))
    {
	die("This file you specify - $cl_pl - does not exist - typo? network error? disk error?\n")
    }
    if ($cl_pl !~ /k31/)
    {
	die("Your cleaned pool $cl_pl does not have a filename endingin containing the text \"k31\" - are you sure you passed in the right thing?\n");
    }
    elsif ($cl_pl !~ /ctx/)
    {
	die("Your cleaned pool $cl_pl does not have a filename endingin containing the text \"ctx\" - are you sure you passed in the right thing?\n");
    }
    elsif ($cl_pl !~ /thresh/)
    {
	die("Your cleaned pool does not have a filename containing the text \"thresh\" - either you passed in the wrong thing, or you are not sticking to the filenaming conventions specified in Zam's recipe\n");
    }
    

    if ($threshold==-1)
    {
	die("You must specify --thresh - it should be the threshold you used when cleaning the merged pool with --remove_low_coverage_supernodes\n");
    }
}



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
	
    if (-e $dir."bin/cortex_var_".$maxk."_c2")
    {
	return $dir."bin/cortex_var_".$maxk."_c2";
    }
    else
    {
	die("Please go and compile Cortex in $cortex_dir for k=$kmer with 2 colours");
    }
}

