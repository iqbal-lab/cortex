#!/usr/bin/perl -w
use strict;
use File::Basename;

use Getopt::Long;

## Pass in a cleaned sample graph, and this will intersect it with the bubbles, and dump a diddy graph

my $cortex_dir = "/ddn/projects3/mcvean_res/zam/phase2_cortex/dummyrun/code/CORTEX_release_v1.0.5.8";
my $sample_cleaned_graph="";
my $mem_height = 23;
my $mem_width = 100;
my $kmer= 31;
my $sample="";
my $thresh=-1;
my $bubble_graph="";

&GetOptions(
    'sample:s'                            =>\$sample,
    'bubble_graph:s'                         => \$bubble_graph,
    'sample_cleaned_graph:s'                 =>\$sample_cleaned_graph,
    'mem_height:i'                          =>\$mem_height,
    'mem_width:i'                           =>\$mem_width,
    'thresh:s'                              =>\$thresh,
#        'kmer|k:i'                               => \$kmer,
    );




check_args($bubble_graph, $sample_cleaned_graph, $mem_height, $mem_width, $sample, $thresh);

if ($cortex_dir !~ /\/$/)
{
    $cortex_dir=$cortex_dir.'/';
}


## Now intersect your sample
my $ctx_binary = check_cortex_compiled_2colours($cortex_dir, $kmer);
my $colour_list = make_colourlist($sample_cleaned_graph, $sample);
my $suffix = "intersect_bubbles_thresh".$thresh;

my $cmd = $ctx_binary." --kmer_size $kmer --mem_height $mem_height --mem_width $mem_width --multicolour_bin $bubble_graph --colour_list $colour_list --load_colours_only_where_overlap_clean_colour 0 --successively_dump_cleaned_colours $suffix";
print "$cmd\n";
my $ret = qx{$cmd};
print "$ret\n";
print "Finished intersecting sample $sample with the bubble graph $bubble_graph\n";


sub check_args
{
    my ($bub_g, $sample_graph, $h, $w, $sam, $thre) = @_;
    
    if ($bub_g eq "")
    {
	die("You must specify the bubble graph, using --bubble_graph");
    }
    if (!(-e $bub_g))
    {
	die("You have specified a non existent file $bub_g");
    }
    if ($sam eq "")
    {
	die("You must specify sample id with --sample");
    }
    if ($sample_graph eq "")
    {
	die("You must specify the cleaned sample graph with --sample_cleaned_graph");
    }
    if ($sample_graph !~ /$sam/)
    {
	die("Suspicious - your sample graph $sample_graph does not pattern match with the name of the sample $sam - you must be using nonstandard naming. Please stick to th convention");
    }
    if ($sample_graph !~ /clean\d+.ctx$/)
    {
	die("You have specified a sample graph $sample_graph which should be a cleaned graph of a sample, but the filename does not contain the text cleanT (for some number T) - please stick to the defined file naming convention");
    }
    if ($sample_graph !~ /.ctx/)
    {
	die("You have specified a sample graph $sample_graph which should be a cleaned graph of a sample, but the filename does not contain\".ctx\" - please stick to the standard naming convention as defined in the recipe PDF ");
    }

    my $num_kmers = (2**$h) * $w;
    if ($num_kmers>=800000000)
    {
	print "You have specified mem height $h and mem_width $w, which will support $num_kmers kmers\n";
	print "This will require ";
	print 16 * $num_kmers/1000000000;
	print "  Gb\n";
    }
    else
    {
	die( "Ideally I would want to support over 800 million to be safe\n");
    }
    if ($thre==-1)
    {
	die("You must specify --thresh - it should be the threshold you used when cleaning the merged pool with --remove_low_coverage_supernodes\n");
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
	die("Please go and compile Cortex in $cortex_dir for k=$kmer with 2 colours");
    }
}


sub make_colourlist
{
    my ($graph, $id) = @_;

    my $bname = basename($graph);
    my $dir = dirname($graph);

    if ($dir !~ /\/$/)
    {
	$dir = $dir.'/';
    }
    my $col_list = $dir.$id."_colourlist_for_intersection";
    my $list_this_binary=$graph.".filelist";
    my $list_this_binary_nodir = basename($list_this_binary); 
    open(COL, ">".$col_list)||die("Cannot open $col_list");
    open(FLIST, ">".$list_this_binary)||die("Cannot open $list_this_binary");

    print FLIST "$bname\n";
    close(FLIST);
    print COL "$list_this_binary_nodir\n";
    close(FLIST);
    close(COL);

    return $col_list;
}
