#!/usr/bin/perl -w
use strict;

use Getopt::Long;
use File::Basename;
## Pass in a list of CLEANED sample  binaries

my $cortex_dir = "/home/zam/dev/hg/bitbucket/CORTEX_mainline/";
my $mem_height = 23;
my $mem_width = 100;
my $pop="";
my $kmer= 31;
my $bubble_calls="";
my $ref = "";
&GetOptions(
    'ref:s'                                  =>\$ref,
    'bubble_calls:s'                         =>\$bubble_calls,
    'mem_height:i'                           =>\$mem_height,
    'mem_width:i'                            =>\$mem_width,
#        'kmer|k:i'                               => \$kmer,
    );




check_args($mem_height, $mem_width, $bubble_calls, $ref);

if ($cortex_dir !~ /\/$/)
{
    $cortex_dir=$cortex_dir.'/';
}


## make the bubble-branch fasta and dump a graph

my $branches = $bubble_calls.".branches.fasta";
my $branches_log = $branches.".log";
my $cmd1 = "perl $cortex_dir"."scripts/analyse_variants/make_branch_fasta.pl --callfile $bubble_calls --kmer 31 > $branches_log 2>&1 ";
print "$cmd1\n";
my $ret1 = qx{$cmd1};
print "$ret1\n";

if (!(-e $branches))
{
    die("Failed to create $branches. See logfile $branches_log\n");
}
my $branches_list = $branches.".list";
my $cmd2 = "ls $branches > $branches_list";
qx{$cmd2};


my $ctx_binary1 = check_cortex_compiled_1colour($cortex_dir, $kmer);
my $bubble_graph = $branches;
$bubble_graph =~ s/fasta/k31.ctx/;
my $bubble_graph_log = $bubble_graph.".log";
my $cmd3 = $ctx_binary1." --se_list $branches_list --mem_height $mem_height --mem_width $mem_width --kmer_size 31  --dump_binary $bubble_graph  > $bubble_graph_log 2>&1";
qx{$cmd3};

print "Finished building a graph just of the bubble branches/alleles. Now intersect the ref binary\n";


my $ctx_binary2 = check_cortex_compiled_2colours($cortex_dir, $kmer);
my $suffix = "intersect_bubbles";
my $ref_intersect_log = basename($ref.".intersect_with_bubbles.log");
my $ref_col_list=get_ref_col_list($ref);

my $cmd4 = $ctx_binary2." --kmer_size 31 --multicolour_bin $bubble_graph --mem_height 23 --mem_width 100 --colour_list $ref_col_list  --load_colours_only_where_overlap_clean_colour 0 --successively_dump_cleaned_colours $suffix  > $ref_intersect_log 2>&1";
qx{$cmd4};




sub get_ref_col_list
{
    my ($ref) = @_;
    my $ref_col_list = basename($ref.".colour_list");
    my $ref_list = basename($ref.".list");
    my $c1 = "ls $ref > $ref_list";
    qx{$c1};
    my $c2 = "ls $ref_list  > $ref_col_list";
    qx{$c2};
    return $ref_col_list;
}

sub check_args
{
    my ($h, $w, $calls, $r) = @_;

    if ($r eq "")
    {
	die("You must specify the reference binary with --ref, so we can dump the intersection of this with the bubbles\n");
    }
#    if ($r !~ /human_g1k_v37.proper_chroms.k31.ctx/)
#    {
#	die("You have specifed a ref binary which does not pattern match with human_g1k_v37.proper_chroms.k31.ctx - please don't change the standardised naming\n");
#    }
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
	#die( "Ideally I would want to support over 800 million to be safe\n");
    }

    if ($calls eq "")
    {
	die("You must specifty --calls, and give as argument the output file from the Bubble Caller");
    }
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

