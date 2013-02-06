#!/usr/bin/perl -w
use strict;

use Getopt::Long;

## Pass in a list of uncleaned binaries - these will be cleaned against the pool, one after another.
## If you have enough memory/high-mem servers to run many of these at the same time,
##then you can of course parallelise. Each sample is cleaned completely independently of the others.
## However you do need to load the cleaned pool into memory, so you're not going to be able to do this on a 
## small cluster node.

my $cortex_dir = "/ddn/projects3/mcvean_res/zam/phase2_cortex/dummyrun/code/CORTEX_release_v1.0.5.8";
my $list_sample_uncleaned_bins="";
my $mem_height = 25;
my $mem_width = 150;
my $pop="";
my $kmer= 31;
my $cleaned_pool = "";
my $thresh=-1;
&GetOptions(
        'cleaned_pool:s'                         =>\$cleaned_pool,
        'list:s'                                 => \$list_sample_uncleaned_bins,
        'pop:s'                                  => \$pop, #population identifier, eg LWK
         'mem_height:i'                          =>\$mem_height,
         'mem_width:i'                           =>\$mem_width,
         'thresh:i'                              => \$thresh, #what threshold did you use when cleaning the merged pool?
#        'kmer|k:i'                               => \$kmer,
);




check_args($pop, $list_sample_uncleaned_bins, $mem_height, $mem_width, $cleaned_pool, $thresh);

if ($cortex_dir !~ /\/$/)
{
    $cortex_dir=$cortex_dir.'/';
}

my $ctx_binary = check_cortex_compiled($cortex_dir, $kmer);
my $colour_list = make_colourlist($list_sample_uncleaned_bins, $pop);
my $log = $colour_list.".log";
my $suffix = "clean".$thresh;
my $cmd = $ctx_binary." --kmer_size $kmer --mem_height $mem_height --mem_width $mem_width --multicolour_bin $cleaned_pool --colour_list $colour_list --load_colours_only_where_overlap_clean_colour 0 --successively_dump_cleaned_colours $suffix > $log 2>&1 ";
print "$cmd\n";
my $ret = qx{$cmd};
print "$ret\n";
print "Finished cleaning samples in $list_sample_uncleaned_bins\n";







sub check_args
{
    my ($p, $l, $h, $w, $cl_pl, $threshold) = @_;
    if ($p eq "")
    {
	die("You must specify --pop");
    }
    if ($cl_pl eq "")
    {
	die("You must specify --cleaned_pool and give the cleaned population pool binary\n");
    }
    if (!(-e $l))
    {
	die("You must specify --list , giving a single file listing all the uncleaned sample binaries that you want to clean - these are all cleaned independently, one after another, so you can parallelise"); 
    }
    open(L, $l)||die("Cannot open your list of uncleaned binaries, $l");
    while (<L>)
    {
	my $f = $_;
	chomp $f;
	if (!(-e $f))
	{
	    #die("Your list $l, contains a file $f that we cannot find on the filesystem - some kind of typo? Or a network.disk issue?");
	}
	elsif ($f !~ /uncleaned.q10.k31.ctx/)
	{
	    die("This list of uncleaned binaries contains this file $f, which does not match this text: uncleaned.q10.k31.ctx - is it really an uncleaned binary?");
	}
    }
    close(L);


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
	die("Your cleaned pool $cl_pl does not have a filename containing the text \"k31\" - are you sure you passed in the right thing?\n");
    }
    elsif ($cl_pl !~ /ctx/)
    {
	die("Your cleaned pool $cl_pl does not have a filename containing the text \"ctx\" - are you sure you passed in the right thing?\n");
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


sub make_colourlist
{
    my ($list, $pop) = @_;

    my $col_list = $pop.$list.".colourlist_for_cleaning_samples";
    open(COL, ">".$col_list)||die("Cannot open $col_list");
    open(LIST, $list)||die("Cannot open your filelist $list");
    while(<LIST>)
    {
	my $line = $_;
	chomp $line;
	my $list_this_binary=$line.".filelist";
	open(FLIST, ">".$list_this_binary)||die("Cannot open $list_this_binary");
	print FLIST "$line\n";
	close(FLIST);
	print COL "$list_this_binary\n";
    }
    close(LIST);
    close(COL);
    return $col_list;
}
