#!/usr/bin/perl -w
use strict;
use File::Basename;
use File::Spec;
use Getopt::Long;
use Benchmark;

my $invcf = "";
my $bubble_callfile = "";

##the following are set by the command-line arguments
my $cortex_dir = "";
my $sample_cleaned_graph="";
my $ref_overlap_bubble_graph="";
my $mem_height = 16;
my $mem_width = 100;
my $kmer= 31;
my $sample="";
my $sample_graph="";
my $bubble_graph="";
my $overlap_log="";
my $outdir="";
my $filename_of_sample_overlap_bubble_binary="";
my $max_allele=100000; ##longest allele which will be genotyped
my $g_size = 0; #genome size
&GetOptions(
    ##mandatory args
    'cortex_dir:s'                        =>\$cortex_dir,
    'invcf:s'                             =>\$invcf,
    'bubble_graph:s'                      =>\$bubble_graph,
    'bubble_callfile:s'                   =>\$bubble_callfile,
    'outdir:s'                            =>\$outdir,
    'sample:s'                            =>\$sample,##sample ID
    'ref_overlap_bubble_graph:s'          =>\$ref_overlap_bubble_graph,
    'genome_size:i'                       =>\$g_size,
    ##if you have already overlapped the sample with the bubble graph
    'sample_overlap_bubbles:s'            =>\$filename_of_sample_overlap_bubble_binary,

    ##if you just have the sample graph and you want this script to overlap it with the bubbles
    'sample_graph:s'                      =>\$sample_graph,##whole genome graph of sample
    'overlap_log:s'                       =>\$overlap_log,#specify name for log file for dumping overlap of sample with bubbles 

    ## there are default mem height and width set, but you can also set them here
    'mem_height:i'                        =>\$mem_height,
    'mem_width:i'                         =>\$mem_width,
    'max_allele:i'                        =>\$max_allele
    );




if ($cortex_dir !~ /\/$/)
{
    $cortex_dir=$cortex_dir.'/';
}
if ($outdir !~ /\/$/)
{
    $outdir=$outdir.'/';
}

my $analyse_variants_dir = $cortex_dir."scripts/analyse_variants/";
check_args($cortex_dir, $bubble_graph, $mem_height, $mem_width, $sample, $overlap_log, $invcf, $g_size);




##if you have not passed in the overlap binary of sample with bubbles, then do it now:

my $ctx_binary = check_cortex_compiled_2colours($cortex_dir, $kmer);

if ($filename_of_sample_overlap_bubble_binary eq "")
{
## Now intersect your sample
    print("\n*************************\n");
    print("First overlap the sample with the bubbles:\n");

### note start time
    my $time_start_overlap = new Benchmark;
    

    my $suffix = "intersect_sites";
    my $colour_list;
    ($colour_list, $filename_of_sample_overlap_bubble_binary) = make_colourlist($sample_graph, $sample, $suffix, $outdir);
    
    my $cmd = $ctx_binary." --kmer_size $kmer --mem_height $mem_height --mem_width $mem_width --multicolour_bin $bubble_graph --colour_list $colour_list --load_colours_only_where_overlap_clean_colour 0 --successively_dump_cleaned_colours $suffix > $overlap_log 2>&1";
    print "$cmd\n";
    my $ret = qx{$cmd};
    print "$ret\n";
    print "Finished intersecting sample $sample with the bubble graph $bubble_graph\n";
    my $time_end_overlap=new Benchmark;
    my $time_taken_overlap=timediff($time_end_overlap,$time_start_overlap);
    print "Time taken to overlap sample with bubbles is ", timestr($time_taken_overlap), "\n";

}


print("\n*************************\n");
print "Now, genotype sample:\n";
### note start time
my $time_start_gt = new Benchmark;

my $bn = basename($bubble_callfile);

my $gt_output = $outdir.$bn.".".$sample.".genotyped";
my $gt_log = $outdir.$sample."_gt.log ";
my $g_colour_list = make_2colourlist($filename_of_sample_overlap_bubble_binary, $ref_overlap_bubble_graph, $sample);
my $gt_cmd = $ctx_binary." --kmer_size $kmer --mem_height $mem_height --mem_width $mem_width --colour_list $g_colour_list --max_read_len $max_allele --gt ".$bubble_callfile.",".$gt_output.",BC --genome_size $g_size --experiment_type EachColourADiploidSampleExceptTheRefColour --print_median_covg_only  --estimated_error_rate 0.01 --ref_colour 0  > $gt_log 2>&1";
print "$gt_cmd\n";
my $gt_ret = qx{$gt_cmd};
print "$gt_ret\n";
print "Genotyping done\n";
my $time_end_gt=new Benchmark;
my $time_taken_gt=timediff($time_end_gt,$time_start_gt);
print "Time taken to gt sample is ", timestr($time_taken_gt), "\n";



#### DUMP VCF
print("\n*************************\n");
print("Now dump a VCF\n");
### note start time
my $time_start_vcf = new Benchmark;


#make sample list
my $samplelist = $outdir."REF_AND_".$sample;
open(O, ">".$samplelist)||die("Unable to open $samplelist");
print O "REF\n$sample\n";
close(O);

my $vcf_cmd = "perl $analyse_variants_dir"."process_calls.pl --callfile $gt_output --callfile_log $gt_log --outvcf $sample --outdir $outdir --samplename_list $samplelist  --num_cols 2  --caller BC --kmer 31 --refcol 0 --ploidy 2  --vcf_which_generated_calls $invcf --print_gls yes > $outdir".$sample.".vcf.log";
print $vcf_cmd;
my $vcfret = qx{$vcf_cmd};
print $vcfret;

print "\nFinished! Sample $sample is genotyped and a VCF has been dumped\n";

my $time_end_vcf=new Benchmark;
my $time_taken_vcf=timediff($time_end_vcf,$time_start_vcf);
print "Time taken to vcf sample with bubbles is ", timestr($time_taken_vcf), "\n";

print("\n*************************\n");




sub check_args
{
    my ($ctx_dir, $bub_g, $h, $w, $sam, $overlap_log, $vcf, $g) = @_;
    if ($ctx_dir eq "")
    {
	die("You must specify the Cortex install directory with --cortex_dir");
    }
    if ($g==0)
    {
	die("You must specify the genome size with --genome_size, which is used to convert bases-loaded into depth loaded\n");
    }
    if ($overlap_log eq "")
    {
	die("You must specify the name of the log file which will be generated when overlapping the sample with the bubbles\n");
    }
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
    my ($graph, $id, $suffix, $odir) = @_;

    $graph = File::Spec->rel2abs($graph);
    my $bname = basename($graph);
    my $dir = dirname($graph);

    if ($dir !~ /\/$/)
    {
	$dir = $dir.'/';
    }
    if ($odir !~ /\/$/)
    {
	$odir = $odir.'/';
    }
    
    my $col_list = $odir.$id."_colourlist_for_intersection";
    my $list_this_binary=$odir.$bname.".filelist";
    my $list_this_binary_nodir = basename($list_this_binary); 
    open(COL, ">".$col_list)||die("Cannot open $col_list");
    open(FLIST, ">".$list_this_binary)||die("Cannot open $list_this_binary");

#    print FLIST "$bname\n";
    print FLIST "$graph\n";
    close(FLIST);
    print COL "$list_this_binary_nodir\n";
    #print COL "$list_this_binary\n";
    close(FLIST);
    close(COL);

    return ($col_list, $list_this_binary."_".$suffix.".ctx");
}


##ref and sample
#we just want a file that says "REF\nSAMPLE_ID\n";
#where 
#  >cat REF
#  /path/to/reference_overlap_bubbles_binary
#  >cat SAMPLE_ID
#  /path/to/sample_ontersect_bubbles_binary
sub make_2colourlist
{
    my ($graph, $ref_binary, $id) = @_;

    $graph = File::Spec->rel2abs($graph);
    $ref_binary = File::Spec->rel2abs($ref_binary);
    my $bname = basename($graph);
    my $ref_bname =  basename($ref_binary);
    my $dir = dirname($graph);
    my $refdir = dirname($ref_binary);
    if ($dir !~ /\/$/)
    {
	$dir = $dir.'/';
    }
    if ($refdir !~ /\/$/)
    {
	$refdir = $refdir.'/';
    }
    my $col_list = $dir.$id."_colourlist_for_genotyping";
    my $list_this_binary=$graph.".filelist";
    my $list_ref_binary =$ref_binary.".filelist";

    my $list_this_binary_nodir = basename($list_this_binary); 
    my $list_ref_binary_nodir = basename($list_ref_binary); 
    open(COL, ">".$col_list)||die("Cannot open $col_list");
    open(FLIST, ">".$list_this_binary)||die("Cannot open $list_this_binary");
    open(RLIST, ">".$list_ref_binary)||die("Cannot open $list_ref_binary");

    print FLIST "$bname\n";
    close(FLIST);
    print RLIST "$ref_bname\n";
    close(RLIST);
    print COL "$list_ref_binary\n$list_this_binary\n";
    close(COL);

    return $col_list;
}
