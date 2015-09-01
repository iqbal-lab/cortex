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
    $cortex_dir =~ s/scripts\/calling\/gt_1sample.pl//;
    push ( @INC, $cortex_dir."/scripts/calling/");
}
use BasicUtils qw ( add_slash );
use ConfigMgmt qw( get_from_config_if_undef print_to_config get_all_info_from_config_if_undef);

my %vars = ( "invcf" => "",
	     "sample_cleaned_graph"=> "",
	     "mem_height" => "",
	     "mem_width" => "",
	     "sample" => "",
	     "sample_graph" => "",
	     "outdir" => "",
	     "filename_of_sample_overlap_bubble_binary"=>"",
	     "genome_size"=>"",
	     "kmer"=>"",
	     "bubble_callfile" => "",
	     "max_allele" => "",
	     "bubble_graph" => "",
	     "ref_overlap_bubble_graph" => "");


my $config ="";

&GetOptions(
    ##mandatory args
    'invcf:s'                             =>\$vars{"invcf"},
    'config:s'                             =>\$vars{"config"},

    ##optional args - usually handled automatically by config files from previous scripts
    'outdir:s'                            =>\$vars{"outdir"},
    'sample:s'                            =>\$vars{"sample"},##sample ID
    'genome_size:i'                       =>\$vars{"genome_size"},

    ##if you just have the sample graph and you want this script to overlap it with the bubbles
    'sample_graph:s'                      =>\$vars{"sample_graph"},##whole genome graph of sample

    ## config files will set mem height and width set, but you can also set them here
    'mem_height:i'                        =>\$vars{"mem_height"},
    'mem_width:i'                         =>\$vars{"mem_width"},

    );


check_mandatory_args(\%vars);
get_all_info_from_config_if_undef(\%vars, $vars{"config"});


$cortex_dir = BasicUtils::add_slash($cortex_dir);
$vars{"outdir"} =  BasicUtils::add_slash($vars{"outdir"});
my $analyse_variants_dir = $cortex_dir."scripts/analyse_variants/";
check_args(\%vars);





##if you have not passed in the overlap binary of sample with bubbles, then do it now:

my $ctx_binary = check_cortex_compiled_2colours($cortex_dir, $vars{"kmer"});


## Now intersect your sample
print("\n*************************\n");
print("First overlap the sample with the bubbles:\n");

### note start time
my $time_start_overlap = new Benchmark;


my $suffix = "intersect_sites";
my $colour_list;

($colour_list, $vars{"filename_of_sample_overlap_bubble_binary"}) 
    = make_colourlist($vars{"sample_graph"}, 
		      $vars{"sample"}, 
		      $suffix, 
		      $vars{"outdir"});

my $overlap_log=$vars{"outdir"}."log_overlap_of_".$vars{"sample"}."_with_bubbles.txt";
my $cmd = $ctx_binary." --kmer_size ".$vars{"kmer"};
$cmd .= " --mem_height ".$vars{"mem_height"};
$cmd .= " --mem_width ".$vars{"mem_width"};
$cmd .= " --multicolour_bin ".$vars{"bubble_graph"};
$cmd .= " --colour_list $colour_list --load_colours_only_where_overlap_clean_colour 0 ";
$cmd .= " --successively_dump_cleaned_colours $suffix > $overlap_log 2>&1";

print "$cmd\n";
my $ret = qx{$cmd};
print "$ret\n";

print "Finished intersecting sample ".$vars{"sample"}." with the bubble graph ".$vars{"bubble_graph"}."\n";
my $time_end_overlap=new Benchmark;
my $time_taken_overlap=timediff($time_end_overlap,$time_start_overlap);
print "Time taken to overlap sample with bubbles is ", timestr($time_taken_overlap), "\n";


print("\n*************************\n");
print "Now, genotype sample:\n";
### note start time
my $time_start_gt = new Benchmark;

my $bn = basename($vars{"bubble_callfile"});

my $gt_output = $vars{"outdir"}.$bn.".".$vars{"sample"}.".genotyped";
my $gt_log = $vars{"outdir"}.$vars{"sample"}."_gt.log ";
my $g_colour_list = make_2colourlist($vars{"filename_of_sample_overlap_bubble_binary"}, 
				     $vars{"ref_overlap_bubble_graph"}, 
				     $vars{"sample"});
my $gt_cmd = $ctx_binary." --kmer_size ".$vars{"kmer"};
$gt_cmd .= " --mem_height ".$vars{"mem_height"};
$gt_cmd .= " --mem_width ".$vars{"mem_width"};
$gt_cmd .= " --colour_list $g_colour_list ";
$gt_cmd .= " --max_read_len ".$vars{"max_allele"};
$gt_cmd .= " --gt ".$vars{"bubble_callfile"}.",".$gt_output.",BC ";
$gt_cmd .= " --genome_size ".$vars{"genome_size"};
$gt_cmd .= " --experiment_type EachColourADiploidSampleExceptTheRefColour --print_median_covg_only ";
$gt_cmd .= " --estimated_error_rate 0.01 --ref_colour 0  > $gt_log 2>&1";
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
my $samplelist = $vars{"outdir"}."REF_AND_".$vars{"sample"};
open(O, ">".$samplelist)||die("Unable to open $samplelist");
print O "REF\n".$vars{"sample"}."\n";
close(O);

my $vcf_cmd = "perl $analyse_variants_dir"."process_calls.pl --callfile $gt_output --callfile_log $gt_log --outvcf ".$vars{"sample"};
$vcf_cmd .= " --outdir ".$vars{"outdir"};
$vcf_cmd .= " --samplename_list $samplelist --num_cols 2  --caller BC ";
$vcf_cmd .= " --kmer ".$vars{"kmer"};
$vcf_cmd .= " --refcol 0 --ploidy 2  --vcf_which_generated_calls ".$vars{"invcf"};
$vcf_cmd .= " --print_gls yes > ".$vars{"outdir"}.".".$vars{"sample"}.".vcf.log 2>&1";
print $vcf_cmd;
my $vcfret = qx{$vcf_cmd};
print $vcfret;

print "\nFinished! Sample ".$vars{"sample"}." is genotyped and a VCF has been dumped\n";

my $time_end_vcf=new Benchmark;
my $time_taken_vcf=timediff($time_end_vcf,$time_start_vcf);
print "Time taken to vcf sample with bubbles is ", timestr($time_taken_vcf), "\n";

print("\n*************************\n");



sub check_mandatory_args
{
    my ($hashref) = @_;
    if ($hashref->{"config"} eq "")
    {
	die("Must specify --config, the config file output by combine_vcfs.pl\n");
    }
    elsif (!(-e $hashref->{"config"}))
    {
	die("Cannot open the user specified file $config\n");
    }


    if ($hashref->{"invcf"} eq "")
    {
	die("Must specify --invcf, the VCF of sites produced by combine_vcfs.pl\n");
    }
    elsif (!(-e $hashref->{"invcf"}))
    {
	die("Cannot open the user specified file ".$hashref->{"invcf"});
    }


}
sub check_args
{
    my ($hashref) = @_;



    if ($hashref->{"kmer"} eq "")
    {
	my $errstr ="Kmer not specified on command-line of gt_1sample.pl,\n";
	$errstr .= "nor in config file ".$hashref->{"config"};
	$errstr .=".\n";
	die($errstr)
    }
    if ($hashref->{"bubble_callfile"} eq "")
    {
	my $errstr ="bubble_callfile not specified on command-line of gt_1sample.pl,\n";
	$errstr .= "nor in config file ".$hashref->{"config"};
	$errstr .=".\n";
	die($errstr)
    }
    if ($hashref->{"max_allele"} eq "")
    {
	my $errstr ="max_allele not specified on command-line of gt_1sample.pl,\n";
	$errstr .= "nor in config file ".$hashref->{"config"};
	$errstr .=".\n";
	die($errstr)
    }
    if ($hashref->{"bubble_graph"} eq "")
    {
	my $errstr ="bubble_graph not specified on command-line of gt_1sample.pl,\n";
	$errstr .= "nor in config file ".$hashref->{"config"};
	$errstr .=".\n";
	die($errstr)
    }
    if ($hashref->{"ref_overlap_bubble_graph"} eq "")
    {
	my $errstr ="ref_overlap_bubble_graph not specified on command-line of gt_1sample.pl,\n";
	$errstr .= "nor in config file ".$hashref->{"config"};
	$errstr .=".\n";
	die($errstr)
    }
    if ($hashref->{"genome_size"} eq "")
    {
	my $errstr ="genome_size not specified on command-line of gt_1sample.pl,\n";
	$errstr .= "nor in config file ".$hashref->{"config"};
	$errstr .=".\n If you use prepare.pl this is handled automatically\n";
	die($errstr)
    }

    if ($hashref->{"sample"} eq "")
    {
	die("You must specify sample id with --sample");
    }
    if ($hashref->{"sample_graph"} eq "")
    {
	die("You must specify the sample graph with --sample_graph");
    }

#$mem_height, $mem_width, $sample, $invcf, $g_size);
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
	die("Please go and compile Cortex in $cortex_dir for k=$k with 2 colours");
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
    my $list_this_binary=$odir.$id.".filelist";
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

sub get_args_from_config
{
    my ($cfile) = @_;
    open(CFILE, $cfile)||die("Cannot open $cfile\n");
    my @arr=();
    while (<CFILE>)
    {
	my $ln = $_;
	chomp $ln;
	my @sp = split("\t", $ln);

	if ($sp[0] eq "kmer")
	{
	    $arr[0]=$sp[1];
	}
	elsif ($sp[0] eq "bubble_callfile")
	{
	    $arr[1]=$sp[1];
	}
	elsif ($sp[0] eq "max_allele")
	{
	    $arr[2]=$sp[1];
	}
	elsif ($sp[0] eq "bubble_graph")
	{
	    $arr[3]=$sp[1];
	}
	elsif ($sp[0] eq "ref_overlap_bubble_graph")
	{
	    $arr[4]=$sp[1];
	}
	else
	{
	    die("Malformed config file\n");
	}


    }
    close(CFILE);
    return @arr;
}
