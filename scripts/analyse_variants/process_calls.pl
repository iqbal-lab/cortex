#!/usr/bin/perl -w

use strict;
use File::Basename;
use Getopt::Long;


#### Dear User - this following line is the only one in this script you may want to edit
my $stampy_bin =
  "/home/zam/installed_apps/stampy-1.0.13/stampy.py";    #### <<<<<<< can also use command line parameter to specify this
my $vcftools_dir="/path/to/vcftools/dir";### <<<< can also use command line param

#### no need to modify anything below this line

# However you may decide you want to increase this threshold
my $mapping_qual_thresh = 40;    #  - demand 5prime flank maps with quality>= 40


my $script_dir;
my $cortex_dir;
my $isaac_bioinf_dir;

my $str_input_args = "*********** Command-line used was : *********************\nperl ".$0." ".join(" ",@ARGV)."\n*********************************************************\n";

BEGIN
{
	use FindBin;
	$script_dir = $FindBin::Bin;
	if ($script_dir !~ /\/$/)
	{
	    $script_dir=$script_dir.'/';
	}
	$isaac_bioinf_dir = $script_dir."bioinf-perl/";
	$cortex_dir = $script_dir . '/../../';
	$cortex_dir = $script_dir;
	$cortex_dir =~ s/scripts\/analyse_variants//;
	
	push( @INC,
		$script_dir
		  . "/perl_modules/Statistics-Descriptive-2.6",
	      $isaac_bioinf_dir."lib/"
	);
}

my $check_perl5 = "echo \$PERL5LIB";
my $check_perl5_ret = qx{$check_perl5};
my $isaac_libdir = $isaac_bioinf_dir."lib";
if ($check_perl5_ret !~ /$isaac_libdir/)
{
    $ENV{PERL5LIB} .= ":$isaac_libdir";
}

use lib $cortex_dir."/scripts/analyse_variants/perl_modules/Statistics-Descriptive-2.6";
use Descriptive;


###******** Description of process_calls.pl ##################################################################
# process_call.pl Takes as input a set of variant calls made by Cortex, maps the flanks, aligns
# branches against each other to determine variant type, splits out SNPs from clusters so we have precise loci,
# 
##############################################################################################################



##make sure there is forward slash on the end:
if ( $cortex_dir !~ /\/$/ )
{
	$cortex_dir = $cortex_dir . '/';
}

my (
	$callfile,             $outdir,
	$outvcf_filename_stub, $colours,
	$number_of_colours,    $reference_colour,
	$classif,              $prefix,
	$ploidy,               $vcf_on_which_calls_based,##if we build branches from a VCF and then genotype, and then pass to process_calls, use vcf_on_which_calls_based
	$stampy_hash_stub,     $ref_fasta, $unioncalls, $caller_type,
        $callfile_log,         $kmer, $global_var_ctr,
        $print_gls,            $callfile_is_median_format
);

#set defaults
my $pooled_colour = -1;    #deprecated
$classif                             = -1;
#$require_one_allele_is_ref           = "yes";
$prefix                              = "";
$outdir                              = ".";
$outvcf_filename_stub                = '';
$callfile                            = "unspecified";
$colours                             = '';
$number_of_colours                   = 0;
$reference_colour                    = -1;
$ploidy                              = 2;
#$apply_filter_one_allele_must_be_ref = "unknown";
$stampy_hash_stub                    = "";
$ref_fasta                           = "unspecified";
$unioncalls                          = "unspecified";
$caller_type                         = "unspecified";
$callfile_log                        = "unspecified";
$vcf_on_which_calls_based            = "unspecified";
$kmer                                = -1;
$global_var_ctr                      =0;
$print_gls                           ='';
$callfile_is_median_format           =0;
my $help = '';    #default false


&GetOptions(
	'callfile|f:s'                         => \$callfile,
	'callfile_log|i:s'                     => \$callfile_log,
	'outdir|o:s'                           => \$outdir,
	'outvcf|v:s'                           => \$outvcf_filename_stub,
	'samplename_list|s:s'                  =>  \$colours,    #list of names of colours/samples, one per line. 
	'num_cols|n:i'                         => \$number_of_colours,
	'refcol|r:i'                           => \$reference_colour,
	                                          # ignore this colour for the VCF - dont 
                                                  #print out anything. if there is no reference in your colours, use -1
	'vcf_which_generated_calls|a:s'        => \$vcf_on_which_calls_based,
                                                ##so "bubbles" always have ref allele first, and we already know chr/pos/type of var, everything except genotypes in fact
	'pop_classifier|c:s'                   => \$classif,
	                                          ## file containing output of the population filter (classifier.R), or -1 if not used (default)
	'prefix|p:s'                           => \$prefix,            ## this allows you to add a prefix to any var name. By default, it uses what is in the cortex callfile
	'ploidy|y:i'                           => \$ploidy,            ## must be 1 or 2
	'stampy_hash|t:s'                      => \$stampy_hash_stub,
	                                            #stampy creates blah.sthash and blah.stidx. You should enter --stampy_hash blah
	'stampy_bin|b:s'                       => \$stampy_bin,    ## must be 1 or 2
        'ref_fasta|e:s'                        => \$ref_fasta, ## optional, but if there, will remove calls which have been misplaced by mapper.
        'vcftools_dir|d:s'                     => \$vcftools_dir,    #mandatory
        'unioncalls|u:s'                       => \$unioncalls, ## purely for internal use - the union callset which is passed in for genotyping (which generates $callfile),
                                                                ##  is sometimes annotated with KMER, which we may want to preserve
        'caller|l:s'                           => \$caller_type,# BC or PD
        'kmer|k:i'                             => \$kmer,
        'global_var_ctr:i'                      => \$global_var_ctr,## for internal use only, used by run_calls.pl
        'print_gls|q:s'                         => \$print_gls,
	'help'                                 => \$help,
);

if ($help)
{
	print "\n\n";
	print "Usage:\n********  mandatory arguments *********  :\n";
	print
"--callfile                    : file of calls output by Cortex (may be from Bubble or Path Divergence caller, but you MUST have used --print_colour_coverages\n";	print
"--callfile_log                : file in which you have saved the Cortex output (stuff it prints to screen when it runs)\n";
	print
"--outvcf                      : the output VCF files will have filenames starting with this\n";
	print
"--outdir                      : all output will go here. Default is current working directory\n";
	print
"--samplename_list             : file listing names of each colour/sample, one line per colour. These names end up in the header line of the VCF\n";
	print
"--num_cols                    : nUmber of colours in your graph (and callfile output)\n";
	print
"--stampy_bin                  : Path to stampy bin - default is set via variable in script. \n";

	print
"--stampy_hash                 : VCF format needs a chromosome and position. You have two options:\n";
	print
"                                 1. Use a reference genome to place yout variants. If you do this, process_calls.pl will\n";
	print
"                                    map your calls to it using Stampy. You need to first build a Stampy hash. \n";
	print
"                                    Stampy creates blah.sthash and blah.stidx. You should enter --stampy_hash blah\n";
	print
"                                 2. You have no reference, and don't care about coordinates, but want to know what the\n";
	print
"                                    variants are (SNPs, indels, complex), and who has what allele (genotypes)\n";
	print
"                                    In that case we already know we are going to be abusing VCF slightly - this is fine, we are pragmatic.\n";
	print
"                                    In this case (suppose you have N variant calls), create a pseudo-reference which has N chromosomes\n";
	print
"                                    Chromosome i is the 5prime flank, branch1 and then 3prime flank of variant i\n";
	print
"                                    Then build a stampy hash of this pseudo-reference and pass it in with the --stampy_hash argument as in 1 above\n";
	print
"                                 3. You have no decent reference but you have a rough consensus draft assembly - same as case 1\n";
	print
"                                 In all cases, you must use Stampy (http://www.well.ox.ac.uk/project-stampy)\n";

	print
 "--vcftools_dir                : VCFtools is used to generate VCFs - mandatory to either specify this on cmd-line, or manually edit the path at the top of this script\n";
	print 
	    "                                 This should be the VCFtools root directory, which has subdirectories called: bin,  cpp,  lib ..\n";
	print
	    "--caller                      : Which caller generated this callset. Acceptable arguments are BC or PD.\n";
	print
	    "--kmer                        : What kmer value was used in generating this callfile\n";
	print "*****   Optional arguments ******* \n";
	print
"--refcol                      : if one of the colours is the reference genome, specify this colour number. Default is -1 (meaning no reference present), but if you \n";
	print
"                                do have a reference genome and are producing a VCF with respect to it, you will ONLY get the correct ref-allele in the VCF\n";
	print
"                                if you have put the reference into the graph and specify which colour using --refcol\n";
	print
"--pop_classifier              : If you used classifier.R, give the filename of the output file\n";
	print
"--ploidy                      : Acceptable values are 1 and 2. Default is 2.\n";
	print
"--prefix                      : String prefix which will go in the front of any variant names. e.g --prefix ZAM will produce variants ZAM_var_1, ZAM_var_2, etc\n";
	print
"--ref_fasta                   : Stampy maps calls to a reference with a mapping quality. We use a threshold of 40 by default, so 1 in 10000 are wrongly placed on the reference\n";
	print
	    "                                If you pass in the name of the reference fasta here, this script will check the VCF and remove misplaced variants\n";
	print
"--vcf_which_generated_calls   : Sometimes we make calls on many samples, and generate a union sitelist VCF. We then make a branches file from this, and genotype all\n";
	print
	    "                                the samples at these sites. In this case there is no need to map the flanks, align branches etc etc - the VCF is already built, all\n";
	print
	    "                                we need to do is add the genotypes\n";
	print "\n\n\n";
	exit();
}


##print out what the user typed in
print "\n\n$str_input_args\n";

### checks
if ($callfile eq "unspecified")
{
    die("You must specify --callfile");
}
$callfile_is_median_format=check_callfile_format($callfile);
if ($callfile_log eq "unspecified")
{
    die("You must specify --callfile_log");
}
if ($caller_type eq "unspecified")
{
    die("You must specify --caller, as BC or PD\n");
}
if (!(-e $callfile))
{
    die("Cannot open callfile $callfile");
}
my $legacy_callfile = check_if_callfile_in_legacy_format($callfile);
if ($legacy_callfile==-1)
{
    die("Failed to run check on callfile, to see if it is in legacy format\n");
}
elsif ($legacy_callfile==1)
{
    print "This looks like the old format in which Cortex used to output calls. With things like > branch_200_1 for the first branch of var 200. \n";
    print " Should be no problem, have tried to get process_calls.pl to handle the old legacy format. This is not an error message\n";
}

if ( ($vcftools_dir eq "/path/to/vcftools/dir") && ($vcf_on_which_calls_based eq "unspecified") )
{
    die("You must specify the VCFTools directory, either on the commandline with --vcftools_dir, or by editing process_calls.pl manually (it's highlighted for you at the top of the file)\n");
}

if ($vcf_on_which_calls_based ne "unspecified")
{
    if (!(-e $vcf_on_which_calls_based))
    {
	die("Cannot find the VCF file on which calls are based:  $vcf_on_which_calls_based\n");
    }
}

if ( ( !( -e $stampy_bin )) && ($reference_colour != -1) && ($vcf_on_which_calls_based eq "unspecified") )
{ 
    print(
	"Since you have specified a reference colour (--refcol), I assume you are using a reference. In this case, this script requires Stampy be installed, but it cannot find it. I tlooked for $stampy_bin which is either the default (Zam's path) or the file you entered as the argument of --stampy_bin. Please can you re-run, entering a valid path to Stampy. eg /home/myname/path/to/stampy.py\n"
	);
	die();
}
if (!( -d $cortex_dir ) )
{
    die("I have tried to get the path to your Cortex directory, but I seem to have failed. I think it is $cortex_dir, but this does not exist (or your fileserver suddenly went down). Please email zam and tell him what happened\n");
}

if ( $number_of_colours == 0 )
{
	print("You must specify --num_cols\n");
	die();
}

if (   ( $outvcf_filename_stub eq '' )
	|| ( $callfile eq '' )
	|| ( $colours  eq '' )
        )
{
	die(
"You must specify --outvcf and --callfile and --samplename_list  "
	);
	exit(1);
}
if ( !( -e $callfile ) )
{
	die("Cannot find this callfile $callfile");
}
if ( ( !( -e $colours ) ) && ( $colours !~ /,/ ) )
{
	die("Cannot find this sample list $colours");
}
if ( ( $reference_colour < 0 ) && ( $reference_colour != -1 ) )
{
	die(
"Reference colour must be -1 (meaning there is no reference), or >=0 - you entered $reference_colour"
	);
}
if ( ( $reference_colour!=-1 ) && ( $stampy_hash_stub eq "" ) )
{
    if ($vcf_on_which_calls_based eq "unspecified")
    {
	die(
	    "If you are using a reference in the graph, then you must provide a stampy hash\n"
	    );
    }
    else
    {
	#everything is based on an input VCF which already has coords - no need to map
    }
}

if ( ( $reference_colour==-1 ) && ( $stampy_hash_stub ne "" ) )
{
    print("You can either use a reference, or not. If you use a reference, specify --refcol AND --stampy_hash. \n");
    die("If you have not used a reference in the Cortex graph, then do not specify --refcol and do not specify --stampy_hash\n");
}

if ( ( $ploidy != 1 ) && ( $ploidy != 2 ) )
{
	die("Must have ploidy 1 or 2");
}

my $flank_bin = $cortex_dir . "/scripts/analyse_variants/make_5p_flank_file.pl";
if ( !( -e $flank_bin ) )
{
	die(
"Cannot find make_5p_flank_file.pl, looked here: $flank_bin. Maybe you mis-entered the path to your Cortex directory when you manually edited process_calls?"
	);
}

my $needleman_wunsch_bin = $cortex_dir
  . "/scripts/analyse_variants/needleman_wunsch/needleman_wunsch";
if ( !( -e $needleman_wunsch_bin ) )
{
	print("Cannot find the needleman_wunsch binary. Either:\n");
	print(
"1. you mis-entered the path to your Cortex directory when you manually edited process_calls.pl, or\n"
	);
	print(
"2. you forgot to cd into scripts/analyse_variants/needleman_wunsch and type \"make\", as documented in the compilation section of the Manual\n"
	);
	die();
}


if (   ( $caller_type eq "PD"  )
	&& ( $reference_colour == -1 ) )
{
	die(
"If you have specified the caller as PD, you must specify the reference colour"
	);
}

########
#### checks
########

if ( !( -e $flank_bin ) )
{
	die("Cannot find $flank_bin");
}

#if (!(-e $process_bubbles_bin))
#{
#    die("Cannot find $process_bubbles_bin");
#}
if ( !( -e $needleman_wunsch_bin ) )
{
	die("Cannot find $needleman_wunsch_bin");
}

if ( !( -e $callfile ) )
{
	die("Cannot find $callfile");
}

if ( ( $classif ne "-1" ) && ( !( -e $classif ) ) )
{
	die(
"Cannot find $classif, which you have entered as the population filter file. Either enter the correct file (and path) or use -1 to signify you are not using it"
	);
}

if ( $outdir !~ /\/$/ )
{
	$outdir = $outdir . '/';
}
if ( !( -e $outdir ) )
{
    my $c = "mkdir $outdir";
    qx{$c};
    if ( !( -e $outdir ) )
    {
	die("Output directory $outdir does not exist, and I cannot create it. Permissions? Out of disk?");
    }
}

## You will often need to print a newline after the last colour which is not the reference or pool.
my $last_sample_col = $number_of_colours - 1;
my $f;
for ( $f = $number_of_colours - 1 ; $f >= 0 ; $f-- )
{
	if ( ( $f != $reference_colour ) && ( $f != $pooled_colour ) )
	{
		$last_sample_col = $f;
		last;
	}
}

## if you have passed in a unioncalls file, get the KMER information from it
my %call_to_extra_info=();
if ($unioncalls ne "unspecified")
{
    get_extra_info_for_calls(\%call_to_extra_info, $unioncalls);
}

###### Get the read length and covg info from the caller log file
my %colour_to_readlen=();
my %colour_to_totalseq=();
get_covg_and_total_seq_from_log($callfile_log, \%colour_to_readlen, \%colour_to_totalseq);


## 1. Map 5p flanks
my $bname     = basename($callfile);
my $flankfile = $outdir . $bname . ".5pflanks";
my $mapped_flanks = $flankfile . ".sam";
my %var_name_to_cut_flank = ();
my %var_name_to_flank_mq_filter           = ();

## Only map the flanks if there is a reference to map to!
## ALSO - do not map if we already have a VCF (on which these call branches were based)
if ( ($reference_colour !=-1) && ($vcf_on_which_calls_based eq "unspecified") )
{
    if ( !-e $flankfile )
    {
	my $flank_cmd = "perl $flank_bin $callfile > $flankfile";
	print "$flank_cmd\n";
	my $flank_ret = qx{$flank_cmd};
	print "$flank_ret\n";
    }
    
    if ( ( !( -e $flankfile ) ) || ( -z $flankfile ) )
    {
	die( "$flankfile either failed to be created, or was created with size zero"
	    );
    }
    
    
    if ( -e $mapped_flanks )
    {
	print("$mapped_flanks already exists, so wont re-do it\n");
    }
    else
    {
	my $map_flanks_cmd =
	    "$stampy_bin -g $stampy_hash_stub -h $stampy_hash_stub --norefoutput --inputformat=fasta -M $flankfile -o $mapped_flanks";
	print "$map_flanks_cmd\n";
	my $map_flanks_ret = qx{$map_flanks_cmd};
	print "$map_flanks_ret\n";
    }
    get_list_vars_with_cut_flanks( $mapped_flanks, \%var_name_to_cut_flank );
    filter_by_flank_mapqual( $mapped_flanks, \%var_name_to_flank_mq_filter );
}



## 2. Align branches against each other (but not if we already have a VCF on which the call branches were based)
my $proc_bub_output = $outdir . $bname . ".aligned_branches";

if ( ( !( -e $proc_bub_output ) ) && ($vcf_on_which_calls_based eq "unspecified") )
{
	wrap_needleman( $needleman_wunsch_bin, $callfile, $prefix,
		$proc_bub_output );
}
else
{
    if (-e $proc_bub_output )
    {
	print "Not aligning branches against each other, as $proc_bub_output already exists\n";
    }
    else
    {
	print "Not aligning branches against each other, as these calls are based on this VCF $vcf_on_which_calls_based, which already has implicitly done this\n";
    }
}

##  3. Apply filters, and collect a list of good calls.


my %var_name_to_covg_and_branch_filter    = ();
my %var_name_to_pd_filter                 = ();
my %var_name_to_combined_filtering_result = ();
my %pop_classifier                        = ();
my %pop_classifier_confidence             = ();

if ( $classif ne "-1" )
{
	get_pop_filter_info( $classif, \%pop_classifier,
			     \%pop_classifier_confidence );
}

if ($caller_type eq "PD")
{
    filter_calls_where_entire_ref_allele_is_multicopy(\%var_name_to_pd_filter, $callfile, $reference_colour);
}
combine_all_filters(
    \%var_name_to_covg_and_branch_filter,    \%var_name_to_flank_mq_filter,
    \%var_name_to_pd_filter,
    \%var_name_to_combined_filtering_result, \%pop_classifier
);

my $fh_calls;
my $fh_map_flanks;
my $fh_proc_bub;
my $fh_simple_vcf;
my $fh_decomp_vcf;

## open file handles
open( $fh_calls,      $callfile )        || die("Cannot open $callfile");
if ($vcf_on_which_calls_based eq "unspecified")
{
    open( $fh_proc_bub,   $proc_bub_output ) || die("Cannot open $proc_bub_output");
}

my $simple_vcf_name;
my $decomp_vcf_name;
if ($vcf_on_which_calls_based eq "unspecified")
{
    $simple_vcf_name = $outdir . $outvcf_filename_stub . ".raw.vcf.uncleaned";
    $decomp_vcf_name = $outdir . $outvcf_filename_stub . ".decomp.vcf.uncleaned";
    open( $fh_decomp_vcf, "> $decomp_vcf_name" )  || die("Cannot open $decomp_vcf_name");
}
else
{
    ## will directly print one single final file
    $simple_vcf_name = $outdir . $outvcf_filename_stub . ".vcf";
}

open( $fh_simple_vcf, "> $simple_vcf_name" ) || die("Cannot open $simple_vcf_name");


if ( ($reference_colour!=-1) && ($vcf_on_which_calls_based eq "unspecified") )
{
    open( $fh_map_flanks, $mapped_flanks )   || die("Cannot open $mapped_flanks");
}


##print vcf header 
my $header = get_vcf_header($colours);
print $fh_simple_vcf $header;

if ($vcf_on_which_calls_based eq "unspecified")
{
    print $fh_decomp_vcf $header;
}

my $ret = 1;
if ($vcf_on_which_calls_based eq "unspecified")
{
    while ( $ret == 1 )
    {
	
	$ret = print_next_vcf_entry_for_easy_and_decomposed_vcfs(
	    $fh_calls,               $fh_map_flanks,
	    $fh_proc_bub,            \%var_name_to_combined_filtering_result,
	    \%var_name_to_cut_flank, 1,
	    1,                       $fh_simple_vcf,
	    $fh_decomp_vcf,          \%pop_classifier_confidence,
	    $reference_colour
	    );
    }
    close($fh_calls);
    close($fh_proc_bub);
    if ($reference_colour !=-1)
    {
	close($fh_map_flanks);
    }
    
   ### cleanup
    
   ## for those non-SNPs where the called variant was in the reverse direction, and the 3p flank was zero-length
   ## we have had to put a Z in place of the base before the variant, in both ref and alt alleles. We need now
   ## to go through and fix these. Parse the VCF once and collect a list of chrom_pos where we have Z's. We need to do
   ## this for both raw and decomp VCFs. 
    my %Z_posns_raw=();## will be chr_pos --> 1 - at these positions the VCFs have a Z
    my %Z_posns_decomp=();
    
   #get_Z_positions(\%Z_posns_raw,    $simple_vcf_name);
   #get_Z_positions(\%Z_posns_decomp, $decomp_vcf_name);
    
    
   #get_characters_to_replace_Z(\%Z_posns_raw,    $ref_fasta);
   #get_characters_to_replace_Z(\%Z_posns_decomp, $ref_fasta);
    
    
   #my $fixed_Z_simple = $simple_vcf_name;
   #$fixed_Z_simple =~ s/\.uncleaned/.fixed_Z/;
   #my $fixed_Z_decomp = $decomp_vcf_name;
   #$fixed_Z_decomp =~ s/\.uncleaned/.fixed_Z/;
   #print "Start fix z\n";
   #fix_Z($simple_vcf_name, $fixed_Z_simple, \%Z_posns_raw);
   #fix_Z($decomp_vcf_name, $fixed_Z_decomp, \%Z_posns_decomp);
   #print "end fix z\n";
    
   #my $final_simple = $fixed_Z_simple;
   #$final_simple =~ s/\.fixed_Z//;
   #my $final_decomp = $fixed_Z_decomp;
   #$final_decomp =~ s/\.fixed_Z//;
    my $final_simple = $simple_vcf_name;
    $final_simple =~ s/.raw.vcf.uncleaned/.raw.vcf/;
    my $final_decomp = $decomp_vcf_name;
    $final_decomp =~ s/.decomp.vcf.uncleaned/.decomp.vcf/;
    
    if ($ref_fasta eq "unspecified")
    {
	## just sort the file and PV tag it
	
	my $cmd1 = "cat $simple_vcf_name   | $vcftools_dir/perl/vcf-sort | $isaac_bioinf_dir"."vcf_scripts/vcf_remove_dupes.pl --take_first --pass --filter_txt DUP_CALL  | $isaac_bioinf_dir"."vcf_scripts/vcf_remove_overlaps.pl --pass --filter_txt OVERLAPPING_SITE > $final_simple";
	print "$cmd1\n";
	my $ret1 = qx{$cmd1};
	print "$ret1\n";
	
	my $cmd2 = "cat $decomp_vcf_name   | $vcftools_dir/perl/vcf-sort | $isaac_bioinf_dir"."vcf_scripts/vcf_remove_dupes.pl --take_first --pass --filter_txt DUP_CALL   | $isaac_bioinf_dir"."vcf_scripts/vcf_remove_overlaps.pl --pass  --filter_txt OVERLAPPING_SITE > $final_decomp ";
	print "$cmd2\n";
	my $ret2 = qx{$cmd2};
	print "$ret2\n";
	
    }
    else # sort and PV tag and remove ref mismatches
    {
	
	#### RAW vcf
	my $tmp1 = $simple_vcf_name.".corrected_ref_mismatch";
	print "Switch ref/alt bases in raw vcf on those sites where we know we have placed them back to front:\n";
	
	my $cmd1 = $isaac_bioinf_dir."vcf_scripts/vcf_correct_strand.pl --keep_phased --filter_mismatches MISMAPPED_UNPLACEABLE $simple_vcf_name $ref_fasta > $tmp1";
	print "$cmd1\n";
	my $ret1 = qx{$cmd1};
	print "$ret1\n";
	
	print "Remove sites where Stampy has placed variant in wrong place\n";
	my $cmd2 = $isaac_bioinf_dir."vcf_scripts/vcf_align.pl  -p --tag PV LEFT $tmp1 $ref_fasta  | $vcftools_dir/perl/vcf-sort | $isaac_bioinf_dir"."vcf_scripts/vcf_remove_dupes.pl  --take_first --pass  | $isaac_bioinf_dir"."vcf_scripts/vcf_remove_overlaps.pl  --pass --filter_txt OVERLAPPING_SITE > $final_simple ";
	print "$cmd1\n";
	my $ret2 = qx{$cmd2};
	print "$ret2\n";
	
	### DECOMP VCF
	
	my $tmp3 = $decomp_vcf_name.".corrected_ref_mismatch";
	print "Switch ref/alt bases in decomp vcf on those sites where we know we have placed them back to front:\n";
	
	my $cmd3 = $isaac_bioinf_dir."vcf_scripts/vcf_correct_strand.pl --keep_phased --filter_mismatches MISMAPPED_UNPLACEABLE $decomp_vcf_name $ref_fasta > $tmp3 ";
	print "$cmd3\n";
	my $ret3 = qx{$cmd3};
	print "$ret3\n";
	
	
	my $cmd4 = $isaac_bioinf_dir."vcf_scripts/vcf_align.pl -p  --tag PV LEFT $tmp3 $ref_fasta  | $vcftools_dir/perl/vcf-sort | $isaac_bioinf_dir"."vcf_scripts/vcf_remove_dupes.pl --take_first --pass | $isaac_bioinf_dir"."vcf_scripts/vcf_remove_overlaps.pl  --pass --filter_txt OVERLAPPING_SITE  > $final_decomp ";
	print "$cmd4\n";
	my $ret4 = qx{$cmd4};
	print "$ret4\n";
	
    }
}## end of "if normal case"
else
{
    ## calls were based on a VCF, from which we made branches, and genotyped samples using those branches.
    
    ## first parse the callfile, and for each variant, collect the confidences and likelihoods
    my %GL_info=();
    my %COV_info=();
    my %GT_info=();
    collect_confidences_and_likelihoods($fh_calls, \%GL_info, \%COV_info, \%GT_info, $ploidy);
    print_final_vcf_by_modifying_sites_vcf($fh_simple_vcf, \%GL_info, \%COV_info, \%GT_info, $vcf_on_which_calls_based);
}


sub check_callfile_format
{
    my ($cf) = @_;
    open(CF, $cf)||die("Cannot open callfile $cf\n");
    my $flag = -1;
    while ($flag==-1)
    {
	my $line = <CF>;
	chomp $line;
	if ($line =~ /3p_flank/)
	{
	    <CF>;#seq
	    <CF>;#blank line
	    $line = <CF>;
	    if ($line =~ /median/)
	    {
		$flag=1;
	    }
	    else
	    {
		$flag=0;
	    }
	}
    }
    close(CF);
    return $flag;
}

##sites VCF has no sample columns
##in this case we KNOW the first allele is the ref allele.
sub print_final_vcf_by_modifying_sites_vcf
{
    my ($final_fh, $href_gl, $href_cov, $href_gt, $in_vcf) = @_;

    if ($callfile_is_median_format==0)
    {
	die("Do not expect to call this with the print_colour_coverages output format callfile\n");
    }
    open(IN, $in_vcf)||die();
    while (<IN>)
    {
	my $line = $_;
	chomp $line;


	if ($line =~ /^\#/)
	{
	    #do nothing - header alreafy printed header
	}
	else
	{
	    my @sp = split(/\s+/, $line);
	    my $name = $sp[2];
	    if ( (!exists $href_gt->{$name})
		 ||
		 (!exists $href_gl->{$name})
		 ||
		 (!exists $href_cov->{$name}) )
	    {
		die("Call Zam. Cannot match variant name $name\n");
	    }
	    print $final_fh join("\t", @sp[0..7]);
	    printf $final_fh "\t";
	    print $final_fh "GT:COV:GT_CONF";
	    if ($print_gls eq "yes")
	    {
		print $final_fh ":GL";
	    }
	    print $final_fh "\t";
	    my $col;
	    for ($col=0; $col<$number_of_colours; $col++)
	    {

		if ($col != $reference_colour)
		{
		    print $final_fh $href_gt->{$name}->{$col};
		    print $final_fh ":";
		    print $final_fh $href_cov->{$name}->{$col};
		    my $gtconf = get_gt_conf_from_commasep_gls($href_gl->{$name}->{$col});
		    print $final_fh ":$gtconf:";
		    print $final_fh $href_gl->{$name}->{$col};
		    if ($col<$number_of_colours-1)
		    {
			print $final_fh "\t";
		    }
		    else
		    {
			print $final_fh "\n";
		    }
		}
	    }
	}
    }
    close(IN);

}


sub get_gt_conf_from_commasep_gls
{
    my ($str) = @_;
    my @sp = split(/,/, $str);
    my @s = sort { $a <=> $b } @sp;##sorting log likelihoods. Find difference between max and one below

    my $rounded_conf = sprintf "%.2f", $s[scalar(@sp)-1] - $s[scalar(@sp)-2];
    return $rounded_conf;
}

##this assumes the first allele is the ref allele, so no switching needed
sub collect_confidences_and_likelihoods
{
    my ($call_fh, $href_gl, $href_cov, $href_gt, $pl) = @_;
    if ($pl==2)
    {
	my %local_hash_gl=();
	my %local_hash_gt=();
	my $temp_name = 0;
	my $realname="";
	while (<$call_fh>)
	{

	    #Colour/sample   GT_call llk_hom_br1 llk_het    llk_hom_br2
	    #0=REF   NO_CALL 0       0      0
	    #1       HOM2    -803.00 -300.00  -14.89
	    #2       HOM2    -923.99 -400.00  -12.99

	    my $line = $_;
	    chomp $line;

	    if ($line =~ /Colour\/sample/)
	    {
		$temp_name++;
		my $num_gls=3;

		my $ct = 0;
		while ($ct<$number_of_colours)
		{
		    $line = <$call_fh>;
		    chomp $line;
		    my @sp = split(/\s+/, $line);
		    my $gt = $sp[1];
		    my @gls = @sp;
		    shift @gls;
		    shift @gls;
		    $num_gls = scalar(@gls);
		    my $z;
		    my $colour = $sp[0];

		    if ($colour =~ /=REF/)
		    {
			$colour =~ s/=REF//;
		    }
		    
		    for ($z=0; $z<$num_gls; $z++)
		    {
			$local_hash_gl{$temp_name}{$colour}{$z}=$gls[$z];
		    }
		    if ($gt eq "HOM1")
		    {
			$local_hash_gt{$temp_name}{$colour}="0/0";
		    }
		    elsif ($gt eq "HET")
		    {
			$local_hash_gt{$temp_name}{$colour}="0/1";
		    }
		    elsif ($gt eq "HOM2")
		    {
			$local_hash_gt{$temp_name}{$colour}="1/1";
		    }
		    else
		    {
			$local_hash_gt{$temp_name}{$colour}="NO_CALL";
		    }
		    $ct++;
		}
		##now all the GT and GL info is saved. Next comes the fasta ibit
		$line = <$call_fh>;
		chomp $line;
		
		##eg >UNION_BC_k61_var_1_5p_flank
		if ($line =~ />(\S*var_\d+)_5p_flank/)
		{
		    $realname = $1;

		    ## we have stored all the GL info for all colours in our local hash
		    ##transfer to the global hash/passed-in hash
		    my $k;
		    my $p;
		    for ($k=0; $k<$number_of_colours; $k++)
		    {
			$href_gl->{$realname}->{$k}=$local_hash_gl{$temp_name}{$k}{0};
			for ($p=1; $p<$num_gls; $p++)
			{
			    $href_gl->{$realname}->{$k} = ($href_gl->{$realname}->{$k}).",";
			    $href_gl->{$realname}->{$k} = ($href_gl->{$realname}->{$k}).($local_hash_gl{$temp_name}{$k}{$p});
			}
			$href_gt->{$realname}->{$k}=$local_hash_gt{$temp_name}{$k};
		    }
		    
		    ### now we can get the coverage information
		    <$call_fh>;#5pseq
		    <$call_fh>;#br1 readid
		    <$call_fh>;#br1 seq
		    <$call_fh>;#br2 readid
		    <$call_fh>;#br2 seq
		    <$call_fh>;#3p readid
		    <$call_fh>;#3p seq
		    <$call_fh>;#blank line
		    $line = <$call_fh>;
		    if ($line !~ /br1_median_covg/)
		    {
			die("Parsing problem with this callfile - I get $line, but I expect to see \"Colour  br1_median_covg br2_median_covg\"");
		    }
		    
		    for ($k=0; $k<$number_of_colours; $k++)
		    {
			$line = <$call_fh>;
			chomp $line;
			my @ar = split(/\t/, $line);
			my $covg1 = $ar[1];
			my $covg2 = $ar[2];
			$href_cov->{$realname}->{$k}="$covg1,$covg2";
			if ( ($covg1==0) && ($covg2==0) )
			{
			    $href_gt->{$realname}->{$k}="./.";
			}
		    }
		    
		}
		else
		{
		    die("Parsing problem - expected to see 5p flank but got $line\n");
		}
	    }
	}
	close($call_fh);
    }

}
sub get_covg_and_total_seq_from_log
{
    my ($log, $href_rlen, $href_seq) = @_;

    open(LOG, $log)||die("Unable to open caller log file $log\n");

    my $outer_done=0;
    while ($outer_done==0)
    {
	my $ln = <LOG>;
	chomp $ln;
	if ($ln =~ /SUMMARY/)
	{
	    <LOG>;
	    my $done = 0;
	    while ($done==0)
	    {
		$ln = <LOG>;
		chomp $ln;
		if ($ln !~ /\*/)
		{
		    my @sp = split(/\t/, $ln);
		    my $colour = $sp[0];
		    my $readlen = $sp[1];
		    my $seq = $sp[2];
		    $href_rlen->{$colour}=$readlen;
		    $href_seq->{$colour} = $seq;
		}
		else
		{
		    $done=1;
		    $outer_done=1;
		}
	    }

	}
	
    }
    close(LOG);
}

sub get_extra_info_for_calls
{
    my ($href, $file) = @_;

    open(UN, $file)||die();
    while (<UN>)
    {
	my $line = $_;
	chomp $line;


	if ( $line =~ /(\w*var_\d+)_5p_flank/ )
	{
	    my $varname = $1;
	    if ($prefix ne "")
	    {
		$varname = $prefix."_".$varname;
	    }
	    if ( $varname eq "" )
	    {
		die(
		    "get_extra_info_for_calls - Found this line :$line in the unionfile $file  - expected it to be of the form \\w+var_<NUMBER>_5p_flank\n"
		    );
	    }

	    if ($line =~ /\w*var_\d+_5p_flank\s+INFO:(\S+)/)
	    {
		my $extra_info=$1;
		$href->{$varname} = $extra_info;
	    }

	}
    }
    close(UN);

}


sub fix_Z
{
    my ($oldvcf, $outvcfname, $href) = @_;
    open(OLD, $oldvcf)||die("Cannot open $oldvcf");
    open(NEW, ">".$outvcfname)||die("Cannot open $outvcfname");

    while(<OLD>)
    {
	my $bobble = $_;
	chomp $bobble;
	if ($bobble =~ /^\#/)
	{
	    print NEW "$bobble\n";
	}
	else
	{
	    my @sp = split(/\t/, $bobble);
	    my $chr = $sp[0];
	    my $pos = $sp[1];
	    my $ref = $sp[3];
	    my $alt = $sp[4];
	    if (scalar(@sp)<4)
	    {
		die("Less than 4 fields in $bobble\n");
	    }

	    if (!defined($ref))
	    {
		print "ZAM! ref is not define on this line:\n$bobble\nof vcf $oldvcf\n";
	    }

	    if ($ref =~ /^Z/)
	    {
		if (!exists $href->{$chr."_".$pos})
		{
		    die("Coding error Zam - $chr $pos is not in your hash for Z replacement. The VCf line is\n$bobble\n");
		}
		my $fix = $href->{$chr."_".$pos};
		$sp[3] =~ s/Z/$fix/;
		$sp[4] =~ s/Z/$fix/;
		print NEW join("\t", @sp);
		print NEW "\n";

	    }
	}
    }
    close(OLD);
    close(NEW);
}
sub get_characters_to_replace_Z
{
    my ($href, $reference_fa) = @_;
    if ($reference_fa eq "unspecified")
    {
	foreach my $k (keys %$href)
	{
	    $href->{$k} = -1;
	}
    }
    else
    {
	my %localhash=();
	foreach my $k (keys %$href)
        {
	    my $c; my $p;
	    if ($k =~ /(\S+)_(\S+)/)
	    {
		$c = $1;
		$p = $2;
		$localhash{$c}{$p}=1;
	    }
	    else
	    {
		die("Problem parsing chromosome name in my hash? k is $k, does not have an underscore. call zam");
	    }
	}
	## read the ref fasta once and collect bases
	open(REFFO, $reference_fa)||die("Cannot open $reference_fa");
	my $curr_chro = "-1";
	my $curr_pos=-1;
	while(<REFFO>)
	{
	    my $lion = $_;
	    chomp $lion;
	    if ($lion =~ /^>(\S+)/)
	    {
		$curr_chro=$1;
		$curr_pos=0;
	    }
	    else
	    {
		my @sp = split(//, $lion);
		my $i;
		for ($i=0; $i<scalar(@sp); $i++)
		{
		    $curr_pos++;
		    if (exists $localhash{$curr_chro}{$curr_pos})
		    {
			$href->{$curr_chro."_".$curr_pos}=$sp[$i];
		    }
		}
	    }
	}
	close(REFFO);
    }
}

sub get_Z_positions
{
    my ($href, $vcf) = @_;
    open(V, $vcf)||die();
    while (<V>)
    {
	my $ly = $_;
	chomp $ly;
	if ($ly !~ /^\#/)
	{
	    my @sp = split(/\t/, $ly);
	    my $chr = $sp[0];
	    my $pos = $sp[1];
	    my $ref = $sp[3];
	    my $alt = $sp[4];
	    if (!defined($ref))
	    {
		print "ZAM! ref is not define on this line:\n$ly\nof vcf $vcf\n";
	    }
	    if (scalar(@sp)<4)
	    {
		die("Less than 4 fields in $ly\n");
	    }
	    if ($ref =~ /^Z/)
	    {
		$href->{$chr."_".$pos}=1;
	    }
	}
    }
    close(V);
}
sub check_if_callfile_in_legacy_format
{
    my ($file) = @_;
    open (CHECK, $file)||die("Cannot open $file");
    my $done = -1;
    while ($done ==-1)
    {
	my $ln = <CHECK>;
	if ($ln =~ /branch/)
	{
	    if ($ln =~ /branch\_\d+\_(1|2)/)
	    {
		$done=1;#legacy
	    }
	    else
	    {
		$done = 0;#not legacy
	    }
	}
    }
    close(CHECK);
    return $done;
}
sub get_vcf_header
{
	my ($colourfile) = @_;

	my $check_file_cmd = "wc -l $colourfile";
	my $check_ret = qx{$check_file_cmd};
	if ($check_ret =~ /^(\d+)\s/)
	{
	    my $nc = $1;
	    if ($nc==$number_of_colours-1)
	    {
		die("Sample list file has one less line ($nc) than there are colours ($number_of_colours). Is it missing a carriage return on the final line? Or is there one too few line? Please fix and rerun");
	    }
	    elsif ($nc==$number_of_colours)
	    {
		#fine
	    }
	    else
	    {
		die("You must pass in a samplename_list file which contains one line per colour, each line ending with a carriage return (ie hit enter at the end of each line, or use \\n if you are printing with perl");
	    }
	}

	my $date_cmd = "date \'+\%d\/\%m\/\%y\'";
	my $date     = qx{$date_cmd};
	chomp $date;

	my $head = "";

	$head = $head . "##fileformat=VCFv4.0\n";
	$head = $head . "##fileDate=$date\n";
	$head = $head
	  . "##phasing=none, though some calls involve phasing clustered variants\n";
	$head = $head
	  . "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	$head = $head
	  . "##FORMAT=<ID=COV,Number=2,Type=Integer,Description=\"Number of reads on ref and alt alleles\">\n";
	$head = $head
	  . "##FORMAT=<ID=GT_CONF,Number=1,Type=Float,Description=\"Genotype confidence. Difference in log likelihood of most likely and next most likely genotype\">\n";
	$head = $head
	  . "##FORMAT=<ID=SITE_CONF,Number=1,Type=Float,Description=\"Probabilitic site classification confidence. Difference in log likelihood of most likely and next most likely model (models are variant, repeat and error)\">\n";
	if ($print_gls eq "yes")
	{
	    $head = $head."##FORMAT=<ID=GL,Number=.,Type=Float,Description=\"Genotype Likelihood. Log likelihoods of all possible genotypes\">\n";
	}
	$head = $head
	  . "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n";
	$head = $head
	  . "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of variant\">\n";
	$head = $head . "##ALT=<ID=SNP,Description=\"SNP\">\n";
	$head = $head
	  . "##ALT=<ID=SNP_FROM_COMPLEX,Description=\"SNP called from a cluster of phased SNPs or complex SNP/indel , split out for easier comparison with other SNP call sets\">\n";
	$head = $head . "##ALT=<ID=DEL,Description=\"Deletion\">\n";
	$head =
	  $head . "##ALT=<ID=INS,Description=\"Insertion of novel sequence\">\n";
	$head = $head . "##ALT=<ID=INDEL,Description=\"Insertion-deletion\">\n";
	$head = $head . "##ALT=<ID=INV,Description=\"Inversion\">\n";
	$head = $head
	  . "##ALT=<ID=INV_INDEL,Description=\"Inversion+indel - this script overcalls these, so worth checking\">\n";
	$head = $head . "##ALT=<ID=DEL_INV,Description=\"Deletion + Inversion\">\n";
	$head =
	  $head . "##ALT=<ID=INS_INV,Description=\"Insertion + Inversion\">\n";
	$head = $head . "##ALT=<ID=PH_SNPS,Description=\"Phased SNPs\">\n";
	$head = $head
	  . "##ALT=<ID=COMPLEX,Description=\"Complex variant, collection of SNPs and indels\">\n";
	$head = $head
	  . "##FILTER=<ID=MAPQ,Description=\"5prime flank maps to reference with mapping quality below $mapping_qual_thresh\">\n";
	$head = $head
	  . "##FILTER=<ID=MISMAPPED_UNPLACEABLE,Description=\"Stampy mapped the variant (using the 5p-flank) confidently (mapqual> $mapping_qual_thresh) to a place where the ref-allele does not match\">\n";
	$head = $head
	  . "##FILTER=<ID=MULTIALLELIC, Description=\"Cortex does not call multiallelic sites, but combining run_calls VCFs can produce them. Filtered as current genotyper assumes biallelic.\">\n";
	$head = $head
	  . "##FILTER=<ID=OVERLAPPING_SITE, Description=\"If Stampy (or combining VCFs) has placed two biallelic variants overlapping, they are filtered\">\n";

	$head = $head . "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
	#if ( $colourfile !~ /,/ )
	#{
	open( COLOURS, $colourfile ) || die("Cannot open $colourfile");
	my $z;
	for ( $z = 0 ; $z < $number_of_colours ; $z++ )
	{
	    my $line = <COLOURS>;
	    chomp $line;
	    
	    if ( ( $z == $reference_colour ) || ( $z == $pooled_colour ) )
	    {
		next;
	    }
	    
	    $head = $head ."$line";
	    if ( $z == $last_sample_col )
	    {
		$head = $head . "\n";
	    }
	    else
	    {
		$head = $head . "\t";
	    }
	}
	close(COLOURS);

	#else
	#{
	#	my @colors = split( /,/, $colourfile );
	#	@colors = grep { $_ ne "REF" } @colors;
	#	$head .= join( "\t", @colors ) . "\n";

	#}
	return $head;
}

sub  filter_calls_where_entire_ref_allele_is_multicopy
{
    my ($href, $callfile, $ref_colour) = @_;
    
    my $fh;
    open($fh, $callfile)||die("Cannot open the callfile $callfile\n");

    my $line = "";
    my $varname;
    my $flank5p;
    my $br1;
    my $br2;
    my $flank3p;

    my $classification = "VARIANT";
    my $class_llk_rep  = -10;
    my $class_llk_var  = -1;
    

    my @arr_br1_covgs = (); 
    my @arr_br2_covgs = (); 

    while (( $line !~ /(\w*var_\d+)_5p_flank/ )
	   && ( $line !~ /PASSES/ )
	   && ( $line !~ /FAILS/ )
	   && ( $line !~ /Colour\/sample/ ) )
    {
	if ( eof($fh) )
	{
 	    close($fh);
	    return;
	}
	else
	{
	    $line = <$fh>;
	}
    }
    
    if ( $line =~ /PASSES/ )
    {
	$classification = "VARIANT";
	if ( eof($fh) )
	{
	    close($fh);
	    return;
	}
	
	$line = <$fh>;
	chomp $line;
	
	if ( $line =~ /llk_var:(\S+).+llk_rep:(\S+)/ )
	{
	    $class_llk_rep = $1;
	    $class_llk_var = $2;
	}
	else
	{
	    die(
		"Parsing issue. Found model selecion result but cant find likelihoods on $line\n"
		);
	}
	if ( eof($fh) )
	{
	    close($fh);
	    return ;
	}
	$line = <$fh>;
	
    }
    elsif ( $line =~ /FAILS/ )
    {
	$classification = "REPEAT";
	if ( eof($fh) )
	{
	    close($fh);
	    return;
	}
	
	$line = <$fh>;
	chomp $line;
	
	if ( $line =~ /llk_var:(\S+).+llk_rep:(\S+)/ )
	{
	    $class_llk_rep = $1;
	    $class_llk_var = $2;
	}
	else
	{
	    die(
		"Parsing issue. Found model selecion result but cant find likelihoods on $line\n"
		);
	}
	if ( eof($fh) )
	{
	    close($fh);
	    return;
	}
	$line = <$fh>;
    }
    
    if ( $line =~ /Colour\/sample/ )
    {
	my $j;
	for ( $j = 0 ; $j < $number_of_colours ; $j++ )
	{
	    $line = <$fh>;#### all the genotype calls
	    chomp $line;
	}
	$line = <$fh>;
	
    }
    
    if ( $line =~ /(\w*var_\d+)_5p_flank/ )
    {
	$varname = $1;
	if ($prefix ne "")
	{
	    $varname = $prefix."_".$varname;
	}
	
	if ( $varname eq "" )
	{
	    die(
		"Found this line :$line in the callfile - expected it to be of the form var_<NUMBER>_5p_flank\n"
		);
	}
	
	
	$flank5p = <$fh>;
	chomp $flank5p;
	<$fh>;    #ignore br1 read id
	$br1 = <$fh>;
	chomp $br1;
	<$fh>;    #ignore br2 read id
	$br2 = <$fh>;
	chomp $br2;
	<$fh>;    #ignore 3p flank read id
	$flank3p = <$fh>;
	chomp $flank3p;
	$line    = <$fh>;
	
	if ( $line =~ /extra information/ )
	{
	    <$fh>;
	}
	$line = <$fh>;
	chomp $line;
	if ( $line !~ /branch1 coverages/ )
	{
	    die("Expected to see \"branch 1 coverages\" but instead saw $line");
	}
	my $z;
	for ( $z = 0 ; $z < $number_of_colours ; $z++ )
	{
	    $line = <$fh>;
	    chomp $line;
	    if ( $line !~ /Covg in Colour/ )
	    {
		die("Expected to see \"Covg in Colour\" but instead saw $line");
	    }
	    $line = <$fh>;
	    chomp $line;
	    if ($z==$ref_colour)
	    {
		my @br1 = split( /\s+/, $line );
		
		my $min1 =  get_min_covg(\@br1);
		if ($min1>1)
		{
		    ### Every kmer on this allele occurs >1 time in reference. We need to filter this call
		    $href->{$varname}=1;
		}
	    }
	}
	$line = <$fh>;
	chomp $line;
	if ( $line !~ /branch2 coverages/ )
	{
	    die("Expected to see \"branch 2 coverages\" but instead saw $line");
	}
	
	for ( $z = 0 ; $z < $number_of_colours ; $z++ )
	{
	    $line = <$fh>;
	    chomp $line;
	    
	    if ( $line !~ /Covg in Colour/ )
	    {
		die("Expected to see \"Covg in Colour\" but instead saw $line");
	    }
	    $line = <$fh>;
	    ### branch 2 is the alt allele
	}
	
	close($fh);
	return;
	
    }
    else
    {
	print
	    "Unexpected error on  $line. Expected the read-id of 5prime flank in callfile to be of the form (optional_text)var_(number)_5p_flank. This bug is probably due to Zam's changes to support the new format of callfiles but also stay legacy format compatible - contact zam\@well.ox.ac.uk\n";
	die();
    }
    close($fh);
    return;
}

sub filter_by_flank_mapqual
{
	my ( $file, $href ) = @_;

	open( FILE, $file ) || die();

	while (<FILE>)
	{
		my $line = $_;

		if ( $line =~ /^@/ )
		{
		}
		else
		{
			my @sp = split( /\t/, $line );

			if ( scalar @sp < 4 )
			{
				die("problem parsing $line");
			}

			my $query = $sp[0];
			my $varname;
			if ( $query =~ /(\w*var_\d+)_5p_flank/ )
			{
			    $varname = $1;
			    if ($prefix ne "")
			    {
				$varname = $prefix."_".$varname;
			    }


			}
			else
			{
				die(
				    "Expected query name in sam file would be of form var_1_5p_flank, but is $query"
				);
			}
			if  ( ( $sp[4] >= $mapping_qual_thresh )
			      &&
			      ( $sp[5] =~ /^\d+M$/) )
			    #  &&
			    #  ( $sp[9] =~ /^[=]+$/) )
			{
				$href->{$varname} = "PASS";
			}
			else
			{
				$href->{$varname} = "FAIL";
				if (
				    ($sp[4] >= $mapping_qual_thresh)
				    &&
				    (
				      ( $sp[5] !~ /^\d+M$/)
				      ||
				      ( $sp[9] !~ /^[=]+$/) 
				     ) 
				    )
				{
				    print "$varname passes mq filter, but flank is contains insertions/deletions, so filter this out (would get wrong coord)\n";
				}
			}
		}
	}

	close(FILE);

}

## Call this as pass through the callfile and proc_bubble file and mapped flanks file, will work out the
## next entries for vcf\s and print them out
## Only call this for
sub print_next_vcf_entry_for_easy_and_decomposed_vcfs
{

	my (
		$file_handle_calls,              $file_handle_map_flanks,
		$file_handle_proc_bubble_output, $href_varname_to_passing_all_filters,
		$href_var_name_to_cut_flank,     $print_easy_vcf,
		$print_decomp_vcf,               $fh_easy,
		$fh_decomp,                      $href_pop_classifier_confidence,
	        $colour_ref
	) = @_;

	##1. Get next var info from all three files
	my (
		$eof,          $var_name,       $flank5p,       $br1_seq,
		$br2_seq,      $flank3p,        $aref_br1_cov,  $aref_br2_cov,
		$which_is_ref, $classification, $class_llk_rep, $class_llk_var,
		$genotype,     $llk_hom1,       $llk_het,       $llk_hom2, $extra_info_fields
	) = get_next_var_from_callfile( $file_handle_calls, $ploidy );


	my ( $eof2, $var_name2, $strand, $chr, $coord );
	if ($colour_ref !=-1)
	{
	    ## then we have already mapped the flanks to the reference:
	    ( $eof2, $var_name2, $strand, $chr, $coord )=
		get_next_var_from_flank_mapfile($file_handle_map_flanks);
	}
	else
	{
	    ( $eof2, $var_name2, $strand, $chr, $coord )
		= ($eof, $var_name, 0, $global_var_ctr, 1000000);#give it fake position of 1million, so can add/subtract flank and still be >0.
	    $global_var_ctr++;
	}

	my (
		$eof3,
		$var_name3,
		$align_num_bases_agreement_at_start,
		$align_num_bases_agreement_at_end,
		$align_direction,
		$num_snps_align,
		$num_indels_align,
		$alignment_br1,
		$alignment_middle,
		$alignment_br2,
		$possible_inversion,
		$clean_indel,
		$aref_snp_coords,
		$aref_snp_alleles,
		$aref_coords_of_indels,
		$aref_alleles_of_indels,
		$aref_indel_needs_extra_base
	  )

	  = get_next_alignment_from_procfile( $file_handle_proc_bubble_output,
		$strand, $which_is_ref );

	## if one but not all reach the end of file
	if ( ( ( $eof eq "EOF" ) || ( $eof2 eq "EOF" ) || ( $eof3 eq "EOF" ) )
		&& (
			!( ( $eof eq "EOF" ) && ( $eof2 eq "EOF" ) && ( $eof3 eq "EOF" ) ) )
	  )
	{
		print join( "\t", "eof",  ( $eof  eq "EOF" ), $callfile ),        "\n";
		print join( "\t", "eof2", ( $eof2 eq "EOF" ), $mapped_flanks ),   "\n";
		print join( "\t", "eof3", ( $eof3 eq "EOF" ), $proc_bub_output ), "\n";

		die("SOme of these files end before the others");
	}
	##if all hit the end of file
	if ( ( $eof eq "EOF" ) && ( $eof2 eq "EOF" ) && ( $eof3 eq "EOF" ) )
	{
		return 0;
	}


	if ( ( $which_is_ref eq "neither" ) || ( $which_is_ref eq "b" ) )
	{
	    $which_is_ref = 1;# this will be caught by post-processing scripts.
	}
	

	if ( ( $var_name ne $var_name2 ) || ( $var_name ne $var_name3 ) )
	{
		die(
"Some ordering problem, next var in the three files are $var_name, $var_name2, $var_name3\n\n"
			  . join( "\n", $callfile, $mapped_flanks, $proc_bub_output ) );
	}

	my $filter_result;
	if ( exists $href_varname_to_passing_all_filters->{$var_name} )
	{
		$filter_result = $href_varname_to_passing_all_filters->{$var_name};
	}
	else
	{
		$filter_result = "PASS";
	}

	if ( $print_easy_vcf == 1 )
	{
		my $split_phased_snps = 0;
		my @empty             = ();
		print_vcf_entry(
			$filter_result,
			$fh_easy,
			$var_name,
			$flank5p,
			$br1_seq,
			$br2_seq,
			$flank3p,
			$strand,
			$chr,
			$coord,
			$which_is_ref,
			$align_num_bases_agreement_at_start,
			$align_num_bases_agreement_at_end,
			$align_direction,
			$num_snps_align,
			$num_indels_align,
			$alignment_br1,
			$alignment_middle,
			$alignment_br2,
			$aref_br1_cov,
			$aref_br2_cov,
			$href_var_name_to_cut_flank,
			$possible_inversion,
			$clean_indel,
			$split_phased_snps,
			0,
			0,
			0,
			0,
			\@empty,
			$classification,
			$class_llk_rep,
			$class_llk_var,
			$genotype,
			$llk_hom1,
			$llk_het,
			$llk_hom2,
			$href_pop_classifier_confidence,
		        $extra_info_fields
		);
	}
	if ( $print_decomp_vcf == 1 )
	{
		my $split_phased_snps = 1;
		print_vcf_entry(
			$filter_result,
			$fh_decomp,
			$var_name,
			$flank5p,
			$br1_seq,
			$br2_seq,
			$flank3p,
			$strand,
			$chr,
			$coord,
			$which_is_ref,
			$align_num_bases_agreement_at_start,
			$align_num_bases_agreement_at_end,
			$align_direction,
			$num_snps_align,
			$num_indels_align,
			$alignment_br1,
			$alignment_middle,
			$alignment_br2,
			$aref_br1_cov,
			$aref_br2_cov,
			$href_var_name_to_cut_flank,
			$possible_inversion,
			$clean_indel,
			$split_phased_snps,
			$aref_snp_coords,
			$aref_snp_alleles,
			$aref_coords_of_indels,
			$aref_alleles_of_indels,
			$aref_indel_needs_extra_base,
			$classification,
			$class_llk_rep,
			$class_llk_var,
			$genotype,
			$llk_hom1,
			$llk_het,
			$llk_hom2,
			$href_pop_classifier_confidence, $extra_info_fields
		);

	}
	return 1;

}

sub print_vcf_entry
{
	my (
		$filter_result,
		$fh_output_vcf,
		$var_name,
		$flank5p,
		$br1_seq,
		$br2_seq,
		$flank3p,
		$strand,
		$chr,
		$coord,
		$which_is_ref,
		$align_num_bases_agreement_at_start,
		$align_num_bases_agreement_at_end,
		$align_direction,
		$num_snps_align,
		$num_indels_align,
		$alignment_br1,
		$alignment_middle,
		$alignment_br2,
		$aref_br1_cov,
		$aref_br2_cov,
		$href_var_name_to_cut_flank,
		$poss_inversion,
		$clean_indel,
		$split_phased_snps,
		$aref_snp_coords,
		$aref_snp_alleles,
		$aref_indel_coords,
		$aref_indel_alleles,
		$aref_indel_needs_extra_base,
		$classification,
		$class_llk_rep,
		$class_llk_var,
		$genotype,
		$llk_hom1,
		$llk_het,
		$llk_hom2,
		$href_pop_classifier_conf,
	        $extra_info
	) = @_;

	my $vcf_entry_chr = $chr;
	my $vcf_entry_pos;
	my $vcf_entry_ref_allele;
	my $vcf_entry_alt_allele;
	my $error;

	if ( $classification eq "REPEAT" )
	{
		$filter_result = $filter_result . "CLASSIF_REPEAT";
	}

	( $vcf_entry_pos, $vcf_entry_ref_allele, $vcf_entry_alt_allele, $error ) =
	  get_simple_vcf_entry_pos_and_alleles(
		$var_name,
		$strand,
		$coord,
		$flank5p,
		$br1_seq,
		$br2_seq,
		$flank3p,
		$which_is_ref,
		$align_num_bases_agreement_at_start,
		$align_num_bases_agreement_at_end,
		$align_direction,
		$num_snps_align,
		$num_indels_align,
		$href_var_name_to_cut_flank
	  );

	if ( $error ne "0" )
	{

		print ("Ignore this $var_name - due to this error $error\n");
		return;
	}

	## quick check

	my $count_num_indels_needing_extra_base = 0;
	foreach my $in (@$aref_indel_needs_extra_base)
	{
		if ( $in == 1 )
		{
			$count_num_indels_needing_extra_base++;
		}
	}
	if ( $count_num_indels_needing_extra_base > 1 )
	{
		die(
"More than one indel in $var_name nees an extra base - namely $count_num_indels_needing_extra_base"
		);
	}

	my $svlen = length($vcf_entry_alt_allele) - length($vcf_entry_ref_allele);
	my $svtype;

	if ( ( $num_snps_align == 1 ) && ( $num_indels_align == 0 ) )
	{
		$svtype = "SNP";
	}
	elsif ( $poss_inversion == 1 )
	{
		$svtype = "INV";
	}
	elsif ( $poss_inversion == 2 )
	{
		$svtype = "INV_INDEL";
	}
	elsif (( $num_snps_align > 1 )
		&& ( $num_indels_align == 0 )
	  )    ## this can only happen if the two branches are the same length :-)
	{
		$svtype = "PH_SNPS";
	}
	elsif ( ( $clean_indel == 1 ) && ( $svlen > 0 ) )
	{
		$svtype = "INS";
	}
	elsif ( ( $clean_indel == 1 ) && ( $svlen < 0 ) )
	{
		$svtype = "DEL";
	}
	elsif ( ( $num_snps_align == 0 ) && ( $num_indels_align == 1 ) )
	{
		$svtype = "INDEL";
	}
	else
	{

		#$svtype = "PH_SNPS_INDELS";
		$svtype = "COMPLEX";
	}

	my $info = "SVTYPE=$svtype;SVLEN=$svlen";
	if ($extra_info ne "")
	{
	    $info = $info.";$extra_info";
	}
	
	if (   ( ( $svtype !~ /COMPLEX/ ) && ( $svtype !~ /PH_SNPS/ ) )
		|| ( $split_phased_snps == 0 ) )
	{
		if ( $vcf_entry_ref_allele =~ /([^ACGT]+)([ACGT]+)/ )
		{
			my $new = $2;

#print "Removed $1 from $var_name 's sequence. Used to be $vcf_entry_ref_allele and now is $new\n";
			$vcf_entry_ref_allele = $new;
		}
		if ( $vcf_entry_ref_allele =~ /([ACGT]+)([^ACGT]+)/ )
		{
			my $new = $1;

#print "Removed $2 from $var_name 's sequence. Used to be $vcf_entry_ref_allele and now is $new\n";
			$vcf_entry_ref_allele = $new;
		}

		print $fh_output_vcf
"$vcf_entry_chr\t$vcf_entry_pos\t$var_name\t$vcf_entry_ref_allele\t$vcf_entry_alt_allele\t.\t$filter_result\t$info\t";

		my $have_called_gt       = 0;
		my $have_used_classifier = 0;
		if (   ( scalar(@$genotype) > 0 )
			&& ( scalar( keys %$href_pop_classifier_conf ) > 0 ) )
		{
			print $fh_output_vcf "GT:COV:GT_CONF:SITE_CONF\t";
			$have_called_gt       = 1;
			$have_used_classifier = 1;
		}
		elsif (
			scalar(@$genotype) >
			0 )   ## genotypes called, but no site confidences/no pop classifier
		{
			print $fh_output_vcf "GT:COV:GT_CONF\t";
			$have_called_gt = 1;
		}
		elsif (
			scalar( keys %$href_pop_classifier_conf ) > 0
		  )   ##pop classifier, but no genotypes called (seems impossible to me)
		{
			print $fh_output_vcf "GT:COV:SITE_CONF\t";
			$have_used_classifier = 1;
		}
		else
		{
			print $fh_output_vcf "COV\t";
		}

		print_all_genotypes_and_covgs(
			$fh_output_vcf,  $aref_br1_cov,
			$aref_br2_cov,   $which_is_ref,
			$classification, $class_llk_rep,
			$class_llk_var,  $genotype,
			$llk_hom1,       $llk_het,
			$llk_hom2,       $href_pop_classifier_conf,
			$var_name
		);

	}
	else
	{

		## double check
		if ( scalar(@$aref_snp_coords) != scalar(@$aref_snp_alleles) )
		{
			print("SNP coord and Allele arrays of different lengths: ");
			print scalar(@$aref_snp_coords);
			print " ";
			print scalar(@$aref_snp_alleles);
			die("\n");
		}

		## we have coordinates of phased SNPs with respect to fw_middle.
		## we will inherit the $genotype and $cov data from the overall variant.
		## So all that changes are: POS, NAME=var_name_snp_NUMBER REF, ALT, sv_type=SNP_FROM_PHASED_SNP_CALL, svlen=0,genotype,
		my $cnt = 0;
		foreach my $c (@$aref_snp_coords)
		{
			$cnt++;
			my $two_alleles = $aref_snp_alleles->[ $cnt - 1 ];
			my $this_snp_ref_allele;
			my $this_snp_alt_allele;
			if ( $two_alleles =~ /^(\S)_(\S)$/ )
			{
				$this_snp_ref_allele = $1;
				$this_snp_alt_allele = $2;
			}
			else
			{
				die("Bad formatting for snp alleles: $two_alleles");
			}

#we have already got the snp coords with respect to the start of the overall-ref allele, AND the ref/alt alleles for each snp
#in the right direction
			my $this_snp_chr = $vcf_entry_chr;

			## take care here. $vcf_entry_pos is the right coordinate if we are going to call an indel - ie one base BEFORE the first variant.
			## if we want to print SNPs, then bloody vcf wants them to be the SAME base as the SNP, not th base before.
			my $this_snp_pos = $vcf_entry_pos + 1 + $c;

			my $this_snp_name = $var_name . "_sub_snp_" . $cnt;
			my $this_snp_info = "SVTYPE=SNP_FROM_COMPLEX;SVLEN=0";
			if ($extra_info ne "")
			{
			    $this_snp_info = $this_snp_info.";$extra_info";
			}

			print $fh_output_vcf
"$this_snp_chr\t$this_snp_pos\t$this_snp_name\t$this_snp_ref_allele\t$this_snp_alt_allele\t.\t$filter_result\t$this_snp_info\t";
			my $have_called_gt       = 0;
			my $have_used_classifier = 0;

			if (   ( scalar(@$genotype) > 0 )
				&& ( scalar( keys %$href_pop_classifier_conf ) > 0 ) )
			{
				print $fh_output_vcf "GT:COV:GT_CONF:SITE_CONF\t";
				$have_called_gt       = 1;
				$have_used_classifier = 1;

			}
			elsif ( scalar(@$genotype) > 0 )   ##just genotype and no pop filter
			{
				print $fh_output_vcf "GT:COV:GT_CONF\t";
				$have_called_gt = 1;

			}
			elsif ( scalar( keys %$href_pop_classifier_conf ) > 0 )
			{
				print $fh_output_vcf "GT:COV:SITE_CONF\t";
				$have_used_classifier = 1;
			}
			else
			{
				print $fh_output_vcf "COV\t";
			}
			print_all_genotypes_and_covgs(
				$fh_output_vcf,  $aref_br1_cov,
				$aref_br2_cov,   $which_is_ref,
				$classification, $class_llk_rep,
				$class_llk_var,  $genotype,
				$llk_hom1,       $llk_het,
				$llk_hom2,       $href_pop_classifier_conf,
				$var_name
			);

		}
		$cnt = 0;
		my $index_of_indel = 0;
		foreach my $c (@$aref_indel_coords)
		{
			$cnt++;
			my $needs_extra_base =
			  $aref_indel_needs_extra_base->[$index_of_indel];
			my $two_alleles = $aref_indel_alleles->[ $cnt - 1 ];
			my $this_indel_ref_allele;
			my $this_indel_alt_allele;
			if ( $two_alleles =~ /^([ACGTN]+)_([ACGTN]+)$/ )
			{
				$this_indel_ref_allele = $1;
				$this_indel_alt_allele = $2;
			}
			elsif ( $two_alleles =~ /^_([ACGTN]+)$/ )
			{
				$this_indel_ref_allele = "";
				$this_indel_alt_allele = $1;
			}
			elsif ( $two_alleles =~ /^([ACGTN]+)_$/ )
			{
				$this_indel_ref_allele = $1;
				$this_indel_alt_allele = "";
			}
			else
			{
				die("Bad formatting for indel alleles: $two_alleles");
			}

			if ( $needs_extra_base == 1 )
			{
				$this_indel_ref_allele = substr( $vcf_entry_ref_allele, 0, 1 )
				  . $this_indel_ref_allele;
				$this_indel_alt_allele = substr( $vcf_entry_alt_allele, 0, 1 )
				  . $this_indel_alt_allele;
			}

#we have already got the snp coords with respect to the start of the overall-ref allele, AND the ref/alt alleles for each snp
#in the right direction

			my $this_indel_chr = $vcf_entry_chr;

			## take care here. $vcf_entry_pos is the right coordinate if we are going to call an indel for the whole shebang - ie one base BEFORE the first variant.
			## hence add 1, as $c has already also subtracted 1
			my $this_indel_pos = $vcf_entry_pos + $c + 1;

			my $this_indel_name = $var_name . "_sub_indel_" . $cnt;
			my $svlen =
			  length($this_indel_ref_allele) - length($this_indel_alt_allele);
			my $this_indel_info = "SVTYPE=INDEL_FROM_COMPLEX;SVLEN=$svlen";
			if ($extra_info ne "")
			{
			    $this_indel_info = $this_indel_info.";$extra_info";
			}

			print $fh_output_vcf
"$this_indel_chr\t$this_indel_pos\t$this_indel_name\t$this_indel_ref_allele\t$this_indel_alt_allele\t.\t$filter_result\t$this_indel_info\t";

			my $have_called_gt       = 0;
			my $have_used_classifier = 0;
			if (   ( scalar @$genotype > 0 )
				&& ( scalar( keys %$href_pop_classifier_conf ) > 0 ) )
			{
				print $fh_output_vcf "GT:COV:GT_CONF:SITE_CONF\t";
				$have_called_gt       = 1;
				$have_used_classifier = 1;
			}
			elsif ( scalar @$genotype > 0 )
			{
				print $fh_output_vcf "GT:COV:GT_CONF\t";
				$have_called_gt = 1;
			}
			elsif ( scalar( keys %$href_pop_classifier_conf ) > 0 )
			{
				print $fh_output_vcf "GT:COV:SITE_CONF\t";
				$have_used_classifier = 1;
			}
			else
			{
				print $fh_output_vcf "COV\t";
			}

			print_all_genotypes_and_covgs(
				$fh_output_vcf,  $aref_br1_cov,
				$aref_br2_cov,   $which_is_ref,
				$classification, $class_llk_rep,
				$class_llk_var,  $genotype,
				$llk_hom1,       $llk_het,
				$llk_hom2,       $href_pop_classifier_conf,
				$var_name
			);

			$index_of_indel++;
		}
	}

}

sub log10
{
	my ($n) = (@_);
	my $tmp = log($n) / log(10);
	if ( $tmp < -1000000 )
	{
		$tmp = -1000000;
	}
	return $tmp;
}

sub log10_factorial
{
	my ($n) = @_;
	my $total = 0;
	my $i;
	for ( $i = 1 ; $i <= $n ; $i++ )
	{
		$total += log10($i);
	}
	return $total;
}

## returns for each colour: log prob (b1b1|data), log prob(b1b2|data), log prob(b2b2|data), genotype (max likelihood)
sub new_genotyper
{
	## each array has one entry per colour - reember to ignore the ref and pool colours
	my ( $aref_genotype, $aref_llk_hom1, $aref_llk_het, $aref_llk_hom2 ) = @_;

	my $i;
	my @big_results = ();
	for ( $i = 0 ; $i < $number_of_colours ; $i++ )
	{
		my @colour_results = ();
		if ( ( $i == $reference_colour ) || ( $i == $pooled_colour ) )
		{
			push @colour_results, 0;
			push @colour_results, 0;
			push @colour_results, 0;
			push @colour_results, "IGNORE";

		}
		else
		{
			if ( scalar @$aref_genotype == $number_of_colours )
			{
				push @colour_results, $aref_llk_hom1->[$i];
				push @colour_results, $aref_llk_het->[$i];
				push @colour_results, $aref_llk_hom2->[$i];
				push @colour_results,
				  $aref_genotype->[$i];    #br1/br2 `NOT recenessarily ref/alt
			}
			else
			{
				push @colour_results, -1;
				push @colour_results, -1;
				push @colour_results, -1;
				push @colour_results, "0/1";

			}
		}
		push @big_results, \@colour_results;
	}

	return \@big_results;

}

sub get_confidence_haploid
{

#args are 3 log likelihoods of genotypes. we want difference between max and next
	my ( $A, $B ) = @_;

	my @arr;
	push @arr, $A;
	push @arr, $B;

	my @sorted_arr = sort { $a <=> $b } @arr;     #sorted in ascending order
	my $conf = $sorted_arr[1] - $sorted_arr[0];
	my $rounded_conf = sprintf( "%.2f", $conf );

	return $rounded_conf;
}

sub get_confidence_diploid
{

#args are 3 log likelihoods of genotypes. we want difference between max and next
	my ( $A, $B, $C ) = @_;

	my @arr;
	push @arr, $A;
	push @arr, $B;
	push @arr, $C;

	my @sorted_arr = sort { $a <=> $b } @arr;     #sorted in ascending order
	my $conf = $sorted_arr[2] - $sorted_arr[1];
	my $rounded_conf = sprintf( "%.2f", $conf );

	return $rounded_conf;
}

sub print_all_genotypes_and_covgs
{
	my (
		$fh,             $aref_br1_cov,  $aref_br2_cov,  $which_is_ref,
		$classification, $class_llk_rep, $class_llk_var, $genotype,
		$llk_hom1,       $llk_het,       $llk_hom2,      $href_pop_conf,
		$name
	) = @_;

	my $do_we_have_genotypes;
	my $num_samples = $number_of_colours;
	if ( $reference_colour >= 0 )
	{
		$num_samples = $num_samples - 1;
	}
	if ( $pooled_colour >= 0 )
	{
		$num_samples = $num_samples - 1;
	}

	my $do_we_have_pop_filter = 0;
	if ( scalar keys(%$href_pop_conf) > 0 )
	{
		$do_we_have_pop_filter = 1;
	}

	if ( scalar(@$genotype) == 0 )
	{
		$do_we_have_genotypes = 0;
	}
	elsif ( scalar(@$genotype) >= $num_samples )
	{
		$do_we_have_genotypes = 1;
	}
	else
	{
		print(
"Programming error - we have $num_samples samples to genotype but this function: print_all_genotypes_and_covgs has only been given this many: "
		);
		print scalar(@$genotype);
		die("Contact Zam\n");
	}

	my $j;

	for ( $j = 0 ; $j < $number_of_colours ; $j++ )
	{
		if ( ( $j == $reference_colour ) || ( $j == $pooled_colour ) )
		{
			if ( $j == $last_sample_col )
			{
				print $fh "\n";
			}
			next;
		}

#my $aref_lik_and_genotypes_br1_v_br2 = new_genotyper($genotype, $llk_hom1, $llk_het, $llk_hom2);
		my $cov;
		my $gtype;
		my $confidence = -99999;
		my $no_data=0;
		my $gls = "";
		if ( $do_we_have_genotypes == 1 )
		{
		    if ( $ploidy == 2 )
		    {
			$confidence =
			    get_confidence_diploid( $llk_hom1->[$j], $llk_het->[$j],
						    $llk_hom2->[$j] );

			if ($which_is_ref==1)
			{
			    $gls = $llk_hom1->[$j].",".$llk_het->[$j].",".$llk_hom2->[$j];
			}
			else
			{
			    $gls = $llk_hom2->[$j].",".$llk_het->[$j].",".$llk_hom1->[$j];#switched
			}

		    }
		    else
		    {
			$confidence =
			    get_confidence_haploid( $llk_hom1->[$j], $llk_hom2->[$j] );

			if ($which_is_ref==1)
			{
			    $gls = $llk_hom1->[$j].",".$llk_hom2->[$j];
			}
			else
			{
			    $gls = $llk_hom2->[$j].",".$llk_hom1->[$j];#switched
			}

		    }
		    
		    
		    if ( ($aref_br1_cov->[$j]==0) && ($aref_br2_cov->[$j]==0) )
		    {
			$no_data=1;
		    }
		}
		my $site_conf = -99999;
		my $stripped_name = $name;
		my $prefix_under;
		if ($prefix ne "")
		{
		    $prefix_under = $prefix."_";
		    $stripped_name =~ s/$prefix_under//;
		}
		if ( $do_we_have_pop_filter == 1 )
		{
		    if (exists $href_pop_conf->{$name})
		    {
			$site_conf = sprintf( "%.2f", $href_pop_conf->{$name} );
		    }
		    elsif (exists $href_pop_conf->{$stripped_name})
		    {
			$site_conf = sprintf( "%.2f", $href_pop_conf->{$stripped_name});
		    }
		    else
		    {
			$site_conf = "ERROR";
			print "Cannot find site conf for $name or $stripped_name\n";
			print "keys are ";
			print join("\n", keys %$href_pop_conf);
			print "There is some problem with your input file. Contact Zam with your commandline and describe your situation\n";
			die();
		    }
		}
		if ( $which_is_ref == 1 )
		{
			if ( ( $do_we_have_genotypes == 1 ) && ($no_data==0) )
			{
			    $gtype = $genotype->[$j];
			}
			else
			{
			    $gtype="./.";
			}
			$cov = $aref_br1_cov->[$j] . "," . $aref_br2_cov->[$j];
		}
		else
		{
			if ( ($do_we_have_genotypes == 1 ) && ($no_data==0) )
			{
				$gtype = switch_genotype( $genotype->[$j] );
			}
			else
			{
			    $gtype="./.";
			}
			if ( scalar @$aref_br2_cov != $number_of_colours )
			{
				print("Parsing problem. Don't have enough covgs for branch2 - just have ");
				print scalar @$aref_br2_cov;
				die();
			}

			if ( scalar @$aref_br1_cov != $number_of_colours )
			{
				print("Parsing problem. Don't have enough covgs for branch1 - just have ");
				print scalar @$aref_br1_cov;
				die();

			}
			$cov = $aref_br2_cov->[$j] . "," . $aref_br1_cov->[$j];
		}

		if ( ( $do_we_have_pop_filter == 1 ) && ( $do_we_have_genotypes == 1 ) )
		{
			print $fh "$gtype:$cov:$confidence:$site_conf";
		}
		elsif ( $do_we_have_genotypes == 1 )
		{
			print $fh "$gtype:$cov:$confidence";
		}
		elsif ( $do_we_have_pop_filter == 1 )
		{
			print $fh "$gtype:$cov:$site_conf";
		}
		else
		{
			print $fh "$cov";
		}  
		if ($print_gls eq "yes")
		{
		    print $fh ":$gls";
		}

		
		if ( $j == $last_sample_col )
		{
			print $fh "\n";
		}
		else
		{
			print $fh "\t";
		}
	}

}

sub get_simple_vcf_entry_pos_and_alleles

{

	my (
		$name,                            $str,
		$pos,                             $flank5p,
		$br1_seq,                         $br2_seq,
		$flank3p,                         $which_br_ref,
		$align_num_bp_agreement_at_start, $align_num_bp_agreement_at_end,
		$align_dir,                       $num_snps_align,
		$num_indels_align,                $href_var_name_to_cut_flank
	) = @_;

	## if the 5p flank is >1000bp, then stampy fails, so we cut the 5p flank and only take the last 1000bp.
	if ( exists $href_var_name_to_cut_flank->{$name} )
	{
		##the flank5p passedin is directly taken from the callfile
		$flank5p = substr( $flank5p, -1000 );
		if ( length($flank5p) != 1000 )
		{
			die("perl issue with $flank5p");
		}
	}

	## snp for single SNP, else not_snp
	my $what_type = determine_type(
		$num_snps_align,  $num_indels_align,
		length($br1_seq), length($br2_seq)
	);

	my $vcf_entry_ref_allele;
	my $vcf_entry_alt_allele;
	my $vcf_entry_pos;
	## if 5p flank maps in + direction
	if ( $str == 0 )
	{
		##bloody vcf format wants the ref allele to start at the first variant base for a SNP, and the base BEFORE that for a non-snp
		if ( $what_type eq "snp" )
		{
			$vcf_entry_pos = $pos + length($flank5p);
		}
		else
		{
			$vcf_entry_pos = $pos + length($flank5p) - 1;
		}

		if ( $which_br_ref eq "1" )
		{
			## if br1 and br2 are the same at the end, cut off that bit

			if ( $what_type eq "snp" )
			{
				$vcf_entry_ref_allele =
				  substr( $br1_seq, 0,
					length($br1_seq) - $align_num_bp_agreement_at_end );
				$vcf_entry_alt_allele =
				  substr( $br2_seq, 0,
					length($br2_seq) - $align_num_bp_agreement_at_end );

			}
			else
			{
				$vcf_entry_ref_allele = substr( $flank5p, -1 )
				  . substr( $br1_seq, 0,
					length($br1_seq) - $align_num_bp_agreement_at_end );

				$vcf_entry_alt_allele = substr( $flank5p, -1 )
				  . substr( $br2_seq, 0,
					length($br2_seq) - $align_num_bp_agreement_at_end );

			}
		}
		elsif ( $which_br_ref eq "2" )
		{
			## if br1 and br2 are the same at the end, cut off that bit

			if ( $what_type eq "snp" )
			{
				$vcf_entry_ref_allele =
				  substr( $br2_seq, 0,
					length($br2_seq) - $align_num_bp_agreement_at_end );
				$vcf_entry_alt_allele =
				  substr( $br1_seq, 0,
					length($br1_seq) - $align_num_bp_agreement_at_end );

			}
			else
			{
				$vcf_entry_ref_allele = substr( $flank5p, -1 )
				  . substr( $br2_seq, 0,
					length($br2_seq) - $align_num_bp_agreement_at_end );

				$vcf_entry_alt_allele = substr( $flank5p, -1 )
				  . substr( $br1_seq, 0,
					length($br1_seq) - $align_num_bp_agreement_at_end );

			}
		}
		else
		{
			$vcf_entry_ref_allele =
			  substr( $br1_seq, 0,
				length($br1_seq) - $align_num_bp_agreement_at_end );

			$vcf_entry_alt_allele =
			  substr( $br2_seq, 0,
				length($br2_seq) - $align_num_bp_agreement_at_end );


#return (0,0,0,"This var has both alleles entirely in the reference: $which_br_ref");
		}
	}
	##else, 5prime flank mapped in the reverse direction
	elsif ( $str == 16 )
	{
	    ## for non-SNPs, you will need the last base before the variant. In some cases, one needs to get this from the 
	    ## (reverse complement of) the 3p flank. However in some of these cases, the 3p flank is zero!
	    ## when that happens, use a placement character Z. At a later stage, when we have a sorted VCF, we will go 
	    ## through the reference fasta once, and replace the Z characters with the relevant correct characters
	    ## from the reference (one base before the variant in each case).
	    my $last_base_of_rev_comp_of_3p_flank = "Z";
	    if (length($flank3p)>0)
	    {
		#chomp $flank3p;
		$last_base_of_rev_comp_of_3p_flank =  substr( rev_comp($flank3p), -1 );
	    } 

		my $br1_excepting_agreement_at_end =
		  substr( $br1_seq, 0,
			length($br1_seq) - $align_num_bp_agreement_at_end );
		my $br2_excepting_agreement_at_end =
		  substr( $br2_seq, 0,
			length($br2_seq) - $align_num_bp_agreement_at_end );

		my $last_base_in_branch_before_variant;
		if ( $align_num_bp_agreement_at_end > 0 )
		{
			$last_base_in_branch_before_variant =
			  substr( $br1_seq,
				length($br1_seq) - $align_num_bp_agreement_at_end, 1 );
		}

		if ( $which_br_ref eq "1" )
		{
			## if br1 and br2 are the same at the end, cut off that bit

			if ( $what_type eq "snp" )
			{
				$vcf_entry_pos =
				  $pos - length($br1_seq) + $align_num_bp_agreement_at_end;
			}
			else
			{
				$vcf_entry_pos =
				  $pos - 1 - length($br1_seq) + $align_num_bp_agreement_at_end;
			}

			if ( $what_type eq "snp" )
			{
				$vcf_entry_ref_allele =
				  rev_comp($br1_excepting_agreement_at_end);
				$vcf_entry_alt_allele =
				  rev_comp($br2_excepting_agreement_at_end);

#$vcf_entry_ref_allele =  rev_comp(substr($br1_seq, - (length($br1_seq) - $align_num_bp_agreement_at_start) ) );
#$vcf_entry_alt_allele =  rev_comp(substr($br2_seq, - (length($br2_seq) - $align_num_bp_agreement_at_start) ) );
			}
			else
			{
				if ( $align_num_bp_agreement_at_end == 0 )
				{
					$vcf_entry_ref_allele =
					  $last_base_of_rev_comp_of_3p_flank
					  . rev_comp($br1_excepting_agreement_at_end);
					$vcf_entry_alt_allele =
					  $last_base_of_rev_comp_of_3p_flank
					  . rev_comp($br2_excepting_agreement_at_end);

				}
				else
				{
					$vcf_entry_ref_allele =
					    rev_comp($last_base_in_branch_before_variant)
					  . rev_comp($br1_excepting_agreement_at_end);
					$vcf_entry_alt_allele =
					    rev_comp($last_base_in_branch_before_variant)
					  . rev_comp($br2_excepting_agreement_at_end);

				}

#$vcf_entry_ref_allele =  substr( rev_comp($flank3p), -1)  .rev_comp(substr($br1_seq, - (length($br1_seq) - $align_num_bp_agreement_at_start) ) );
#$vcf_entry_alt_allele =  substr( rev_comp($flank3p), -1)  .rev_comp(substr($br2_seq, - (length($br2_seq) - $align_num_bp_agreement_at_start) ) );
			}
		}
		elsif ( $which_br_ref eq "2" )
		{
			## if br1 and br2 are the same at the end, cut off that bit

			if ( $what_type eq "snp" )
			{
				$vcf_entry_pos =
				  $pos - length($br2_seq) + $align_num_bp_agreement_at_end;
				## use if mapping 5p and br1     $vcf_entry_pos        =  $pos+length($br1_seq) -length($br2_seq) + $align_num_bp_agreement_at_end;
			}
			else
			{
				$vcf_entry_pos =
				  $pos - length($br2_seq) - 1 + $align_num_bp_agreement_at_end;
				## use if mapping 5p and br1  $vcf_entry_pos        =  $pos+length($br1_seq) -length($br2_seq)-1 + $align_num_bp_agreement_at_end;
			}

			if ( $what_type eq "snp" )
			{
				$vcf_entry_ref_allele =
				  rev_comp($br2_excepting_agreement_at_end);
				$vcf_entry_alt_allele =
				  rev_comp($br1_excepting_agreement_at_end);

			}
			else
			{
				if ( $align_num_bp_agreement_at_end == 0 )
				{
				    $vcf_entry_ref_allele =
					$last_base_of_rev_comp_of_3p_flank
					. rev_comp($br2_excepting_agreement_at_end);
				    $vcf_entry_alt_allele =
					$last_base_of_rev_comp_of_3p_flank
					. rev_comp($br1_excepting_agreement_at_end);

				}
				else
				{
					$vcf_entry_ref_allele =
					    rev_comp($last_base_in_branch_before_variant)
					  . rev_comp($br2_excepting_agreement_at_end);
					$vcf_entry_alt_allele =
					    rev_comp($last_base_in_branch_before_variant)
					  . rev_comp($br1_excepting_agreement_at_end);

				}

			}
		}
		else
		{
			return ( 0, 0, 0,
"This var has both alleles entirely in the reference: $which_br_ref"
			);
		}

	}
	else
	{
		return ( 0, 0, 0, "Did not map" );
	}


	return ( $vcf_entry_pos, $vcf_entry_ref_allele, $vcf_entry_alt_allele,
		"0" );

}

sub rev_comp
{
	my ($seq) = @_;
	
	if ($seq=~/(\S+)\s+/)
	{
	    $seq = $1;
	}
	my $r_seq = reverse($seq);
	$r_seq =~ tr/acgtACGT/tgcaTGCA/;

	# print join(" ",$seq,$r_seq),"\n";
	return $r_seq;

}

sub determine_type
{
	my ( $num_snps_align, $num_indels_align, $lenbr1, $lenbr2 ) = @_;

	my $what_type;
	if (   ( $num_snps_align == 1 )
		&& ( $num_indels_align == 0 )
		&& ( $lenbr1 == $lenbr2 ) )
	{
		$what_type = "snp";
	}
	else
	{
		$what_type = "not_snp";
	}

	return $what_type;
}

sub get_next_var_from_callfile
{
	my ( $fh, $p_loidy ) = @_;

	my $line = "";
	my $varname;
	my $flank5p;
	my $br1;
	my $br2;
	my $flank3p;

	my $classification = "VARIANT";
	my $class_llk_rep  = -10;
	my $class_llk_var  = -1;

	## variables I will return
	my @arr_br1_covgs = ();    ## start covg + jumps - array, one per colour
	my @arr_br2_covgs = ();    ## start covg + jumps - array, one per colour
	my @arr_br1_min_covg=();
	my @arr_br2_min_covg=();
	my @arr_geno      = ();
	my @arr_llk_hom1  = ();
	my @arr_llk_hom2  = ();
	my @arr_llk_het   = ();
	my $extra_info="";
	while (( $line !~ /(\w*var_\d+)_5p_flank/ )
		&& ( $line !~ /PASSES/ )
		&& ( $line !~ /FAILS/ )
		&& ( $line !~ /Colour\/sample/ ) )
	{
		if ( eof($fh) )
		{
			return ( "EOF", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 );
		}
		else
		{
			$line = <$fh>;
		}
	}

	if ( $line =~ /PASSES/ )
	{
		$classification = "VARIANT";
		if ( eof($fh) )
		{
			return ( "EOF", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0 );
		}

		$line = <$fh>;
		chomp $line;

		if ( $line =~ /llk_var:(\S+).+llk_rep:(\S+)/ )
		{
			$class_llk_rep = $1;
			$class_llk_var = $2;
		}
		else
		{
			die(
"Parsing issue. Found model selecion result but cant find likelihoods on $line\n"
			);
		}
		if ( eof($fh) )
		{
			return ( "EOF", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 );
		}
		$line = <$fh>;

	}
	elsif ( $line =~ /FAILS/ )
	{
		$classification = "REPEAT";
		if ( eof($fh) )
		{
			return ( "EOF", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 );
		}

		$line = <$fh>;
		chomp $line;

		if ( $line =~ /llk_var:(\S+).+llk_rep:(\S+)/ )
		{
			$class_llk_rep = $1;
			$class_llk_var = $2;
		}
		else
		{
			die(
"Parsing issue. Found model selecion result but cant find likelihoods on $line\n"
			);
		}
		if ( eof($fh) )
		{
			return ( "EOF", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 );
		}
		$line = <$fh>;
	}

	if ( $line =~ /Colour\/sample/ )
	{
		my $j;
		for ( $j = 0 ; $j < $number_of_colours ; $j++ )
		{
			$line = <$fh>;
			chomp $line;

			my @sp = split( /\t/, $line );
			if ( $sp[1] eq "HOM1" )
			{
				push @arr_geno, "0/0";
			}
			elsif ( $sp[1] eq "HET" )
			{
				push @arr_geno, "0/1";
			}
			elsif ( $sp[1] eq "HOM2" )
			{
				push @arr_geno, "1/1";
			}
			elsif ( $sp[1] eq "NO_CALL" )
			{
				push @arr_geno, "-1/-1"
				  ; ## should never see the light of day, and easy to spot if appears in VCF
			}
			else
			{
				die("Bad genotype on $line");
			}

			if ( $p_loidy == 2 )
			{
				push @arr_llk_hom1, $sp[2];
				push @arr_llk_het,  $sp[3];
				push @arr_llk_hom2, $sp[4];
			}
			else
			{
				push @arr_llk_hom1, $sp[2];
				push @arr_llk_het,  -9999999;
				push @arr_llk_hom2, $sp[3];
			}
		}
		$line = <$fh>;

	}

	if ( $line =~ /(\w*var_\d+)_5p_flank/ )
	{
		$varname = $1;
		if ($prefix ne "")
		{
		    $varname = $prefix."_".$varname;
		}

		if ( $varname eq "" )
		{
			die(
"Found this line :$line in the callfile - expected it to be of the form var_<NUMBER>_5p_flank\n"
			);
		}

		if (exists $call_to_extra_info{$varname})
		{
		    $extra_info=$call_to_extra_info{$varname};
		}
		elsif ($varname =~ /(\S+)_sub/)
		{
		    my $rawname = $1;
		    if (exists $call_to_extra_info{$rawname})
		    {
			$extra_info=$call_to_extra_info{$rawname};
		    }
		}
		else
		{
		    #print "No extra info (eg Kmer value)  for $varname\n";
		}


		$flank5p = <$fh>;
		chomp $flank5p;
		<$fh>;    #ignore br1 read id
		$br1 = <$fh>;
		chomp $br1;
		<$fh>;    #ignore br2 read id
		$br2 = <$fh>;
		chomp $br2;
		<$fh>;    #ignore 3p flank read id
		$flank3p = <$fh>;
		chomp $flank3p;
		$line    = <$fh>;

		if ( $line =~ /extra information/ )
		{
			<$fh>;
		}
		$line = <$fh>;
		chomp $line;
		if ( $line !~ /branch1 coverages/ )
		{
			die("Expected to see \"branch 1 coverages\" but instead saw $line");
		}
		my $z;
		for ( $z = 0 ; $z < $number_of_colours ; $z++ )
		{
			$line = <$fh>;
			chomp $line;
			if ( $line !~ /Covg in Colour/ )
			{
				die("Expected to see \"Covg in Colour\" but instead saw $line");
			}
			$line = <$fh>;
			chomp $line;
			my @br1 = split( /\s+/, $line );

			# default/standard Cortex
			push @arr_br1_covgs, get_num_reads( \@br1 );
			
			## contemplating moving to this, using median 
			#push @arr_br1_covgs, get_num_reads_using_median( $line, $colour_to_readlen{$z}, $kmer); 

			push @arr_br1_min_covg, get_min_covg(\@br1);
		}
		$line = <$fh>;
		chomp $line;
		if ( $line !~ /branch2 coverages/ )
		{
			die("Expected to see \"branch 2 coverages\" but instead saw $line");
		}

		for ( $z = 0 ; $z < $number_of_colours ; $z++ )
		{
			$line = <$fh>;
			chomp $line;

			if ( $line !~ /Covg in Colour/ )
			{
				die("Expected to see \"Covg in Colour\" but instead saw $line");
			}
			$line = <$fh>;
			chomp $line;
			my @br2 = split( /\s+/, $line );
			
			#default/standard
			push @arr_br2_covgs, get_num_reads( \@br2 );

			#contemplating moving to median:
			#push @arr_br2_covgs, get_num_reads_using_median($line, $colour_to_readlen{$z}, $kmer);
			
			push @arr_br2_min_covg, get_min_covg(\@br2);
		}

		##determine which is ref allele
		my $which_is_ref = "b";

		if ( $reference_colour == -1 )
		{
			$which_is_ref = 1;
		}
		##elsif (( $arr_br1_covgs[$reference_colour] >= 1 )
		##	&& ( $arr_br2_covgs[$reference_colour] == 0 ) )
		elsif (( $arr_br1_min_covg[$reference_colour] >= 1 )
			&& ( $arr_br2_min_covg[$reference_colour] == 0 ) )
		{
			$which_is_ref = 1;
		}
		#elsif (( $arr_br1_covgs[$reference_colour] == 0 )
		#	&& ( $arr_br2_covgs[$reference_colour] >= 1 ) )
		elsif (( $arr_br1_min_covg[$reference_colour] == 0 )
			&& ( $arr_br2_min_covg[$reference_colour] >= 1 ) )
		{
			$which_is_ref = 2;
		}
		#elsif (( $arr_br1_covgs[$reference_colour] >= 1 )
		#	&& ( $arr_br2_covgs[$reference_colour] >= 1 ) )
		elsif (( $arr_br1_min_covg[$reference_colour] >= 1 )
			&& ( $arr_br2_min_covg[$reference_colour] >= 1 ) )
		{
			#$which_is_ref = "b";
			$which_is_ref = 1;## ZAm - added during debugging - will get Isaac's scripts to fix up afterwards
		}
		else
		{
			$which_is_ref = "neither";
		}

		my $eof = "";

		return (
			$eof,            $varname,        $flank5p,
			$br1,            $br2,            $flank3p,
			\@arr_br1_covgs, \@arr_br2_covgs, $which_is_ref,
			$classification, $class_llk_rep,  $class_llk_var,
			\@arr_geno,      \@arr_llk_hom1,  \@arr_llk_het,
			\@arr_llk_hom2, $extra_info
		);

	}
	else
	{
		print
	    "Parsing next line of callfile - unexpected error on  $line. Expected the read-id of 5prime flank in callfile to be of the form (optional_text)var_(number)_5p_flank. This bug is probably due to Zam's changes to support the new format of callfiles but also stay legacy format compatible - contact zam\@well.ox.ac.uk\n";

		die();
	}
}

sub get_num_reads_using_median
{
    my ($line, $read_length, $khmer) = @_;
    my @sp = split(/\s+/, $line);
    my $i;
    my $total = $sp[1];
    if (scalar @sp < 3)
    {
	return 0;
    }

    pop(@sp);
    shift(@sp); ##remove firsdt and last elements
    my $stat = Statistics::Descriptive::Full->new();
    $stat ->add_data(@sp);

    my $median = $stat->median();

    if ($read_length-$khmer+1>scalar(@sp) )## the variant is shorted than the effective read length
    {
	return int($median);
    }
    else
    {
	## total number of kmer-covg = median * length
	my $num_kmercovg = scalar(@sp) * $median;
	## total reads = = number bases/read length
	my $num_reads = int($num_kmercovg/($read_length-$khmer+1));
	#print "Num reads is $num_reads. Median is $median zaz. num kmercovg is $num_kmercovg, r-k is ";
	#print $read_length-$kmer+1;
	#print "\n";
	return $num_reads;
    }
}
sub get_num_reads
{
	my ($aref) = @_;
	my $count;
	if ( scalar(@$aref) > 2 )
	{
		$count = $aref->[1];
	}
	else
	{
		return 0;
	}
	my $i;
	for ( $i = 2 ; $i < scalar(@$aref) - 1 ; $i++ )
	{

		my $jump      = $aref->[$i] - $aref->[ $i - 1 ];
		my $next_jump = -1;
		if ( $i + 1 < scalar(@$aref) - 1 )
		{
			$next_jump = $aref->[ $i + 1 ] - $aref->[$i];
		}

		if ( ( $jump > 0 ) && ( $next_jump != $jump ) )
		{
			$count += $jump;
		}
	}
	return $count;
}


sub get_next_var_from_flank_mapfile
{
	my ($fh) = @_;

	my $line = "@";
	while ( $line =~ /^\@/ )
	{

		if ( eof($fh) )
		{
			return ( "EOF", 0, 0, 0, 0 );
		}
		$line = <$fh>;
	}
	chomp $line;
	my @sp = split( /\t/, $line );
	if ( scalar(@sp) < 10 )
	{
		die("Unexpected num <10 of fields on $line");
	}
	my $name;
	if ( $sp[0] =~ /(\w*var_\d+)/ )
	{
		$name = $1;
		if ($prefix ne "")
		{
		    $name = $prefix."_".$name;
		}
	}
	else
	{
		die("Unexpected format of first field in $line");
	}
	my $strand = $sp[1];
	my $chr    = $sp[2];
	my $coord  = $sp[3];
	my $eof    = "";
	return ( $eof, $name, $strand, $chr, $coord );
}

sub switch_genotype
{
	my ($str) = @_;
	if ( $str =~ /^([01])\/([01])/ )
	{
		my $a = $1;
		my $b = $2;

		if ( $a == 0 )
		{
			$a = 1;
		}
		else
		{
			$a = 0;
		}
		if ( $b == 0 )
		{
			$b = 1;
		}
		else
		{
			$b = 0;
		}

		## might as well return 0/1 rather than 1/0
		my $ret;
		if ( $a . '/' . $b eq "1/0" )
		{
			$ret = "0/1";
		}
		else
		{
			$ret = $a . '/' . $b;
		}

		return $ret;

	}
	elsif ($str eq "-1/-1")
	{
	    return "-1/-1";
	}
	else
	{
		die("Passed bad genotype $str to switch_genotype");
	}
}

## return info about how the alignment of the two branches goes.
## if it turns out that an alignment of the rev comp of the branches looks good,
## then return that info as a binary datum, so can flag the variant as possible inversion
sub get_next_alignment_from_procfile
{

#we almost do not need $dir_of_alignment_of_5pflank, except at the end when we get coordinates of SNPs in a set pf phased SNPS
#similarly for $which_is_ref. Must handle gracefully if this has value b (both) or neither. In both these cases we will ignore the variant
#and we only call this function to make sure we be at the right point in the file to go to the NEXT one, next time we call this. ie just read past this var.
	my ( $fh, $dir_of_alignment_of_5pflank, $which_is_ref ) = @_;

	if ( !defined $which_is_ref )
	{
		die(
			"Passing undefined which_is_ref to get_next_alignment_from_procfile"
		);
	}

	my $line = "";
	while ( $line !~ /START NEW VAR/ )
	{
		if ( eof($fh) )
		{
			return ( "EOF", 0, 0, 0, 0, 0, 0, 0, 0, 0 );
		}
		$line = <$fh>;
	}

	$line = <$fh>;
	my $name;
	if ( $line =~ /^(\S+)/ )
	{
		$name = $1;
	}
	else
	{
		die("Expected name on $line");
	}
	<$fh>;

	my $fw_br1         = "";
	my $rev_br1        = "";
	my $fw_middle      = "";
	my $rev_middle     = "";
	my $fw_br2         = "";
	my $rev_br2        = "";
	my $fw_num_snps    = 0;
	my $rev_num_snps   = 0;
	my $fw_num_indels  = 0;
	my $rev_num_indels = 0;

	$line = <$fh>;
	if ( $line !~ /FORWARD ALIGNMENT/ )
	{
		die(
		    "format issue - this line $line should have said FORWARD ALIGNMENT. Current var is $name\n"
		);
	}
	$fw_br1 = <$fh>;
	chomp $fw_br1;
	$fw_middle = <$fh>;
	chomp $fw_middle;
	$fw_br2 = <$fh>;
	chomp $fw_br2;
	$fw_br1 =~ s/Br1://;
	$fw_br2 =~ s/Br2://;

	if ( $fw_middle =~ /^\s{4}(.+)/ )
	{
		$fw_middle = $1;
	}
	else
	{
		die(
"bad format of fw middle $fw_middle. which does not have 4 space sat the start - put there so it all aligns visually when printed, but nothing to do with the alignment"
		);
	}

	$line = <$fh>;
	chomp $line;

	if ( $line =~ /(\d+)\s+(\d+)/ )
	{
		$fw_num_snps   = $1;
		$fw_num_indels = $2;
	}
	else
	{
		die("Bad format of $line - expect num snps and indels");
	}
	<$fh>;

	my $there_is_rev_alignment = 0;
	$line = <$fh>;
	if ( $line =~ /NO REVERSE ALIGNMENT/ )
	{

	}
	elsif ( $line =~ /REVERSE ALIGNMENT/ )
	{
		$there_is_rev_alignment = 1;
		$rev_br1                = <$fh>;
		chomp $rev_br1;
		$rev_middle = <$fh>;
		chomp $rev_middle;
		$rev_br2 = <$fh>;
		chomp $rev_br2;
		$rev_br1 =~ s/Br1://;
		$rev_br2 =~ s/Br2://;

		if ( $rev_middle =~ /^\s{4}(.+)/ )
		{
			$rev_middle = $1;
		}
		else
		{
			die(
"bad format of rev middle $rev_middle. which does not have 4 space sat the start - put there so it all aligns visually when printed, but nothing to do with the alignment"
			);
		}

		my $line = <$fh>;
		chomp $line;

		if ( $line =~ /(\d+)\s+(\d+)/ )
		{
			$rev_num_snps   = $1;
			$rev_num_indels = $2;
		}
		else
		{
			die("Bad format of $line - expect num snps and indels");
		}
	}

	if ( ( $which_is_ref eq "b" ) || ( $which_is_ref eq "neither" ) )
	{
		## we will probably ignore this variant. But in case the alignment is useful, set ref allele to 1
		$which_is_ref = 1;
	}

	##  if in doubt, say there are 0 bases of agreement at start/end, so that no clever trimming is done

	## 0 =no, 1=inversion (pure), 2=deletion+inversion or insertion + inversion
	my $possible_inversion = 0;
	if ( $there_is_rev_alignment == 1 )
	{
		##then we have an inversion + maybe deletion
		if ( $rev_middle =~ /^\s+[\|]+$/ )
		{
			$possible_inversion = 2;
		}
		elsif ( $rev_middle =~ /^[\|]+\s+$/ )
		{
			$possible_inversion = 2;
		}
		elsif ( $rev_middle =~ /^[\|]+$/ )    ##perfect inversion
		{
			$possible_inversion = 1;
		}
		else
		{
			my $max_consec_matches_fw = get_max_consecutive_matches($fw_middle);
			my $max_consec_matches_rev =
			  get_max_consecutive_matches($rev_middle);

			if ( $max_consec_matches_fw < $max_consec_matches_rev )
			{
				$possible_inversion = 2;
			}

		}
	}

	my $clean_indel = 0;
	if ( $fw_middle =~ /^\s+[\|]+$/ )
	{
		$clean_indel = 1;
	}
	elsif ( $fw_middle =~ /^[\|]+\s+$/ )
	{
		$clean_indel = 1;
	}

	my $num_bases_agreement_at_start = 0;
	my $num_bases_agreement_at_end   = 0;
	my $align_direction;
	my $num_snps;
	my $num_indels;

	$align_direction = "+";
	if (   ( length($fw_br1) == length($fw_br2) )
		&& ( $fw_middle =~ /^([\|]+)/ ) )
	{
		if ( $fw_middle =~ /^([\|]+)/ )
		{
			$num_bases_agreement_at_start = length($1);
		}
	}
	if (   ( length($fw_br1) == length($fw_br2) )
		&& ( $fw_middle =~ /([\|]+)$/ ) )
	{
		if ( $fw_middle =~ /([\|]+)$/ )
		{
			$num_bases_agreement_at_end = length($1);
		}
	}
	$num_snps   = $fw_num_snps;
	$num_indels = $fw_num_indels;

	#if this is a bunch of phased SNPs, then get their positions
	my @coords_of_snps  = ();
	my @indices_of_snps = ();
	my @alleles_of_snps = ();

	my @indices_of_indels           = ();
	my @indices_of_indel_ends       = ();
	my @coords_of_indels            = ();
	my @alleles_of_indels           = ();
	my @indel_flag_needs_extra_base = ()
	  ; #will be 1 if that indel needs an extra base on the front. Should only be possible for ONE of th indels

	my $warning_need_to_add_flank_base = 0;

	if (   ( $num_snps > 0 )
		&& ( $num_indels == 0 )
	  )    ## this can only happen if the two branches are the same length
	{

		if ( $dir_of_alignment_of_5pflank == 0 )
		{
			### first get the SNP coords
			my $j;
			my @sp = split( //, $fw_middle );
			for ( $j = 0 ; $j < length($fw_middle) ; $j++ )
			{
				if ( $sp[$j] eq "\*" )
				{
					push @coords_of_snps, $j;
				}
			}
			##then get the SNP alleles, always in form REF_allele,ALT_allele. Last argument is ignored as in fw direction
			get_snp_alleles( \@coords_of_snps, \@alleles_of_snps, $fw_br1,
				$fw_middle, $fw_br2, "fw", $which_is_ref, -1 );

		}
		elsif ( $dir_of_alignment_of_5pflank == 16 )
		{
			## start at the other end, and ignore the bases of agreement
			## e.g.
			##Br1:CGCCGTTGTTGAGTGTTCTATGGAATTGTCGTTTATTGAGCACAACTACAGCATTT
			##    *|||*||||||||||||||*||||||||||||||*|||||||||||||||||||||
			##Br2:TGCCCTTGTTGAGTGTTCTTTGGAATTGTCGTTTTTTGAGCACAACTACAGCATTT
			##
			##                                       <---ignore these--->
			##                                      ^start here and work left

			my $j;
			my @sp = split( //, $fw_middle );
			my $leng = scalar(@sp);
			##first find first base of disagreement starting from the end
			my $fwd_pos_of_last_snp = $leng - 1;

			for ( $j = $leng - 1 ; $j >= 0 ; $j-- )
			{
				if ( $sp[$j] ne "\|" )
				{
					$fwd_pos_of_last_snp = $j;
					last;
				}
			}
			for ( $j = 0 ; $j <= $fwd_pos_of_last_snp ; $j++ )
			{
				if ( $sp[ $fwd_pos_of_last_snp - $j ] eq "\*" )
				{
					push @coords_of_snps, $j;
				}
			}

			##then get the SNP alleles, always in form REF_allele,ALT_allele.  Note we DELIBERATELY pass in fw_br, fw_middle etc AND "rev".
			get_snp_alleles( \@coords_of_snps, \@alleles_of_snps, $fw_br1,
				$fw_middle, $fw_br2, "rev", $which_is_ref,
				$fwd_pos_of_last_snp );

		}
	}

	elsif (( $possible_inversion == 0 )
		&& ( $num_snps > 0 )
		&& ( $num_indels > 0 ) )
	{

		if ( $dir_of_alignment_of_5pflank == 0 )
		{
			### first get the SNP coords and indices
			my $j;
			my @sp        = split( //, $fw_middle );
			my @split_br1 = split( //, $fw_br1 );
			my @split_br2 = split( //, $fw_br2 );

			my $underscore_on_ref_allele_so_far = 0;

			my $min;    ##zam added all this min stuff in debug
			if ( $which_is_ref == 1 )
			{
				if ( length($fw_middle) > scalar(@split_br1) )
				{
					$min = scalar(@split_br1);
				}
				else
				{
					$min = length($fw_middle);
				}
			}
			else
			{
				if ( length($fw_middle) > scalar(@split_br2) )
				{
					$min = scalar(@split_br2);
				}
				else
				{
					$min = length($fw_middle);
				}
			}

			for ( $j = 0 ; $j < $min ; $j++ )

			  # removed in debug for ($j=0; $j<length($fw_middle); $j++)
			{
				if ( ( $which_is_ref == 1 ) && ( $split_br1[$j] eq "_" ) )
				{
					$underscore_on_ref_allele_so_far++;
				}
				elsif ( ( $which_is_ref == 2 ) && ( $split_br2[$j] eq "_" ) )
				{
					$underscore_on_ref_allele_so_far++;
				}

				if ( $sp[$j] eq "\*" )
				{

					#check base before/after is not an indel
					if (
						(
							   ( $j > 0 )
							&& ( $j < length($fw_middle) - 1 )
							&& ( $sp[ $j - 1 ] ne " " )
							&& ( $sp[ $j + 1 ] ne " " )
						)
						|| ( ( $j == 0 ) && ( $sp[ $j + 1 ] ne " " ) )
						|| (   ( $j == length($fw_middle) - 1 )
							&& ( $sp[ $j - 1 ] ne " " ) )
					  )
					{
						push @indices_of_snps, $j;
						push @coords_of_snps,
						  $j - $underscore_on_ref_allele_so_far;
					}
				}

			}
			##then get the SNP alleles, always in form REF_allele,ALT_allele. Last argument is ignored as in fw direction
			get_snp_alleles( \@indices_of_snps, \@alleles_of_snps, $fw_br1,
				$fw_middle, $fw_br2, "fw", $which_is_ref, -1 );

			## now get the indel INDICES of where the indels are in these strings, and hence will get the COORDINATES
			##  eg ig we have AAA_T__C on the ref allele, then the index of C is 7, and the coordinate is 5
			## will get the coord just before the indel in the forward direction along the ref, and for ref/alt alleles give that prior base also.
			## look for gaps in ||| ||||| that are preceded by two || and succeeded also by two.

			$underscore_on_ref_allele_so_far = 0;

			for ( $j = 0 ; $j < length($fw_middle) ; $j++ )
			{

				if ( $sp[$j] eq " " )
				{
					if (
						   ( ( $j == 0 ) )
						|| ( ( $j == 1 ) && ( $sp[ $j - 1 ] ne " " ) )
						|| (   ( $j > 1 )
							&& ( $sp[ $j - 1 ] ne " " )
							&& ( $sp[ $j - 2 ] ne " " ) )
					  )
					{

						my $indel_start_coord =
						  $j -
						  $underscore_on_ref_allele_so_far -
						  1;    ##base before indel
						my $indel_start_index =
						  $j - 1;    ## base before first base of space

						##we have found an indel. Now find the end: - keep going to you see two |, or reach end
						while (
							(
								!(
									   ( $j < scalar(@sp) - 2 )
									&& ( $sp[ $j + 1 ] eq "\|" )
									&& ( $sp[ $j + 2 ] eq "\|" )
								)
								&& !(
									   ( $j == scalar(@sp) - 2 )
									&& ( $sp[ $j + 1 ] eq "\|" )
								)
								&& !( ( $j == scalar(@sp) - 1 ) )
							)
							&& ( $j < scalar(@sp) )
						  )
						{

							if (   ( $which_is_ref == 1 )
								&& ( $split_br1[$j] eq "_" ) )
							{
								$underscore_on_ref_allele_so_far++;
							}
							elsif (( $which_is_ref == 2 )
								&& ( $split_br2[$j] eq "_" ) )
							{
								$underscore_on_ref_allele_so_far++;
							}

							$j++;
						}

						push @indices_of_indels,     $indel_start_index;
						push @coords_of_indels,      $indel_start_coord;
						push @indices_of_indel_ends, $j;

					}
				}

				if ( ( $which_is_ref == 1 ) && ( $split_br1[$j] eq "_" ) )
				{
					$underscore_on_ref_allele_so_far++;
				}
				elsif ( ( $which_is_ref == 2 ) && ( $split_br2[$j] eq "_" ) )
				{
					$underscore_on_ref_allele_so_far++;
				}

			}

			## now get the indel alleles
			get_indel_alleles(
				\@indices_of_indels, \@indices_of_indel_ends,
				\@alleles_of_indels, \@indel_flag_needs_extra_base,
				$fw_br1,             $fw_middle,
				$fw_br2,             "fw",
				$which_is_ref,       -1
			);

		}

		## get forward INDEX of base just after last space in the indel -  ie when we reverse, will be the base just befroe the indel
		## also, get the coord of that same base.
		elsif ( $dir_of_alignment_of_5pflank == 16 )
		{
			## e.g.
			##Br1:CGCCGTTGTTGAGTGTTCTATGG__TTGTCGTTTATTGAGCACAACTACAGCATTT
			##    *|||*||||||||||||||*|||  |||||||||||||||||||||||||||||||
			##Br2:TGCCCTTGTTGAGTGTTCTTTGGAATTGTCGTTTATTGAGCACAACTACAGCATTT
			##                             ^find this base

			my $j;
			my @sp        = split( //, $fw_middle );
			my @split_br1 = split( //, $fw_br1 );
			my @split_br2 = split( //, $fw_br2 );
			my $underscore_on_ref_allele_so_far = 0;

			my $leng = scalar(@sp);
			##first find first base of disagreement starting from the end
			my $fwd_pos_of_last_non_match = $leng - 1;

			for ( $j = $leng - 1 ; $j >= 0 ; $j-- )
			{
				if ( $sp[$j] ne "\|" )
				{
					$fwd_pos_of_last_non_match = $j;
					last;
				}
			}
			for ( $j = 0 ; $j <= $fwd_pos_of_last_non_match ; $j++ )
			{

				if ( $sp[ $fwd_pos_of_last_non_match - $j ] eq "\*" )
				{
					push @indices_of_snps, $j;
					push @coords_of_snps, $j - $underscore_on_ref_allele_so_far;
				}

				if (   ( $which_is_ref == 1 )
					&& ( $split_br1[ $fwd_pos_of_last_non_match - $j ] eq "_" )
				  )
				{
					$underscore_on_ref_allele_so_far++;
				}
				elsif (( $which_is_ref == 2 )
					&& ( $split_br2[ $fwd_pos_of_last_non_match - $j ] eq "_" )
				  )
				{
					$underscore_on_ref_allele_so_far++;
				}

			}

			##then get the SNP alleles, always in form REF_allele,ALT_allele.  Note we DELIBERATELY pass in fw_br, fw_middle etc AND "rev".
			get_snp_alleles( \@indices_of_snps, \@alleles_of_snps, $fw_br1,
				$fw_middle, $fw_br2, "rev", $which_is_ref,
				$fwd_pos_of_last_non_match );

			##now do indels.

			$underscore_on_ref_allele_so_far = 0;

			for ( $j = 0 ; $j <= $fwd_pos_of_last_non_match ; $j++ )
			{

				if ( $sp[ $fwd_pos_of_last_non_match - $j ] eq " " )
				{
					if (    ##there are two | on the RHS
						( ( $j == 0 ) )
						|| (
							( $j == 1 )
							&& ( $sp[ $fwd_pos_of_last_non_match - $j + 1 ] eq
								"\|" )
						)
						|| (
							( $j > 1 )
							&& ( $sp[ $fwd_pos_of_last_non_match - $j + 1 ] eq
								"\|" )
							&& ( $sp[ $fwd_pos_of_last_non_match - $j + 2 ] eq
								"\|" )
						)
					  )

					{

#                 v one after fwd pos of last non-match - in fact might actually be AFTER the whole string.
#  A___CCCCCCCCCC_C
#      ^- 10 is the index_start_coord for this deletion ___
						##indel_start_coord is the absolute size of the difference between the fwd_pos_of_last_non_match, and the base after (going forwards) the start
						my $indel_start_coord =
						  $j -
						  $underscore_on_ref_allele_so_far -
						  1;    ##base before indel
						my $indel_start_index =
						  $fwd_pos_of_last_non_match -
						  $j + 1;    ## base before first base of space

						##we have found an indel. Now find the end: - keep going to you see two |, or reach end
						while (
							!(
								(
									( $j <= $fwd_pos_of_last_non_match - 2 )
									&& ( $sp[ $fwd_pos_of_last_non_match - $j -
										1 ] eq "\|" )
									&& ( $sp[ $fwd_pos_of_last_non_match - $j -
										2 ] eq "\|" )
								)
								|| (
									( $j == $fwd_pos_of_last_non_match - 1 )
									&& ( $sp[ $fwd_pos_of_last_non_match - $j -
										1 ] eq "\|" )
								)
								|| ( ( $j == $fwd_pos_of_last_non_match ) )
							)
							&& ( $j <= $fwd_pos_of_last_non_match )
						  )
						{

							if (
								( $which_is_ref == 1 )
								&& ( $split_br1[ $fwd_pos_of_last_non_match -
									$j ] eq "_" )
							  )
							{
								$underscore_on_ref_allele_so_far++;
							}
							elsif (
								( $which_is_ref == 2 )
								&& ( $split_br2[ $fwd_pos_of_last_non_match -
									$j ] eq "_" )
							  )
							{
								$underscore_on_ref_allele_so_far++;
							}

							$j++;
						}

						push @indices_of_indels, $indel_start_index;
						push @coords_of_indels,  $indel_start_coord;
						push @indices_of_indel_ends,
						  $fwd_pos_of_last_non_match - $j;

					}
				}

				if (   ( $which_is_ref == 1 )
					&& ( $split_br1[ $fwd_pos_of_last_non_match - $j ] eq "_" )
				  )
				{
					$underscore_on_ref_allele_so_far++;
				}
				elsif (( $which_is_ref == 2 )
					&& ( $split_br2[ $fwd_pos_of_last_non_match - $j ] eq "_" )
				  )
				{
					$underscore_on_ref_allele_so_far++;
				}

			}

			## now get the indel alleles
			get_indel_alleles(
				\@indices_of_indels, \@indices_of_indel_ends,
				\@alleles_of_indels, \@indel_flag_needs_extra_base,
				$fw_br1,             $fw_middle,
				$fw_br2,             "rev",
				$which_is_ref,       $fwd_pos_of_last_non_match
			);

		}

	}

	my $eof = "";

	return (
		$eof,                          $name,
		$num_bases_agreement_at_start, $num_bases_agreement_at_end,
		$align_direction,              $num_snps,
		$num_indels,                   $fw_br1,
		$fw_middle,                    $fw_br2,
		$possible_inversion,           $clean_indel,
		\@coords_of_snps,              \@alleles_of_snps,
		\@coords_of_indels,            \@alleles_of_indels,
		\@indel_flag_needs_extra_base
	);

}

## string of form |||||| |||| | | || |          || || | ||
sub get_max_consecutive_matches
{
	my ($str) = @_;

	my $max = 0;

	my @sp = split( /\s+/, $str );
	foreach my $word (@sp)
	{
		if ( length($word) > $max )
		{
			$max = length($word);
		}
	}
	return $max;
}

##assumes which_is_ref is either 1 or 2
sub get_snp_alleles
{
	###index_of_last_non_match is only used in the case dir=rev. branches 1,2 agree for some number of bases at the end
	### - this gives the index (in normal forward coords) of the first base of difference from the end.
	## eg if we had
	#  Br1:CGCCGTTGT
	#      *|||*||||
	#  Br2:TGCCCTTGT
	# then this should be 4. ie coord of last star!

	my ( $aref_snp_coor, $aref_snp_allel, $b1, $middle, $b2, $dir,
		$which_is_ref, $index_of_last_non_match )
	  = @_;

	##double check
	if (   ( length($b1) != length($b2) )
		|| ( length($b1) != length($middle) ) )
	{
		return ( "", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 );

#die("This alignment has the 3 lines of different lengths.\n$b1\n$middle\n$b2\n");
	}

	if ( $dir eq "fw" )
	{
		my @sp1       = split( //, $b1 );
		my @sp_middle = split( //, $middle );
		my @sp2       = split( //, $b2 );
		my $i;
		for ( $i = 0 ; $i < scalar(@$aref_snp_coor) ; $i++ )
		{
			if ( $which_is_ref eq "1" )
			{
				push @$aref_snp_allel,
				  substr( $b1, $aref_snp_coor->[$i], 1 ) . "_"
				  . substr( $b2, $aref_snp_coor->[$i], 1 );
			}
			elsif ( $which_is_ref eq "2" )
			{
				push @$aref_snp_allel,
				  substr( $b2, $aref_snp_coor->[$i], 1 ) . "_"
				  . substr( $b1, $aref_snp_coor->[$i], 1 );
			}
			else
			{
				## no idea - just put one or the other
				print
"WARNING - no idea which allele is ref as we dont have a reference - just picking one\n";
				push @$aref_snp_allel,
				  substr( $b1, $aref_snp_coor->[$i], 1 ) . "_"
				  . substr( $b2, $aref_snp_coor->[$i], 1 );
			}
		}
	}
	else    ##dir is rev
	{
		my @sp1       = split( //, $b1 );
		my @sp_middle = split( //, $middle );
		my @sp2       = split( //, $b2 );
		my $i;
		for ( $i = 0 ; $i < scalar(@$aref_snp_coor) ; $i++ )
		{
			if ( $which_is_ref eq "1" )
			{
				push @$aref_snp_allel,
				  rev_comp(
					$sp1[ $index_of_last_non_match - ( $aref_snp_coor->[$i] ) ]
				  )
				  . "_"
				  . rev_comp(
					$sp2[ $index_of_last_non_match - ( $aref_snp_coor->[$i] ) ]
				  );
			}
			elsif ( $which_is_ref eq "2" )
			{
				push @$aref_snp_allel,
				  rev_comp(
					$sp2[ $index_of_last_non_match - ( $aref_snp_coor->[$i] ) ]
				  )
				  . "_"
				  . rev_comp(
					$sp1[ $index_of_last_non_match - ( $aref_snp_coor->[$i] ) ]
				  );
			}
			else
			{
				print "WARNING - no idea which is ref\n";
				push @$aref_snp_allel,
				  rev_comp(
					$sp1[ $index_of_last_non_match - ( $aref_snp_coor->[$i] ) ]
				  )
				  . "_"
				  . rev_comp(
					$sp2[ $index_of_last_non_match - ( $aref_snp_coor->[$i] ) ]
				  );

			}
		}
	}
}

## assumes which_is_ref is either 1 or 2
## it is possible, if there is an indel at the start or end, for the indel allele returned to be incpomplete, as you need to add
## a base from the 5prime or 3prime flank.
sub get_indel_alleles
{
	###index_of_last_non_match is only used in the case dir=rev. branches 1,2 agree for some number of bases at the end
	### - this gives the index (in normal forward coords) of the first base of difference from the end.
	## eg if we had
	#  Br1:CGCCGT__TTTT
	#      *| |*|  ||||
	#  Br2:TG_CCTTGTTTT
	# then this should be 7. ie coord of last non |

	my (
		$aref_indel_indices, $aref_indel_end_indices,
		$aref_indel_allel,   $aref_indel_need_base_adding,
		$b1,                 $middle,
		$b2,                 $dir,
		$which_is_ref,       $index_of_last_non_match
	) = @_;

	##double check
	if (   ( length($b1) != length($b2) )
		|| ( length($b1) != length($middle) ) )
	{
		return ( "", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 );

#die("This alignment has the 3 lines of different lengths.\n$b1\n$middle\n$b2\n");
	}

	if ( $dir eq "fw" )
	{

		my @sp1       = split( //, $b1 );
		my @sp_middle = split( //, $middle );
		my @sp2       = split( //, $b2 );
		my $i;
		for ( $i = 0 ; $i < scalar(@$aref_indel_indices) ; $i++ )
		{
			my $flag_needs_extra_base = 0;
			my $ref                   = "";
			my $alt                   = "";
			my $k = $aref_indel_indices->[$i]; ##base before indel's first space

			if ( $k < -1 )
			{
				die("k should not be $k\n");
			}

			if ( $k == -1 )   #removed constraint that index of last nonmatch==0
			{
				$flag_needs_extra_base = 1;
				$k++;
			}

			if ( $which_is_ref eq "1" )
			{
				$ref = $ref . $sp1[$k];
				$alt = $alt . $sp2[$k];
			}
			elsif ( $which_is_ref eq "2" )
			{
				$ref = $ref . $sp2[$k];
				$alt = $alt . $sp1[$k];
			}
			else
			{
				$ref = $ref . $sp1[$k];
				$alt = $alt . $sp2[$k];

			}
			$k++;

			while ( $k <= $aref_indel_end_indices->[$i] )
			{
				if ( $which_is_ref eq "1" )
				{
					$ref = $ref . $sp1[$k];
					$alt = $alt . $sp2[$k];
				}
				elsif ( $which_is_ref eq "2" )
				{
					$ref = $ref . $sp2[$k];
					$alt = $alt . $sp1[$k];
				}
				else
				{
					$ref = $ref . $sp1[$k];
					$alt = $alt . $sp2[$k];
				}

				$k++;
			}
			$ref =~ s/_//g;
			$alt =~ s/_//g;
			push @$aref_indel_allel,            $ref . "_" . $alt;
			push @$aref_indel_need_base_adding, $flag_needs_extra_base;
		}
	}
	else    ##dir is rev
	{

		my @sp1       = split( //, $b1 );
		my @sp_middle = split( //, $middle );
		my @sp2       = split( //, $b2 );
		my $i;
		for ( $i = 0 ; $i < scalar(@$aref_indel_indices) ; $i++ )
		{
			my $flag_needs_extra_base = 0;
			my $ref                   = "";
			my $alt                   = "";
			my $k                     = $aref_indel_indices->[$i];

			if ( $which_is_ref == 1 )
			{

				if ( $k > scalar(@sp1) - 1 )
				{
					$flag_needs_extra_base = 1;
					$k--;
				}

				$ref = $ref . $sp1[$k];
				$alt = $alt . $sp2[$k];
			}
			else
			{

				if ( $k > scalar(@sp2) - 1 )
				{
					$flag_needs_extra_base = 1;
					$k--;
				}

				$ref = $ref . $sp2[$k];
				$alt = $alt . $sp1[$k];

			}

			$k--;

			while ( $k >= $aref_indel_end_indices->[$i] )
			{
				if ( $which_is_ref eq "1" )
				{
					$ref = $ref . $sp1[$k];
					$alt = $alt . $sp2[$k];
				}
				elsif ( $which_is_ref eq "2" )
				{
					$ref = $ref . $sp2[$k];
					$alt = $alt . $sp1[$k];
				}
				else
				{
					$ref = $ref . $sp1[$k];
					$alt = $alt . $sp2[$k];
				}
				$k--;
			}
			$ref =~ s/_//g;
			$alt =~ s/_//g;
			$ref = reverse($ref);
			$alt = reverse($alt);
			$ref = rev_comp($ref);
			$alt = rev_comp($alt);

			push @$aref_indel_allel,            $ref . "_" . $alt;
			push @$aref_indel_need_base_adding, $flag_needs_extra_base;

		}

	}

}

sub combine_all_filters
{
	my (
	    $href_var_name_to_covg_and_branch_filter,
	    $href_var_name_to_flank_mq_filter,
	    $href_pd_filter,
	    $href_var_name_to_combined_filtering_result,
	    $href_pop_classifier
	    ) = @_;

	foreach my $key ( keys %$href_var_name_to_flank_mq_filter )
	{

		#	if ($href_var_name_to_flank_mq_filter->{$key} eq "PASS")#
		#	{
		#		$href_var_name_to_combined_filtering_result->{$key}="PASS";
		#	}
		if ( $href_var_name_to_flank_mq_filter->{$key} ne "PASS" )
		{
			$href_var_name_to_combined_filtering_result->{$key} = "MAPQ";
		}
	}
	foreach my $key (keys %$href_pd_filter)
	{
	    if ( !exists $href_var_name_to_combined_filtering_result->{$key} )
	    {
		$href_var_name_to_combined_filtering_result->{$key} =
		    "PF_FAIL_PD_PARALOG" ;
	    }
	    else
	    {
		$href_var_name_to_combined_filtering_result->{$key} =
		    $href_var_name_to_combined_filtering_result->{$key}
		. ",PF_FAIL_PD_PARALOG";

	    }
	}
	foreach my $key ( keys %$href_pop_classifier )
	{
		my $reason = uc( $href_pop_classifier->{$key} );
		if ( $href_pop_classifier->{$key} ne "variant" )
		{
			if ( !exists $href_var_name_to_combined_filtering_result->{$key} )
			{
				$href_var_name_to_combined_filtering_result->{$key} =
				  "PF_FAIL_" . $reason;
			}
			else
			{
				$href_var_name_to_combined_filtering_result->{$key} =
				    $href_var_name_to_combined_filtering_result->{$key}
				  . ",PF_FAIL_"
				  . $reason;
			}
		}
	}

	foreach my $key ( keys %$href_var_name_to_flank_mq_filter )
	{
		if ( !exists $href_var_name_to_combined_filtering_result->{$key} )
		{
			$href_var_name_to_combined_filtering_result->{$key} = "PASS";
		}
	}

}

sub get_list_vars_with_cut_flanks
{
	my ( $file, $href ) = @_;

	open( FILE, $file ) || die();

	while (<FILE>)
	{
		my $line = $_;
		if ( ( $line =~ /\w*var_\d+_5p_flank/ ) && ( $line =~ /cut_at_1000/ ) )
		{
			if ( $line =~ /(\w*var_\d+)_5p_flank/ )
			{
				my $name = $1;
				if ($prefix ne "")
				{
				    $name = $prefix . "_" .$name;
				}
				$href->{$name} = 1;
			}
			else
			{
				die("programming error on $line");
			}
		}
	}
	close(FILE);
}

sub get_pop_filter_info
{
	my ( $file, $href, $href_conf ) = @_;
	open( FILE, $file )
	  || die(
"Cannot find the file containing output of classifer.parallel.ploidy_aware.R - you entered it as an argument, $file"
	  );
	while (<FILE>)
	{
		my $line = $_;
		chomp $line;
		my @sp = split( /\t/, $line );
		my $name = $sp[0];
		if ($prefix ne "")
		{
		    $name = $prefix . "_" .$name; 
		}
		if ($name =~ /^(\S+)\s+/)
		{
		    $name = $1;
		}

		$href->{$name}      = $sp[1];    ## classification
		$href_conf->{$name} = $sp[2];

	}
	close(FILE);
}

sub wrap_needleman
{

	my ( $bin, $callfile, $prefix, $outfile ) = @_;
	
#	my $printed=0;
	open( OUT, ">" . $outfile ) || die("Cannot open $outfile");

## to add in front of var names to make them globally unique, otherwise all files contain var_1, var_2, etc
	open( FILE, $callfile ) || die("Cnnot open $callfile");
	my %seq   = ();
	my $count = 0;
	my @seq1  = ();
	my @seq2  = ();

	my $printed_at_start_of_var = 0;
	while (<FILE>)
	{

		my $line = $_;
		my $var_name;
		if ( $line =~ /([^>]+)_branch\_(1|2)/ )
		{
			if ( $printed_at_start_of_var == 0 )
			{
				print OUT "\n\nSTART NEW VAR\n";
				$printed_at_start_of_var = 1;
			}
			elsif ( $printed_at_start_of_var == 1 )
			{
				$printed_at_start_of_var = 0;
			}

			$var_name = $1;
			if ($prefix ne "")
			{
			    $var_name = $prefix ."_".$var_name;
			}
			my $which_branch = $2;
			print OUT "$var_name branch $which_branch\n";
			my $a = <FILE>;
			chomp $a;
			$seq{$which_branch} = $a;
			$count++;
		}
		elsif ( $line =~ /branch\_(\d+)\_(1|2)/ )
		{
		    $var_name = "var_" . $1;
		    if ($prefix ne "")
		    {
			$var_name = $prefix."_".$var_name;
		    }
		    my $which = $2;

		    if ( $printed_at_start_of_var == 0 )
		    {
			print OUT "\n\nSTART NEW VAR\n";
			$printed_at_start_of_var = 1;
		    }
		    elsif ( $printed_at_start_of_var == 1 )
		    {
			$printed_at_start_of_var = 0;
		    }

		    my $which_branch=$which;
		    if ( $which eq "trusted" )
		    {
			$which_branch = 1;
		    }
		    elsif ( $which eq "variant" )
		    {
			$which_branch = 2;
		    }

		    print OUT "$var_name branch $which_branch\n";
		    my $a = <FILE>;
		    chomp $a;
		    
		    $seq{$which_branch} = $a;
		    $count++;
		}
		else
		{
		    #die("Unexpected format of line in callfile - contact Zam. Your offending line is :\n$line\n");
		}

		if ( $count == 2 )
		{

			#print $seq{1},"\n";
			#print $seq{2},"\n";

			@seq1 = split //, $seq{1};
			@seq2 = split //, $seq{2};

			#my $cmd1 = "$bin --zam $seq{1} $seq{2}";
			#my $ret1 = qx{$cmd1};
			my $ret1 = run_nw_alignment_handling_long_strings($seq{1},$seq{2}, $bin, $var_name);

			print OUT "FORWARD ALIGNMENT\n";
			my ( $count_snps_f, $count_indels_f ) = get_snp_indel_counts($ret1);
			print OUT "$ret1";

			if (   ( $count_indels_f > @seq1 / 3 )
				|| ( $count_indels_f > @seq2 / 3 ) )
			{
				print OUT "REVERSE ALIGNMENT\n";
				@seq2 = split //, rev_comp( $seq{2} );
				my $rseq2 = rev_comp( $seq{2} );
				#my $cmd2  = "$bin --zam $seq{1} $rseq2";
				#my $ret2  = qx{$cmd2};
				my $rvar_name = $var_name."_rev_alignment";
				my $ret2 = run_nw_alignment_handling_long_strings($seq{1},$rseq2, $bin, $rvar_name);
				my ( $count_snps_r, $count_indels_r ) =
				    get_snp_indel_counts($ret2);
				print OUT "$ret2";

				#print OUT "\n";
			}
			else
			{
				print OUT "NO REVERSE ALIGNMENT\n";
			}
			$count = 0;

		}
	}
	close(FILE);
	close(OUT);
}

sub run_nw_alignment_handling_long_strings
{
    my ($str1, $str2, $bin, $var_id) = @_;

    if (length($str1) + length($str2) < 5000)
    {
	my $cmd = "$bin --zam $str1 $str2";
	my $ret = qx{$cmd};
	return $ret;
    }
    else
    {
	my $alignment_dir;
	if ($outdir !~ /\/$/)
	{
	    $alignment_dir = $outdir."/tmp_alignments/";
	}
	else
	{
	    $alignment_dir = $outdir."tmp_alignments/";
	}
	if (!(-d $alignment_dir))
	{
	    my $c = "mkdir $alignment_dir";
	    qx{$c};
	}
	my $tmpfile = $alignment_dir."alignment_".$var_id;
	open(TMPAL, ">".$tmpfile)||die("Cannot make temp alignment file $tmpfile - permissions problem? Out of disk?\n");
	print TMPAL ">\n$str1\n>\n$str2\n";
	close(TMPAL);
	my $cmd = "$bin --zam --file $tmpfile";
	my $ret = qx{$cmd};
	return $ret;
    }
}
sub get_snp_indel_counts
{
	my ($str) = @_;
	if ( $str =~ /^.+\n.+\n.+\n(\d+) (\d+)/ )
	{
		return ( $1, $2 );
	}
	else
	{
		die("Bad forma of $str in get_snp_indel_counts");
	}
}


sub get_min_covg
{
    my ($aref) = @_;
    my $min = 999999999999999;
    my $i;
    for ($i=0; $i<scalar(@$aref); $i++)
    {
	if ($aref->[$i]<$min)
	{
	    $min = $aref->[$i];
	}
    }
    return $min;
}
