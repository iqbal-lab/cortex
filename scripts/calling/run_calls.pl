#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Long;
use IO::Tee;

###******** Description of run_calls.pl ##################################################################
# This script allows the user to specify a number of kmer values and cleaning thresholds, and then
# will build graphs, clean them, call with both the bubble caller and PD caller, make union callsets
# genotype all samples at the union of all called sites, make a VCF and remove duplicate sites from the VCF
##############################################################################################################



#### Dear User - this following line is the only one in this script you may want to edit
my $stampy_bin ="/path/to/stampy.py";    #### <<<<<<< can also use command line parameter
my $vcftools_dir="/path/to/vcftools/dir";### <<<< can also use command line param
#### ****************  no need to modify anything below this line

# However you may decide you want to increase this threshold
my $mapping_qual_thresh = 40;    #  - demand 5prime flank maps with quality>= 40

my $callingscript_dir;
my $analyse_variants_dir;
#my $script_dir;
my $cortex_dir;
my $isaac_bioinf_dir;


BEGIN
{
	use FindBin;
	$callingscript_dir = $FindBin::Bin;
	$cortex_dir = $callingscript_dir . '/../../';
	$analyse_variants_dir = $cortex_dir."/scripts/analyse_variants/";
	$isaac_bioinf_dir = $analyse_variants_dir."bioinf-perl/";

	push( @INC,
		$cortex_dir
		  . "/scripts/analyse_variants/perl_modules/Statistics-Descriptive-2.6",
	      $isaac_bioinf_dir."lib/"
	    );
}

my $proc_calls = $analyse_variants_dir."process_calls.pl";
if (!(-e $proc_calls))
{
    die("Cannot find process_calls script: $proc_calls\n");
}
my $make_union = $analyse_variants_dir."make_union_varset.pl";
if (!(-e $make_union))
{
    die("Cannot find make_union_varset script: $make_union\n");
}

#use lib $isaac_bioinf_dir;
use Descriptive;



##make sure there is forward slash on the end:
if ( $cortex_dir !~ /\/$/ )
{
	$cortex_dir = $cortex_dir . '/';
}

my (
        $first_kmer,           $last_kmer,      $kmer_step,
        $do_auto_cleaning,     $auto_below,     $auto_above,
        $do_user_spec_cleaning,$user_min_clean, $user_max_clean, $user_clean_step,
        $do_bc,                            $bc_col_args,  $bc_out_stub,
        $do_pd,                            $pd_col_args,  $pd_out_stub,
        $outdir,               $fastaq_index,
	$outvcf_filename_stub, 
        $number_of_colours,    $use_ref,
	                 $apply_filter_one_allele_must_be_ref,
	$apply_classif,  $prefix,
	$ploidy,                $require_one_allele_is_ref,
	$stampy_hash_stub,
        $outdir_binaries,      $outdir_calls,   $outdir_vcfs, 
        $qthresh,              $dups,           $homopol,
        $mem_height,           $mem_width, $binary_colourlist,
        $max_read_len,         $format, $max_var_len, $genome_size, $refbindir, $list_ref_fasta,
        $expt_type,            $do_union, $manual_override_cleaning,
        $build_per_sample_vcfs, $global_logfile,
);

#set defaults
$first_kmer=0;
$last_kmer=0;
$kmer_step=1;
$do_auto_cleaning="no";
$auto_below=0;
$auto_above=0;
$do_user_spec_cleaning="no";
$user_min_clean=0;
$user_max_clean=0;
$user_clean_step=1;
$do_bc="yes";
$bc_col_args="-1/-1";
$bc_out_stub="bubbles";
$do_pd="no";
$pd_col_args="";
$pd_out_stub="PD";
$outdir="cortex_results/";
$fastaq_index="";
$require_one_allele_is_ref           = "yes";
$prefix                              = "cortex";
$outvcf_filename_stub                = "default_vcfname";
#$samplenames                          = '';
#$number_of_colours                   = 0;
$use_ref                             = "yes";
$ploidy                              = 2;
$apply_filter_one_allele_must_be_ref = "unknown";
$apply_classif                       ='';
$qthresh                             =-1;
$dups                                ='';
$homopol                             =-1;
$binary_colourlist                   ="";
$mem_height=-1;
$mem_width=-1;
$max_read_len=-1;
$max_var_len = 10000;
$format="UNPSECIFIED";
$genome_size=0;
$refbindir="";
$list_ref_fasta = "nonexistent_nonsense";
$expt_type = "";
$do_union = "no";
$manual_override_cleaning="no";
$build_per_sample_vcfs="no";
$global_logfile = "default_logfile";
my $pooled_colour = -1;    #deprecated
my $help = '';    #default false




&GetOptions(
    'first_kmer:i'         => \$first_kmer,
    'last_kmer:i'          => \$last_kmer,
    'kmer_step:i'          => \$kmer_step,
    'auto_cleaning:s'   => \$do_auto_cleaning, 
    'auto_below:i'         => \$auto_below,
    'auto_above:i'         => \$auto_above,
    'user_cleaning:s'      => \$do_user_spec_cleaning,
    'user_min_clean:i'        => \$user_min_clean,
    'user_max_clean:i'        => \$user_max_clean,
    'user_clean_step:i'    => \$user_clean_step,
    'bc:s'                 => \$do_bc, # $bc_col_args,
    'pd:s'                 => \$do_pd, # $pd_col_args,
    'outdir:s'             => \$outdir,   
    'outvcf:s'             => \$outvcf_filename_stub,
    'use_ref:s'             => \$use_ref,
    'ploidy:i'             => \$ploidy,
    'require_one_allele_is_ref' => \$apply_filter_one_allele_must_be_ref,
    'apply_pop_classifier'   => \$apply_classif,
    'prefix:s'             => \$prefix,            ## this will prefix any var name
    'stampy_hash:s'        => \$stampy_hash_stub, #stampy creates blah.sthash and blah.stidx. You should enter --stampy_hash blah
    'stampy_bin:s'         => \$stampy_bin,    ## must be 1 or 2
    'fastaq_index:s'        => \$fastaq_index,
    'qthresh:i'            => \$qthresh,
    'dups'                   => \$dups,
    'homopol:i'            => \$homopol,
    'mem_height:i'         => \$mem_height,
    'mem_width:i'          => \$mem_width,
    'colourlist:s'         => \$binary_colourlist,
    'max_read_len:i'       => \$max_read_len,
    'format:s'             =>\$format,
    'max_var_len:i'        => \$max_var_len,
    'help'                   => \$help,
    'genome_size:i'        => \$genome_size,
    'refbindir:s'          =>\$refbindir,
    'expt_type:s'          =>\$expt_type,
    'list_ref_fasta:s'     =>\$list_ref_fasta,
    'vcftools_dir:s'       =>\$vcftools_dir,
    'manual_override_cleaning:s' => \$manual_override_cleaning,
    'do_union:s'         =>\$do_union,
    'build_per_sample_vcfs:s' =>\$build_per_sample_vcfs,
    'logfile:s'           => \$global_logfile,
);

if ($help)
{
	print "\n\n";
	print "Usage:\n:\n";
	print "--first_kmer\t\t\t\tThis script allows you to ru across a range of kmer sizes. This is the lowest. It must be odd\n";
	print "--last_kmer\t\t\t\tIgnore this if you want to run for one kmer only\n";
	print "--kmer_step\t\t\t\tIf you run for many kmers, this is the increment.\n";
	print "--auto_cleaning\t\t\t\tTakes values yes or no. Default no. This looks at covg distribution and chooses a cleaning threshold\n";
	print "--auto_below\t\t\t\tYou can also ask it to run for, say 2 thresholds below the auto-chosen one. By default it wont do this\n";
	print "--auto_above\t\t\t\tYou can ask it to run for, say 3 threshold values above the auto-chosen one (will stay below the expected depth\n";
	print "--user_cleaning\t\t\t\tValid arguments are yes and no. Default is no. Make your own cleanig choices\n";
	print "--user_min_clean\t\t\t\tIf you want to try a range. Use this also if you only want to use one threshold.\n";
	print "--user_max_clean\t\t\t\tIf you want to try a range. Ignore this if you only want to use one threshold\n";
	print "--use_clean_step\t\t\t\tIncrement between cleaning thresholds.\n";
	print "--bc\t\t\t\tMake Bubble Calls. You must enter yes or no. Default (if you don't use --bc) is no.\n";
	print "--pd\t\t\t\tMake Path Divergence Calls. You must enter yes or no. Default (if you don't use --bc) is no.\n";
	print "--outdir\t\t\t\tOutput directory. Everything will go into dsubdirectories of this directory\n";
	print "--outvcf\t\t\t\tVCFs generated will have names that start with the text you enter here\n";
	print "--use_ref\t\t\t\tValid values are yes and no, depending on whether you want to include a reference genome when calling.\n";
	print "--ploidy\t\t\t\tMust be 1 or 2.\n";
	print "--require_one_allele_is_ref\t\t\t\tyes or no. Highly recommended if you want a VCF with ref chromosome positions\n";
	print "--prefix\t\t\t\tIf you want your variant calls to have names with a specific prefix, use this\n";
	print "--stampy_hash\t\t\t\tMANDATORY. Build a hash of your reference genome, and specify here the path to it. if stampy makes blah.stdidx etc then specify blah.\n";
	print "--stampy_bin\t\t\t\tSpecify the path to your Stampy binary. Or manually edit this at the top of the file (it's marked out for you)\n";
	print "--fastaq_index\t\t\t\tMANDATORY. File has format SAMPLE_NAME\tse_list\tpe_list1\tpe_list2. One line per sample\n";
	print "--qthresh\t\t\t\tIf you want Cortex to use a quality score threshold, speify it here\n";
	print "--dups\t\t\t\tIf you want Cortex to remove PCR duplicates, specify this flag (no arguments, just --dups)\n";
	print "--homopol\t\t\t\tIf you want to cut homopolymers, threshold specified here\n";
	print "--mem_height\t\t\t\tFor Cortex\n";
	print "--mem_width\t\t\t\tFor Cortex\n";
	print "--max_read_len\t\t\t\tMax read length\n";
	print "--format\t\t\t\tFASTA or FASTQ\n";
	print "--max_var_len\t\t\t\tSee Cortex manual - max var length to look for\n";
	print "--genome_size\t\t\t\tGenome length in base pairs - needed for genotyping\n";
	print "--refbindir\t\t\t\tDiorectry containing binaries built of the reference at all the kmers you want to use. \n";
	print "\t\t\t\t\tThe binary filename should contain the kmer value, eg refbinary.k31.ctx\n";
	print "--expt_type\t\t\t\tAs in Cortex input\n";
	print "--list_ref_fasta\t\t\t\tFile listing the fasta files (one per chromosome) for the reference. Needed for the PD caller\n";
	print "--vcftools_dir\t\t\t\tVCFtools is used to generate VCFs - mandatory to either specify this on cmd-line, or manually edit the path at the top of this script\n";
	print "--do_union\t\t\t\tHaving made per-sample callsets (per kmer and cleaning), should we combine all calls into a union set, and genotype all samples? Valid values are yes and no. Default is no.\n";
	print "--manual_override_cleaning\t\t\t\tYou can specify specific thresholds for specific samples by giving a file here, each line has three (tab sep) columns: sample name, kmer, and comma-separated thresholds\nDon't use this unless you know what you are doing\n";
	print "--build_per_sample_vcfs\t\t\t\tThis script repeatedly runs Cortex BC and PD callers, calling on each sample separately, and then by default builds one pair (raw/decomp) of VCFs for the union set. If in addition you want VCFs built for each callset, enter \"yes\" here. In general, do not do this, it is very slow.\n";
	print "--help\t\t\t\tprints this\n";
	exit();
}





my %k_to_refbin=();
### checks
run_checks();


## sort out loggin
my $tee = IO::Tee->new(">> $global_logfile", \*STDOUT);


print $tee "*********** Command-line used was : *********************\n";
print $tee "perl ".$0." ".join(" ",@ARGV)."\n";
print $tee "*********************************************************\n";




$outdir_binaries=$outdir."binaries/";
$outdir_calls=$outdir."calls/";
$outdir_vcfs=$outdir."vcfs/";

## Now the process we will follow is
## 1. Build uncleaned binaries and dump covg distributions
## 2. Clean
## 3. Call variants
## 4. Make a union callset
## 5. Genotype all samples on the union callset
## 6. Use pop classifier  <<<<< not incorporated yet
## 7. Build a VCF
## 8. Clean up the vCF (sort it, remove dup lines, remove lines where ref-allele doesnt match reference).


## basic setup
my @samples = ();
get_sample_names($fastaq_index, \@samples);

my $ref_name = "none";
if ($use_ref eq "yes")
{
    $ref_name = "REFERENCE";
}

my @kmers=();
get_kmers(\@kmers, $first_kmer, $last_kmer, $kmer_step);



##################################
## 1. Build uncleaned binaries
##################################

print $tee "********************************************\n";
print $tee "Build uncleaned binaries:\n";
print $tee "********************************************\n";

my %sample_to_uncleaned=(); #samplename --> kmer --> uncleaned bin name
my %sample_to_uncleaned_log=(); #samplename --> kmer --> log file from building uncleaned bin 
my %sample_to_uncleaned_covg_distrib=(); #samplename --> kmer --> covg distribution file
my %sample_to_min_cleaning_thresh=(); #sample -> kmer -> min clean thresh

if ($fastaq_index ne "")
{
    build_all_unclean_binaries($fastaq_index, $outdir_binaries, 
			       \%sample_to_uncleaned, \%sample_to_uncleaned_log,
			       \%sample_to_uncleaned_covg_distrib,
			       $cortex_dir, $qthresh, $dups, $homopol, $mem_height, $mem_width, $max_read_len, $format);
}


my $uncleaned_stats_log = $outdir_binaries."uncleaned/UNCLEANED_BINS_STATS";
print_build_stats(\%sample_to_uncleaned_log, $uncleaned_stats_log); 


##################################
## 2. Clean
##################################
print $tee "********************************************\n";
print $tee "Clean binaries\n";
print $tee "********************************************\n";
my %sample_to_cleaned_bin=(); #sample -> kmer -> cleaning -> $binary
if ( ($do_auto_cleaning eq "yes") || ($do_user_spec_cleaning eq "yes") || ($manual_override_cleaning ne "no") )
{
    build_all_cleaned_binaries(\%sample_to_uncleaned, \%sample_to_uncleaned_log,
			       \%sample_to_uncleaned_covg_distrib,
			       \%sample_to_cleaned_bin,
			       $outdir_binaries, 
			       $cortex_dir, $mem_height, $mem_width, $max_var_len,
			       $do_auto_cleaning, $auto_below, $auto_above,
			       $do_user_spec_cleaning, $user_min_clean, $user_max_clean, $user_clean_step, $genome_size,
	                       $manual_override_cleaning);
}
else
{
    print $tee "No cleaning specified, so will not do any - will call variants on uncleaned binaries\n";

    foreach my $k (@kmers)
    {
	foreach my $sam (@samples)
	{
	   # if ($sam eq $ref_name)
	   # {
#		next;
#	    }
	    $sample_to_cleaned_bin{$sam}{$k}{0}=$sample_to_uncleaned{$sam}{$k};
	}
    }


}




##################################
## 3. Call variants
##################################


## v1 of this script is for calling independently on all non-ref colours, make a union callset, and then genotype all samples on the union set

## How do we do this? Suppose we are doing 3  levels of cleaning per sample, say cl1, cl2, cl3 (these will be different for each sample, but I just want to be
## able to refer to the first, second and third level). 

if ( ($do_bc ne "yes") && ($do_pd ne "yes") )
{
    print $tee "Binaries are built, but you have not specified to run any calling, so will halt now. Good night.\n";
    exit(0);
}

print $tee "********************************************\n";
print $tee "Call variants\n";
print $tee "********************************************\n";
my $dir_for_per_sample_calls=$outdir_calls."per_sample_callsets/";
print $tee "per sample call dir us $dir_for_per_sample_calls\n";
if (!(-d $dir_for_per_sample_calls))
{
    my $c = "mkdir -p $dir_for_per_sample_calls";
    qx{$c};
}


my %sample_to_bc_callfile=();#sample ->kmer ->cleaning ->callfile
my %sample_to_pd_callfile=();

foreach my $k (@kmers)
{
    my $ctx_bin = get_right_binary($k, $cortex_dir,2);
    foreach my $sam (@samples)
    {
	#if ($sam eq $ref_name)
	#{
	#    next;
	#}

	foreach my $cleaning (keys %{$sample_to_cleaned_bin{$sam}{$k}})
	{
	    my $uniqid = "sample_".$sam."_kmer".$k."_cleaning".$cleaning;
	    my $colour_list = make_2sample_filelist($ref_name.".".$uniqid, $sam.".".$uniqid, $k_to_refbin{$k}, $sample_to_cleaned_bin{$sam}{$k}{$cleaning}, $uniqid );
	    ## load reference binary and make calls. 
	    my $cmd = $ctx_bin." --kmer_size $k --mem_height $mem_height --mem_width $mem_width --ref_colour 0 --colour_list $colour_list  --print_colour_coverages ";
	    print $tee "Load reference $k_to_refbin{$k} in colour 0, and sample ";
	    print $tee $sample_to_cleaned_bin{$sam}{$k}{$cleaning};
	    print $tee " into colour 1\n";

	    my $bubble_output = $dir_for_per_sample_calls.$sam."_bubbles_k".$k."_clean".$cleaning.".calling_on_this_sample_only";
	    my $pd_output = $dir_for_per_sample_calls.$sam."_pd_k".$k."_clean".$cleaning.".calling_on_this_sample_only";

	    ##these variables needed to decide whether to run var valling
	    my $bc_already_done="no";
	    my $pd_already_done = "no";
	    if ($do_bc eq "no")
	    {
		$bc_already_done="yes";
	    }
	    if ($do_pd eq "no")
	    {
		$pd_already_done="yes";
	    }

	    if ( ($do_bc eq "yes") && (!(-e $bubble_output)) )## we want to do it and not already done
	    {
		$cmd = $cmd." --detect_bubbles1 -1/-1 --output_bubbles1 $bubble_output --exclude_ref_bubbles --max_var_len $max_var_len ";
		$sample_to_bc_callfile{$sam}{$k}{$cleaning} =  $bubble_output ;
	    }
	    elsif  ( ($do_bc eq "yes") && (-e $bubble_output) )## we want t do it but it has already been done
	    {
		print $tee "Bubble calling already seems to have been done - will not rerun\n";
		$bc_already_done = "yes";
		$sample_to_bc_callfile{$sam}{$k}{$cleaning} =  $bubble_output ;
	    }
	    if ( ($do_pd eq "yes") && (!(-e $pd_output."_pd_calls")) )
	    {
		$cmd = $cmd." --path_divergence_caller 1 --path_divergence_caller_output $pd_output --list_ref_fasta $list_ref_fasta --max_var_len $max_var_len ";
		$sample_to_pd_callfile{$sam}{$k}{$cleaning} =  $pd_output."_pd_calls"
	    }
	    elsif ( ($do_pd eq "yes") && (-e $pd_output."_pd_calls"))
	    {
		print $tee "Path divergence calling already seems to have been done - will not rerun\n";
		$pd_already_done = "yes";
		$sample_to_pd_callfile{$sam}{$k}{$cleaning} =  $pd_output."_pd_calls"
	    }

	    my $log = $dir_for_per_sample_calls.$sam."_varcalling_log_k".$k."_clean".$cleaning.".calling_on_this_sample_only.log";
	    if (-e $log)
	    {
		$log = $log.".again";
	    }

	    if (!( ($bc_already_done eq "yes") && ($pd_already_done eq "yes")))##not all done
		{
		    $cmd = $cmd." > $log 2>&1";
		    print $tee "$cmd\n";
		    my $ret = qx{$cmd};
		    print $tee "$ret\n";
		}


	    ## Add code to check log file for errors here:
		
	}
    }
}




##################################
## 4. Make union callsets
##################################
my %vcfs_needing_post_processing=();
if ($do_union eq "yes")
{
    print $tee "********************************************\n";
    print $tee "Make union variant callsets\n";
    print $tee "********************************************\n";
    
    
    my $tmpdir = $outdir."tmp_filelists";
    
    if (!(-d  $tmpdir))
    {
	my $cmd1 = "mkdir -p $tmpdir";
	qx{$cmd1};
    }
    
    my $union_of_bc_callsets = $outdir_calls."union_all_bc_callsets";
    my $max_read_len_bc_union=0;
    if ($do_bc eq "yes")
    {
	my $bc_call_list  = $tmpdir."/list_bc_callfiles";
	open(BCCALL, ">".$bc_call_list)||die();
	foreach my $k (@kmers)
	{
	    foreach my $sam (@samples)
	    {
		foreach my $cleaning (keys %{$sample_to_cleaned_bin{$sam}{$k}})
		{
		    print BCCALL $sample_to_bc_callfile{$sam}{$k}{$cleaning};
		    print BCCALL "\n";
		}
	    }
	}
	close(BCCALL);
	my $bc_cmd = "perl $make_union --filelist $bc_call_list --varname_stub UNION_BC > $union_of_bc_callsets 2>&1";
	if (! -e($union_of_bc_callsets))
	{
	    print $tee "$bc_cmd\n";
	    my $bc_ret = qx{$bc_cmd};
	}
	$max_read_len_bc_union = get_max_read_len_of_fasta($union_of_bc_callsets)+$last_kmer+10;
    }
    
    my $union_of_pd_callsets = $outdir_calls."union_all_pd_callsets";
    my $max_read_len_pd_union=0;
    if ($do_pd eq "yes")
    {
	
	my $pd_call_list  = $tmpdir."/list_pd_callfiles";
	open(PDCALL, ">".$pd_call_list)||die();
	foreach my $k (@kmers)
	{
	    foreach my $sam (@samples)
	    {
		foreach my $cleaning (keys %{$sample_to_cleaned_bin{$sam}->{$k}})
		{
		    print PDCALL $sample_to_pd_callfile{$sam}{$k}{$cleaning};
		    print PDCALL "\n";
		}
	    }
	}
	close(PDCALL);
	
	my $pd_cmd = "perl $make_union --filelist $pd_call_list --varname_stub UNION_PD > $union_of_pd_callsets 2>&1";
	if (!(-e $union_of_pd_callsets))
	{
	    print $tee "$pd_cmd\n";
	    my $pd_ret = qx{$pd_cmd};
	}
	$max_read_len_pd_union = get_max_read_len_of_fasta($union_of_pd_callsets)+$last_kmer+10;
    }
    
    
    
    
    if ( ($expt_type ne "") && ($genome_size != 0) )
    {
	
######################################################################################################
## 5. Genotype te union callset on ALL samples, using multicoloured graph with lowest k and lowest cleaning
##    Note it has to be the lowest kmer, otherwise some calls will have flanks and branches which are shorter than our kmer
######################################################################################################
	
	print $tee "********************************************\n";
	print $tee "Genotype the union callset\n";
	print $tee "********************************************\n";
	
	
	my @ordered_list_binaries=();
	push @ordered_list_binaries, $k_to_refbin{$first_kmer};
	foreach my $sam (@samples)
	{	
	    my $min =999999999;
	    my $max=0;

	    foreach my $c (keys %{$sample_to_cleaned_bin{$sam}{$first_kmer}})
	    {
		#print $tee "For $sam, $first_kmer, cleanign thresh $c\n";
		#if ($c<$min)
		#{
		#    $min=$c;
		#}
		if ($c>$max)
		{
		    $max=$c;
		}
	    }
	    #push @ordered_list_binaries, $sample_to_cleaned_bin{$sam}{$first_kmer}{$min};
	    push @ordered_list_binaries, $sample_to_cleaned_bin{$sam}{$first_kmer}{$max};
	}
	my $multicolour_list = make_multicol_filelist(\@ordered_list_binaries);
## genotype all of the calls, called with all different kmers, on the max kmer
	my $gt_bc_out="";
	my $multicol_ctx_bin = get_right_binary($first_kmer, $cortex_dir, $number_of_colours); 
	if ($do_bc eq "yes")
	{
	    
	    $gt_bc_out = $outdir_calls.basename($union_of_bc_callsets).".genotyped";
	    print $tee "GT bc out is $gt_bc_out\n";
	    my $gt_bc_log = $gt_bc_out.".log";
	    
	    if (!(-e $gt_bc_out))
	    {
		my $gt_bc_cmd = $multicol_ctx_bin." --colour_list $multicolour_list  --kmer_size $first_kmer --mem_height $mem_height --mem_width $mem_width --experiment_type $expt_type --genome_size $genome_size --ref_colour 0 --gt $union_of_bc_callsets,$gt_bc_out,BC --max_read_len $max_read_len_bc_union  --print_colour_coverages --experiment_type $expt_type --genome_size $genome_size > $gt_bc_log 2>&1  ";
		print $tee "$gt_bc_cmd\n";
		my $gt_bc_ret = qx{$gt_bc_cmd};
		print $tee "$gt_bc_ret\n";
	    }
	    else
	    {
		print $tee "$gt_bc_out already exists so no need to genotype the union of BC calls\n";
	    }
	}
	my $gt_pd_out="";
	if ($do_pd eq "yes")
	{	
	    $gt_pd_out = $outdir_calls.basename($union_of_pd_callsets).".genotyped";
	    my $gt_pd_log = $gt_pd_out."_pd_calls.log";
	    
	    if (!(-e $gt_pd_out))
	    {
		my $gt_pd_cmd = $multicol_ctx_bin." --colour_list $multicolour_list --kmer_size $first_kmer --mem_height $mem_height --mem_width $mem_width --experiment_type $expt_type --genome_size $genome_size --ref_colour 0  --gt $union_of_pd_callsets,$gt_pd_out,PD --max_read_len $max_read_len_pd_union  --print_colour_coverages  --experiment_type $expt_type --genome_size $genome_size > $gt_pd_log 2>&1  ";
		print $tee "$gt_pd_cmd\n";
		my $gt_pd_ret = qx{$gt_pd_cmd};
		print $tee "$gt_pd_ret\n";
	    }
	    else
	    {
		print $tee "$gt_pd_out already exists so no need to genotype the union of PD calls\n";
	    }
	}
	
	
	
######################################################################################################
## 6. Make one VCF for BC and one for PD
######################################################################################################
	
	
	print $tee "********************************************\n";
	print $tee "Build union VCFs\n";
	print $tee "********************************************\n";
	
	
	if (!(-d $outdir_vcfs))
	{
	    my $c = "mkdir -p $outdir_vcfs";
	    qx{$c};
	}
	## build the VCFs of union calls
	
	print $tee "\n******   Build initial (uncleaned) VCFs of union callsets  ***\n";
	
	

	build_vcfs($gt_bc_out, "union_BC_calls", $number_of_colours, $outdir_vcfs."output_proc_union_bc", $outdir_vcfs, \%vcfs_needing_post_processing);
	build_vcfs($gt_pd_out, "union_PD_calls", $number_of_colours, $outdir_vcfs."output_proc_union_pd", $outdir_vcfs, \%vcfs_needing_post_processing);
	
	

	if ($build_per_sample_vcfs eq "yes")
	{
	    print $tee "\n******   Build VCFs for each of the individual callsets   ***\n";
	    build_per_sample_vcfs($outdir_vcfs.'/'."per_sample_vcfs", \%vcfs_needing_post_processing);
	}

######################################################################################################
## 7. Post-process the VCFS - make a union, remove dupliate lines, remove variants where 
##    the ref allele does not match the reference. Check each line has right number of fields.
######################################################################################################

## Now we have a bunch of VCFs, and we need to remove duplicate lines, sort them, 
	print $tee "\n****** Clean the union VCFs  ***\n";
	clean_all_vcfs(\%vcfs_needing_post_processing);
    }
}

else  ##### if not making a union set
{

    if ($build_per_sample_vcfs eq "yes")
    {
	## build the VCFs for each of the individual callsets
	print $tee "\n******   Build VCFs for each of the individual callsets   ***\n";
	build_per_sample_vcfs($outdir_vcfs.'/'."per_sample_vcfs", \%vcfs_needing_post_processing);
    }
    clean_all_vcfs(\%vcfs_needing_post_processing);
}

















#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################
#######################################################################################################


sub clean_all_vcfs
{
    my ($href) = @_;
    foreach my $k (keys %$href)
    {
	clean_vcf($k);
    }
}


sub clean_vcf
{
    my ($file) = @_;

    my $type="";
    if ($file=~ /raw/)
    {
	$type="raw";
    }
    elsif ($file =~ /decomp/)
    {
	$type="decomp";
    }
    else
    {
	die("Trying to clean a VCF that is neither raw nor decomp - $file\n");
    }

    #my $tmpdir = $outdir_vcfs."tmpdir";
    #if (!(-d $tmpdir))
    #{#
	#my $c = "mkdir -p $tmpdir";
	#qx{$c};
    #}

    my $ref_fa = get_ref_fasta($list_ref_fasta); 
    #my $bname = basename($file);
    #my $cleaned_file = $tmpdir.'/'."$bname".".cleaned";
    my $final_file = $file.".clean.sorted";
    my $cmd1 = $isaac_bioinf_dir."vcf_scripts/vcf_align.pl --remove_ref_mismatch --tag PV LEFT $file $ref_fa  | $vcftools_dir/perl/vcf-sort | $isaac_bioinf_dir"."vcf_scripts/vcf_remove_dupes.pl | grep -v vcf_remove_dupes.pl  > $final_file";
    print $tee "$cmd1\n";
    my $ret1 = qx{$cmd1};
    print $tee "$ret1\n";


    ##now sort
    #
    #my $sortcmd = "(cat $cleaned_file  | head -100 | grep ^#; cat $cleaned_file  | grep -v ^# | sort -k1,1d -k2,2n;) > $final_file";
    #qx{$sortcmd};


    ## Now remove variants where the ref allele does not match the ref allele at that position in the ref (mapping error)    
    ## <<<< fill in


}


sub build_per_sample_vcfs
{
    my ($dir, $href) = @_;

    if ($dir !~ /\/$/)
    {
	$dir = $dir.'/';
    }
    if (!(-d $dir))
    {
	my $c = "mkdir -p $dir";
	qx{$c};
    }

    foreach my $k (@kmers)
    {
	foreach my $sam (@samples)
	{
	    my $samdir = $dir.$sam;
	    if (!(-e $samdir))
	    {
		my $s = "mkdir -p $samdir";
		qx{$s};
	    }
	    
	    foreach my $cleaning (keys %{$sample_to_cleaned_bin{$sam}{$k}})
	    {
		if ($do_bc eq "yes")
		{
		    my $bc_file =  $sample_to_bc_callfile{$sam}{$k}{$cleaning};
		    my $stub = "BC_sample_".$sam."_kmer".$k."_cleaning".$cleaning;
		    my $this_log =$samdir. $stub.".log";
		    build_vcfs($bc_file, $stub, 2, $this_log, $samdir, $href);
		}
		if ($do_pd eq "yes")
		{
		    my $pd_file =  $sample_to_pd_callfile{$sam}{$k}{$cleaning};
		    my $stub = "PD_sample_".$sam."_kmer".$k."_cleaning".$cleaning;
		    my $this_log =$samdir. $stub.".log";
		    build_vcfs($pd_file, $stub, 2, $this_log, $samdir, $href);
		}

	    }
	}
    }

}


sub build_vcfs
{
    my ($file, $string, $num, $log, $directory, $href_store_vcf_names) = @_;
    if (!(-e $file))
    {
	return;
    }
    if ($directory !~ /\/$/)
    {
	$directory = $directory.'/';
    }
    my $w = "wc -l $file";
    my $rw = qx{$w};
    if ($rw =~ /^0\s+$file/)
    {
	return;
    }

    my $colournames=$outdir_vcfs.'/'."SAMPLES";
    if (!(-e $colournames))
    {
	open(COL, ">".$colournames)||die();
	if ($use_ref eq "yes")
	{
	    print COL "REF\n";
	}
	my $j;
	for ($j=0; $j<scalar(@samples); $j++)
	{
	    print COL $samples[$j];
	    print COL "\n";
	}
	close(COL);
    }
    my $cmd = "perl ".$proc_calls." --callfile $file --outdir $directory --outvcf $string --samplename_list $colournames --num_cols $num  --ploidy $ploidy --stampy_hash $stampy_hash_stub --stampy_bin $stampy_bin ";
    if ($use_ref eq "yes")
    {
	$cmd = $cmd." --refcol 0 --require_one_allele_is_ref yes"; 
    }
    $cmd = $cmd."  > $log 2>&1" ;


    if ( (!(-e $directory.$string.".decomp.vcf"))  || (!(-e $directory.$string.".raw.vcf")) )
    {
	print $tee "$cmd\n";
	my $ret = qx{$cmd};
	print $tee "$ret\n";
    }
    else
    {
	print $tee $directory.$string.".decomp.vcf";
	print $tee " and ";
	print $tee $directory.$string.".raw.vcf";
	print $tee " both exist, so no need to rebuild\n";
    }
    $href_store_vcf_names->{$directory.$string.".decomp.vcf"}=1;
    $href_store_vcf_names->{$directory.$string.".raw.vcf"}=1;
}




sub make_multicol_filelist
{
    my ($aref_bins) = @_;

    my $tmpdir = $outdir."tmp_filelists";
    if (!(-d  $tmpdir))
    {
	my $cmd1 = "mkdir $tmpdir";
	qx{$cmd1};
    }
    my $colourlist = $tmpdir."/tmp_multicol_col_list";
    open(TMP, ">".$colourlist)||die();
    my $i;
    for ($i=0; $i<scalar(@$aref_bins); $i++)
    {
	print TMP "$tmpdir/colour$i"."_filelist\n";
    }
    close(TMP);
    for ($i=0; $i<scalar(@$aref_bins); $i++)
    {    
	open(TMP, ">"."$tmpdir/colour$i"."_filelist")||die();
	print TMP $aref_bins->[$i];
	print TMP "\n";
	close(TMP);
    }
    return $colourlist;
}

sub make_2sample_filelist
{

    my ($str1, $str2, $bin1, $bin2, $uniq_id) = @_;

    my $tmpdir = $outdir."tmp_filelists";

    if (!(-d  $tmpdir))
    {
	my $cmd1 = "mkdir -p $tmpdir";
	qx{$cmd1};
    }
    my $colourlist = $tmpdir."/tmp_col_list_".$uniq_id;
    open(TMP, ">".$colourlist)||die("Cannot open $colourlist");
    print TMP "$tmpdir/$str1\n$tmpdir/$str2\n";
    close(TMP);
    
    open(TMP, ">".$tmpdir."/$str1")||die();
    print TMP "$bin1\n";
    close(TMP);
    open(TMP, ">".$tmpdir."/$str2")||die();
    print TMP "$bin2\n";
    close(TMP);
    return $colourlist;
}








sub build_all_cleaned_binaries
{
    my ($href_sample_to_uncleaned, 
	$href_sample_to_uncleaned_log,
	$href_sample_to_uncleaned_covg_distrib,
	$href_sample_to_cleaned,
	$outdir_binaries, #,$samplenames
	$cortex_dir, $mem_height, $mem_width, $max_var_len,
	$do_auto_cleaning, $auto_below, $auto_above,
	$do_user_spec_cleaning, $user_min_clean, $user_max_clean, $user_clean_step, $g_size,
	$manual_override_cleaning_file) = @_;



    foreach my $sample (keys %$href_sample_to_uncleaned)
    {
	foreach my $k (@kmers)
	{

	    my @clean_threshes=();


	    get_cleaning_thresholds(
		$sample, $k,
		$href_sample_to_uncleaned_log, $href_sample_to_uncleaned_covg_distrib,
		$do_auto_cleaning, $auto_below, $auto_above,
		$do_user_spec_cleaning, $user_min_clean, $user_max_clean, $user_clean_step,
		\@clean_threshes, $g_size, $manual_override_cleaning_file);
	    

	    my $min = 9999999;
	    foreach my $c (@clean_threshes)
	    {
		if ($c < $min)
		{
		    $min=$c;
		}
		build_clean_binary($sample, $k, $c,
				   $outdir_binaries, 
				   $href_sample_to_uncleaned, 
				   $cortex_dir, $mem_height, $mem_width,
                  		    $href_sample_to_cleaned);
	    }
	    $sample_to_min_cleaning_thresh{$sample}{$k}=$min;
	}
    }    
}


sub build_clean_binary
{
    my ($sample, $kmer, $clean_thresh, $outdir_bins, $hash_sample_to_uncleaned, $cortex_dir, $height, $width, $href_sam_to_cleaned_bin) = @_;
    
    if ($outdir_bins !~ /\/$/)
    {
	$outdir_bins = $outdir_bins.'/';
    }
    my $outdr = $outdir_bins."cleaned/k".$kmer;
    if (!(-d $outdr))
    {
	my $cmd1= "mkdir -p $outdr";
	qx{$cmd1};
    }
    my $uncleaned = $hash_sample_to_uncleaned->{$sample}->{$kmer};
    my $uncleaned_bname = basename($uncleaned);
    
    my $ctx = $outdr.'/'.$uncleaned_bname;
    $ctx =~ s/.unclean//;
    $ctx =~ s/.ctx//;
    $ctx = $ctx."cleaned_".$clean_thresh.".ctx";
    my $log = $ctx.".log";
    $href_sam_to_cleaned_bin->{$sample}->{$kmer}->{$clean_thresh}=$ctx;

    ## skip this if binary already built
    if (-e $ctx)
    {
	if (!(-e $log))
	{
	    print $tee "Binary $ctx exists, so will not rebuild, but the log file is missing. Carrying on nevertheless.\n";
	}
	else
	{
	    print $tee "Binary $ctx already exists, so will not rebuild\n";
	}
	return;
    }


    my $cortex_binary = get_right_binary($kmer, $cortex_dir,1 );##one colour
    my $cmd2 = $cortex_binary." --kmer_size $kmer --mem_height $height --mem_width $width --dump_binary $ctx --remove_low_coverage_supernodes $clean_thresh --multicolour_bin $uncleaned > $log 2>&1";
    my $ret2 = qx{$cmd2};
    print $tee "$cmd2\n$ret2\n";
    $sample_to_cleaned_bin{$sample}{$kmer}{$clean_thresh}=$ctx;
    print $tee "add $sample  $kmer $clean_thresh = $ctx\n";
    if (!(-e $ctx))
    {
	die("Unable to build $ctx\n");
    }
}
sub get_cleaning_thresholds
{
    my ($sampl, $km, $href_log, $href_covg,
	$do_auto, $auto_below, $auto_above, 
	$do_user, $user_min,   $user_max, $user_step,
	$aref, $g_size, $manual_override_file) = @_;

    if ($do_auto eq "yes")
    {
	get_auto_thresholds($sampl, $km, $auto_below, $auto_above, $aref, $href_log, $href_covg, $g_size);
    }
    elsif ($do_user eq "yes")
    {
	my $i;
	for ($i=$user_min; $i<= $user_max; $i+=$user_step)
	{
	    push @$aref, $i;
	}
    }

    ##either way, you can also manually override and add some more thresholds
    get_manually_specified_thresholds_for_this_sample($sampl, $km, $manual_override_file, $aref);

}

sub get_manually_specified_thresholds_for_this_sample
{
    my ($sample_name, $khmer, $file, $array_ref) = @_;
    if ($file eq "no")
    {
	return;
    }

    my %tmp=();## to make sure we don't double-add thresholds that are already there
    foreach my $c (@$array_ref)
    {
	$tmp{$c}=1;
    }

    open(F, $file)||die("Cannot open file $file\n");
    while (<F>)
    {
	my $line = $_;
	chomp $line;
	my @sp = split(/\t/, $line);
	if(  ($sp[0] eq $sample_name) && ($sp[1] eq $khmer) )
	{
	    my @sp2 = split(/,/, $sp[2]);
	    foreach my $thresh (@sp2)
	    {
		if (!exists $tmp{$thresh})##if it isn't already on the list
		{
		    push @$array_ref, $thresh;
		}
	    }
	    last;
	}
    }
    close(F);
    
}

sub get_auto_thresholds
{
    my ($sample, $kmer, $num_below, $num_above, $aref_results, $href_logfiles, $href_covg, $genome_siz) = @_;
    
    my @distrib=();
    my $deff = get_expected_depth_and_cvg_distrib($href_logfiles->{$sample}->{$kmer}, $kmer, $genome_siz);
    my $auto_thresh  = get_cleaning_thresh_and_distrib($href_covg->{$sample}->{$kmer}, $deff,  \@distrib);
    print $tee "Automated cleaning chooses threshold $auto_thresh for sample $sample, kmer $kmer\nExpected depth of covg (D_eff) is $deff\n";
    if ($deff<$auto_thresh)
    {
	die("WARNING - the auto-estimator for cleaning has estimated a cleaning threshold higher than the expected covg. Look at the covg profile - is there something odd?\n");
    }
    ## ok, get results
    my $i;
    push @$aref_results, $auto_thresh;
    ## how many cleaning thresholds were you asked to generate below the auto one?

    ## work in steps which are 1/10 of the distance between the auto-choice, and the deff.
    my $step = int($deff) - $auto_thresh;
    for ($i=0; ($i<$num_below); $i++)
    {
	my $t= $auto_thresh-($i+1)*$step;
	if ($t>0)
	{
	    push @$aref_results, $t;
	    print $tee "Will also clean at threshold $t\n";
	}
    }

    ## how many cleaning thresholds were you asked to generate above the auto one?
    for ($i=0; ($i<$num_above); $i++)
    {
	my $t= $auto_thresh+($i+1)*$step;
	if ($t<$deff)
	{
	    push @$aref_results, $t;
	    print $tee "Will also clean at thresh $t\n";
	}
    }

}

sub get_kmers
{
    my ($aref_kmers, $first_k, $last_k, $step)=@_;
    if ( ($last_k - $first_k) % $step != 0)
    {
	die("You must be able to increment from first kmer, to last in an integer number of steps");
    }
    my $i;
    for ($i=$first_k; $i<=$last_k; $i+=$step)
    {
	push @$aref_kmers, $i;
    }
}

sub get_sample_names
{
    my ($file, $aref) = @_;
    open(FILE, $file)||die("Cannot open $file");
    while (<FILE>)
    {
	my $line = $_;
	chomp $line;
	my @sp = split(/\t/, $line);
	my $name = $sp[0];
	push @$aref, $name;
    }
    close(FILE);
}

sub get_fastq_filelists
{
    my ($href_se_lists, $href_pe_lists, 
	$indx, $outdir_bins)=@_;
    
    ## index format
    # name se_list pe_list1  pe_list2

    open(SAM, $indx)||die("Cannot openfastq  index file $indx");
    while (<SAM>)
    {
	my $line = $_;
	chomp $line;
	my @sp = split(/\t/, $line);
	my $name = $sp[0];
	my $se   = $sp[1];
	if (get_num_lines($se)==0)
	{
	    $href_se_lists->{$name} = "NO";	    
	}
	else
	{
	    $href_se_lists->{$name} = $se;
	}
	my $pe1   = $sp[2];
	my $l1 = get_num_lines($pe1);
	my $pe2   = $sp[3];
	my $l2 = get_num_lines($pe2);
	if ($l1 != $l2)
	{
	    die("These lists of paired end files are not the same length");
	}
	elsif ($l1>0)
	{
	    $href_pe_lists->{$name} = $pe1.",".$pe2;
	}
	else
	{
	    $href_pe_lists->{$name} ="NO";
	}
    }
    close(SAM);
}

sub get_num_lines
{
    my ($file) = @_;
    my $cmd = "wc -l $file";
    my $ret = qx{$cmd};
    if ($ret =~ /(\d+)\s+$file/)
    {
	return $1;
    }
    else
    {
	die("$cmd returns $ret - unexpected")
    }
}
sub build_all_unclean_binaries
{
    my ($index, $odir_bins,  
	$href_sam_to_uncleaned, $href_sam_to_uncleaned_log, $href_sam_to_covg,
	$cortex_directory, $qual, $dup, $hom, $hei, $widt, $max_read, $format) = @_;
    my %se_lists=();
    my %pe_lists=();
    get_fastq_filelists(\%se_lists, \%pe_lists, $index, $odir_bins);
    foreach my $sample (keys %se_lists)
    {
	foreach my $k (@kmers)
	{
	    build_unclean($sample, $se_lists{$sample}, $pe_lists{$sample}, $k, 
			  $outdir_binaries, $href_sam_to_uncleaned, $href_sam_to_uncleaned_log,
			  $href_sam_to_covg,
			  $cortex_directory, $qual, $dup, $hom, $hei, $widt, $max_read, $format);
	}
    }
}


sub build_unclean
{
    my ($name, $se, $pe, $km, $out, $href, $href_log, $href_covg, $cdir, $q, $dupremoval, $hp,
	$height, $width, $max_r, $format) = @_;

    if ($out !~ /\/$/)
    {
	$out = $out.'/';
    }
    my $c1 = "mkdir -p $out"."uncleaned/$km";
    if (!(-d $out."uncleaned/$km"))
    {
	qx{$c1};
    }
    my $ctx = $out."uncleaned/$km/".$name.".unclean.kmer".$km;
    if ($q>0)
    {
	$ctx = $ctx.".q$q";
    }

    if ($hp>0)
    {
	$ctx = $ctx."remv_hp".$hp;
    }
    $ctx = $ctx.".ctx";
    my $log = $ctx.".build_log";
    my $covg = $ctx.".covg";

    $href->{$name}->{$km}=$ctx;
    $href_log->{$name}->{$km}=$log;
    $href_covg->{$name}->{$km}=$covg;

    if (-e $ctx)
    {
	if (!(-e $log))
	{
	    print $tee "Binary $ctx exists, so will not rebuild, but the log file is missing. Carrying on nevertheless.\n";
	}
	if (!(-e $covg))
	{
	    print $tee "Binary $ctx exists, so will not rebuild, but the covg distribution file is missing\n";
	}
	print $tee "Binary $ctx already exists, so will not rebuild\n";
	return;
    }

    my $cortex_binary = get_right_binary($km, $cdir,1 );##one colour
    my $cmd = $cortex_binary." --kmer_size $km --mem_height $height --mem_width $width --dump_binary $ctx --max_read_len $max_r  --format $format --dump_covg_distribution $covg";
    if ($se ne "NO")
    {
	$cmd = $cmd." --se_list $se ";
    }
    if ($pe ne "NO")
    {
	$cmd = $cmd." --pe_list $pe ";
    }
    if ($q>0)
    {
	$cmd=$cmd." --qualty_score_threshold $q";
    }
    if ($hp>0)
    {
	$cmd = $cmd."--cut_homopolymers $hp";
    }
    if ($dupremoval)
    {
	$cmd = $cmd." --remove_pcr_duplicates ";
    }

    $cmd = $cmd." > $log 2>&1";
    my $ret = qx{$cmd};
    print $tee "$cmd\n$ret\n";
    if (!(-e $ctx))
    {
	die("Unable to build $ctx");
    }
}

sub get_right_binary
{
    my ($k, $root, $num_cols) = @_;
    
    if ($root !~ /\/$/)
    {
	$root = $root.'/';
    }
    my $bindir = $root."/bin/";
    my $maxk=(roundup($k/32) * 32 -1);

    my $cortex_binary = $bindir."cortex_var_".$maxk."_c".$num_cols;
    if (!(-e $cortex_binary))
    {
	die("Abort - please go and compile $cortex_binary\n");
    }
    else
    {
	return $cortex_binary;
    }
}




sub get_number_samples
{
    my ($index) = @_;
    open(FILE, $index)||die("Cannot open $index");
    my $num=0;
    while (<FILE>)
    {
	my $line = $_;
	chomp $line;
	my @sp = split(/\t/, $line);
	if (scalar (@sp) != 4)
	{
	    die("Format error in fastaq index $index - each line should be tab separated with 4 fields. Name of sample\tse_list\tpe_list1\tpe_list2\n");
	}
	else
	{
	    $num++;
	}
    }
    close(FILE);
    return $num;
}


sub run_checks
{
    if (-e $global_logfile)
    {
	die("Abort. run_calls.pl will either output its logs to a file called \"default_logfile\", or if you have specified --logfile, then it outputs to whatever you specify. However this file $global_logfile, already exists. Specify another, or delete this one\n");
    }
    if ( ($do_union ne "yes") && ($do_union ne "no") )
    {
	die("If you specify --do_union, you must give the value \"yes\" or \"no\"\n");
    }
    if ($vcftools_dir eq "/path/to/vcftools/dir")
    {
	die("You must specify the VCFTools directory, either on the commandline with --vcftools_dir, or by editing run_calls.pl manually (it's highlighted for you at the top of the file)\n");
    }
    if ($use_ref eq "no")
    {
	$number_of_colours = get_number_samples($fastaq_index);
    }
    else
    {
	$number_of_colours = get_number_samples($fastaq_index)+1;
    }
    if ($expt_type eq "")
    {
	die("You must specify --expt_type as one of EachColourADiploidSample, EachColourADiploidSampleExceptTheRefColour,  EachColourAHaploidSample,EachColourAHaploidSampleExceptTheRefColour\n");
    }
    if ($do_pd eq "yes")
    {
	if ( $list_ref_fasta eq "nonexistent_nonsense" )
	{
	    die("If you specify --pd yes, then you must also specify --list_ref_fasta\n");
	}
	elsif (!(-e $list_ref_fasta ))
	{
	    die("Cannot open the list of ref fasta you specified, $list_ref_fasta\n");
	}
    }
    if ($outdir !~ /\/$/)
    {
	$outdir = $outdir .'/';
    }
    if (!(-d $outdir))
    {
	my $cmd = "mkdir -p $outdir";
	qx{$cmd};
    }
    
    if ( ($use_ref eq "yes") && ($refbindir eq ""))
    {
	die("Since you said yes to --use_ref, You must specify --refbindir");
    }

    if ( ($do_auto_cleaning eq "yes") && ($genome_size==0) )
    {
	die("If you want automatic cleaning, you need tos pecify --genome_size");
    }
    if ($fastaq_index ne "")
    {
	if ($max_read_len==-1)
	{
	    die("You must specify max read length if you are passing in a fasta/q index");
	}
	if ( ($format ne "FASTA") && ($format ne "FASTQ") )
	{
	    die("You must specify --format (FASTA or FASTQ)  if you are passing an index");
	}
    }
    if ( ($fastaq_index eq "") && ($binary_colourlist eq "") )
    {
	die("Either specify a fastq index or a colourlist");
    }
    if ($mem_height==-1)
    {
	die("Must specify --mem_height");
    }
    if ($mem_width==-1)
    {
	die("Must specify --mem_width");
    }
    if ($first_kmer==0)
    {
	die("You must specify at least one kmer. To run for just one kmer, say 31, use --first_kmer 31\n");
    }
    if ($first_kmer % 2 ==0)
    {
	die("Only odd kmers allowed\n");
    }
    if ($last_kmer==0)
    {
	$last_kmer=$first_kmer;
    }
    if ($kmer_step <=0)
    {
	die("--kmer_step must have positive integer argument");
    }

    ## check ref binary exists
    opendir (DIR, $refbindir) or die("Cannot open reference binary directory $refbindir\n");
    my @files=();
    while (my $file = readdir(DIR)) {
	push @files, $file;
    }
    my $z;
    for ($z=$first_kmer; $z<=$last_kmer; $z+=$kmer_step)
    {
	# is there a binary with this k in its name in that directory
	my $found=0;
	foreach my $f (@files)
	{
	    if (($f =~ /$z/) && ($f =~ /.ctx/))
	    {
		$found=1;
		$k_to_refbin{$z}=$refbindir.$f;
	    }
	}
	if ($found==0)
	{
	    die("Cannot find a reference binary fr k=$z in the specified dir $refbindir\n");
	}
    }
    if ( ($do_auto_cleaning ne "yes") && ($do_auto_cleaning ne "no") )
    {
	print $tee ("--auto_cleaning takes an argument that must be \"yes\" or \"no\", but you have given it this value ");
	print $tee $do_auto_cleaning;
	die("\n");
    }
    if ( ($do_user_spec_cleaning ne "yes") && ($do_user_spec_cleaning ne "no") )
    {
	die("--user_cleaning takes an argument that must be \"yes\" or \"no\"");
    }

    if ($do_user_spec_cleaning eq "yes")
    {
	if ( ($user_min_clean==0) && ($user_max_clean==0)  )
	{
	    print $tee ("If you specify --user_clean that means you want to tell the script what cleaning threshold(s) to use.\n");
	    print $tee "If you want to use just one threshold, say 2, use --user_min_clean 2\n";
	    print $tee "If you want to do many then do --user_min_clean 2 --user_max_clean 10 --user_clean_step 2\n";
	    print $tee "This will do 2,4,6,8,10\n";
	    die();
	}
	else
	{
	    print $tee "ZAMZAM $user_min_clean $user_max_clean $user_clean_step\n";
	}
	if ( ($user_max_clean-$user_min_clean) % $user_clean_step !=0)
	{
	    die("The user min and max cleaning thresholds must differ by a multiple of user_clean_step\n");
	}

    }
    if ( ($do_auto_cleaning eq "yes") && ($do_user_spec_cleaning eq "yes") )
    {
	die("You cannot specify both automatic cleaning and user-specified cleaning\n");
    }
    if ($do_user_spec_cleaning eq "yes")
    {
	if ($user_min_clean > $user_max_clean)
	{
	    die("you have specified --user_min_clean above --user_max_clean\n");
	}
    }
    #if ($bc_col_args ne "")
    #{
#	$do_bc="yes";
	#my @sp = split(/\//, $bc_col_args);
	#if (scalar @sp ne 2)
	#{
	#    die("--bc argument should be comma-separated numbers, a slash /, and then more comma-separated numbers");
	#}
	#my @sp0_sp = split(/,/,$sp[0]);
	#my $i;
	#for ($i=0; $i<scalar(@sp0_sp); $i++)
	#{
	#    if ($sp0_sp[$i] !~ /[-]{0,1}\d+/)
	#    {
	#	die("--bc argument should be comma-separated numbers, a slash /, and then more comma-separated numbers");
	#    }
	#}


	#my @sp1_sp = split(/,/,$sp[1]);
	#for ($i=0; $i<scalar(@sp1_sp); $i++)
	#{
	#    if ($sp1_sp[$i] !~ /[-]{0,1}\d+/)
	#    {
	#	die("--bc argument should be comma-separated numbers, a slash /, and then more comma-separated numbers");
	#    }
	#}
    #}

    if ( ($do_bc ne "yes") && ($do_bc ne "no") )
    {
	die("--bc must be yes or no");
    }
    if ( ($do_pd ne "yes") && ($do_pd ne "no") )
    {
	die("--pd must be yes or no");
    }
    if (!(-d $outdir))
    {
	my $cmd1 = "mkdir -p $outdir";
	qx{$cmd1};
	my $cmd2 = "mkdir -p $outdir_binaries";
	qx{$cmd2};
	my $cmd3 = "mkdir -p $outdir_calls";
	qx{$cmd3};
	my $cmd4 = "mkdir -p $outdir_vcfs";
	qx{$cmd4};
	
    }
    if ($fastaq_index ne "")
    {
	if (! (-e $fastaq_index))
	{
	    die("Cannot access this fastaq index: $fastaq_index");
	}
    }


    if ($qthresh!=-1)
    {
	if ($qthresh !~ /^\d+/)
	{
	    die("--qthresh must be numeric");
	}
	if ($qthresh<=0)
	{
	    die("--qthresh expects a positive integer");
	}
    }
    if ($homopol!=-1)
    {
	if ($homopol<=0)
	{
	    die("--homopol expects a positive integer.");
	}
    }
}


sub roundup {
    my $n = shift;
    return(($n == int($n)) ? $n : int($n + 1))
}



sub get_expected_depth_and_cvg_distrib
{
    my ($file, $kmer, $genome_len) = @_;

    open(LOG, $file)||die();
    while (<LOG>)
    {
	my $logline = $_;
	chomp $logline;
	if ($logline =~ /SUMMARY/)
	{
	    <LOG>;
	    $logline = <LOG>;
	    if ($logline =~ /\d+\s+(\d+)\s+(\d+)/)
	    {
		my $mean_read_len = $1;
		my $bases = $2;
		my $depth = $bases/$genome_len;
		my $deff= $depth* ($mean_read_len-$kmer+1)/$mean_read_len;
		close(LOG);
		return $deff;
	    }
	    else
	    {
		die("Bad regexp on line $logline");
	    }
	}
    }
    close(LOG);
    die("Bad parsing");
}



sub get_cleaning_thresh_and_distrib
{
    my ($file, $exp_covg, $aref) = @_;
    open(CLEANINGFILE, $file)||die();
    my @covgs=();
    my $count=1;
    my $min_index=1;
    my $min_val=999999999999999;
    <CLEANINGFILE>;#ignore number with zero covg
    while (<CLEANINGFILE>)
    {
	my $line = $_;
	chomp $line;

	#multiplicity:1  Number:8
	if ($line =~ /multiplicity:(\d+)\s+Number:(\d+)/)
	{
	    my $mult = $1;
	    my $num = $2;
	    if ($mult != $count)
	    {
		die("counting error mult is $mult and count is $count");
	    }
	    push @covgs, $num;
	    push @$aref, $num;

	}
	$count++;
    }
    close(CLEANINGFILE);

#    return $min_index;
    my $thresh=$exp_covg;
    my $i;
    for ($i=int($exp_covg); $i>0; $i--)
    {
	if ($covgs[$i]<$min_val)
	{
	    $min_val=$covgs[$i];
	    $min_index=$i;
	}
    }
    return $min_index;
}


sub get_max_read_len_of_fasta
{
    my ($fasta) = @_;
    open(FA, $fasta)||die("Cannot open $fasta");
    my $max_read_len=0;
    while (<FA>)
    {
	my $line = $_;
	chomp $line;
	if ($line =~ /^>/)
	{
	    $line = <FA>;
	    chomp $line;
	    my $len = length($line);
	    if ($len>$max_read_len)
	    {
		$max_read_len = $len;
	    }
	}
    }
    close(FA);
    return $max_read_len;
}


sub get_max
{
    my ($a, $b)  = @_;
    if ($a>=$b)
    {
	return $a;
    }
    else
    {
	return $b;
    }
}
sub get_min
{
    my ($a, $b)  = @_;
    if ($a<=$b)
    {
	return $a;
    }
    else
    {
	return $b;
    }
}





sub print_build_stats
{
    my ($href_sample_to_uncleaned_log, $file) = @_;
    
}


sub get_ref_fasta
{
    my ($list) = @_;
    my $out = $outdir_vcfs."temp_ref.fa";
    open(LI, $list)||die();
    while (<LI>)
    {
	my $lyn = $_;
	chomp $lyn;
	my $cmd ="cat $lyn >> $out";
	qx{$cmd};
    }
    close(LI);
    return $out;
}
