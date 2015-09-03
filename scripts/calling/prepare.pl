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
    $cortex_dir =~ s/scripts\/calling\/prepare.pl//;
    push ( @INC, $cortex_dir."scripts/calling/");
}

use BasicUtils qw ( check_cortex_runnable add_slash is_fastq count_bases_in_fasta create_dir_if_does_not_exist );



## Script to prepare reference graph binary, Stampy hash of reference etc

my $vcftools_dir = "";
my $stampy_bin="";
my $stampy_hash="";
my $all_samples_index="";
my $ref_fa="";
my $refdir = "";
my $outdir="";
my $kmer="";
my $genome_size="";
my $ref_id = "REF";
my $verbose ='';

&GetOptions(
    'index:s'                     =>\$all_samples_index,
    'ref_fa:s'                    =>\$ref_fa,
    'dir_for_ref_objects:s'       =>\$refdir, 
    'vcftools_dir:s'              =>\$vcftools_dir,
    'outdir:s'                    =>\$outdir,
    'stampy_bin:s'               =>\$stampy_bin,
    #'stampy_hash:s'               =>\$stampy_hash,
    'ref_id:s'                    =>\$ref_id, ##not mandatory
    'kmer:i'                      =>\$kmer,
    'verbose'                    => \$verbose,
    );





## will make outdir if it does not exist, as well as refdir
## it returns absolute paths 
my $num_samples=0;
my $mem_height=0;
my $mem_width=0;
($num_samples, $mem_height, $mem_width, $genome_size,
 $ref_fa, $stampy_bin, $stampy_hash, $outdir )  
    
    ## includes creating missing dirs and checking Stampy, VCFtools
    = check_args($all_samples_index, $ref_fa, $refdir, $vcftools_dir, $outdir, 
		 $stampy_bin, $stampy_hash, $kmer, $ref_id);

$outdir = BasicUtils::add_slash($outdir);



##Variables which I will write to the config file
my %arr_config = ("ref_fa" => abs_path($ref_fa), 
		  "refbindir"=> abs_path($refdir."/ctx_bins"), 
		  "vcftools_dir"=>abs_path($vcftools_dir),
		  "stampy_bin"=>abs_path($stampy_bin), 
		  #"stampy_hash"=>$stampy_hash, 
		  "genome_size"=>$genome_size);


my $cmd = "export LD_LIBRARY_PATH=$cortex_dir"."/libs/gsl-1.15/:$cortex_dir"."/libs/gsl-1.15/blas/:$cortex_dir"."/libs/gsl-1.15/.libs/:$cortex_dir"."/libs/gsl-1.15/blas/.libs/:$cortex_dir"."/libs/gsl-1.15/cblas/.libs/:\$"."LD_LIBRARY_PATH";
qx{$cmd};


### Compile Cortex for 1,2,num_samples, num_samples+1 colours, and return the path to the 1-colour binary
my $cortex_bin = compile_cortex($cortex_dir, $num_samples, $kmer);

### Check can run Cortex
check_cortex_runnable($cortex_bin);


### Build reference binary
my $ref_falist =$outdir."filelists/"."ref_list";
my $cmd_list =  "ls $ref_fa > $ref_falist";
qx{$cmd_list};
my $refbin = $refdir."ctx_bins/".$ref_id.".k".$kmer.".ctx";
my $ref_bin_log = $refbin.".log";
my $cmd_build = "$cortex_bin --se_list $ref_falist --kmer $kmer ";
    $cmd_build .= "--mem_height $mem_height --mem_width $mem_width ";
    $cmd_build .= "--sample_id $ref_id --dump_binary $refbin > $ref_bin_log";
    $cmd_build .= " 2>&1";
if (!(-e $refbin))
{
    qx{$cmd_build};
    if (! -e($refbin))
    {
        die("Failed to build $refbin - check error in logfile $ref_bin_log\n")
    }
    else
    {
	if ($verbose)
	{
	    print "Created reference genome Cortex binary graph file $refbin\n";
	}
    }
}
else
{
    if ($verbose)
    {
	print "reference binary $refbin already exists\n";
    }
}

$arr_config{"list_ref"}=abs_path($ref_falist);
$arr_config{"ref_fa"} =abs_path($ref_fa);
$arr_config{"refbin"}= abs_path($refbin);



### Build stampy hash
my $stampy_stub = $refdir."stampy/".$ref_id;
if ( !( (-e $stampy_stub.".stidx") && (-e $stampy_stub.".sthash") ) )
{
    ##if they don't both exist, kill the one that does
    if (-e $stampy_stub.".stidx")
    {
	my $rmv1 = "rm $stampy_stub".".stidx";
	qx{$rmv1};
    }
    elsif (-e $stampy_stub.".sthash")
    {
        my $rmv1 = "rm $stampy_stub".".sthash";
        qx{$rmv1};
    }
}
if (!(-e $stampy_stub.".stidx"))
{
    my $cmd_stampy1 = "$stampy_bin -G $stampy_stub $ref_fa 2>&1";
    my $ret_stampy1 = qx{$cmd_stampy1};
    my $cmd_stampy2 = "$stampy_bin -g $stampy_stub -H $stampy_stub 2>&1";
    my $ret_stampy2 = qx{$cmd_stampy2};
    if ($verbose)
    {
	print "Built Stampy hash\n";
	#print "$cmd_stampy1\n$ret_stampy1\n$cmd_stampy2\n$ret_stampy2\n";
    }
}
$arr_config{"stampy_hash"}= abs_path($stampy_hash);

if ($outdir ne "no")
{
    ##  Will save useful info in a config file
    my $config = $outdir."config.prep.txt";
    my $c_fh;#c_fh means Config File Handle
    open ($c_fh, ">".$config)||die("Cannot create $config file - file permissions 
    issue?\n");
    
    add_to_config($c_fh, \%arr_config);
    close($c_fh);
}

###################################################################

sub check_args

{
    my ($loc_index, $loc_ref_fa, $loc_refdir, $loc_vcftools_dir, 
        $loc_outdir, $loc_stampy_bin, $loc_stampy_hash, 
        $loc_kmer, $loc_refid) = @_;
    
    if ($loc_index eq "no")
    {
    }
    elsif (($loc_index eq "") || (!(-e $loc_index)))
    {
	my $str = "You must generally specify an index file with --index - the file format is as follows.\n";
	$str .= "Each sample corresponds to one line of this file.\n";
	$str .= "Each line is tab separated with 4 fields.";
	$str .= "Field1 = sample identifier\n";
	$str .= "Field2 = full path to a filelLIST, listing all fastq or bams";
	$str .= "         for this sample. Do NOT put the name of an actual";
	$str .= "         fastq/bam in this field - it expects a LIST.\n";
	$str .= "         These sequence will be treated as single-ended, but\n";
	$str .= "         Cortex only uses paired-end info to remove PCR \n";
	$str .= "         duplicates.\n";
	$str .= "Field3 = As field 2, for left-hand reads from paired ends\n";
	$str .= "         You can just put a dot \".\" to ignore this.\n";
	$str .= "Field3 = As field 3, for right-hand reads from paired ends\n";
	$str .= "         You can just put a dot \".\" to ignore this.\n";
	$str .= "If you actively want to set up reference utils for a species but not\n";
	$str .= "for some specific dataset, use \"--index no\" .\n";
	#die("Please create this file and specify it using --index\n");
    }
    

    ## R must be in path
    my $rcom = "R --version";
    my $rret = qx{$rcom};
    if ($rret !~ /R version/)
    {
	die("R must be in your path\n");
    }

    my $num_samples=0;
    if ($loc_index ne "no")
    {
	my $num_samples_cmd = "wc -l $loc_index";
	$num_samples = qx{$num_samples_cmd};
	chomp $num_samples;
	if ($num_samples =~ /^(\d+)/)
	{
	    $num_samples = $1;
	}
	else
	{
	    $num_samples=0;
	}
    }
    
    if (!(-e $loc_ref_fa))
    {
	die("Unable to open your specified reference fasta $loc_ref_fa\n")
    }
    else 
    {
	my $err = BasicUtils::is_fasta($loc_ref_fa);
	if ($err ne "EFasta")
	{
	    my $str = "Your specified reference fasta file $loc_ref_fa appears ";
	    $str .= "not to be in fasta format\n";
	    die($str."   ".$err);
	}
	$loc_ref_fa = abs_path($loc_ref_fa);
    }
    
    
    
### Calculate genome size and memory parameters
    my $genome_len = BasicUtils::count_bases_in_fasta($loc_ref_fa);
    
    my $mem_w=100;
## 2^height * width = 3g
## => height = log_2(3g/width)
    my $mem_h = int(log(3*$genome_len/100)/log(2) +0.5);
    
    
## Set up the reference directory
    if (!(-d $loc_refdir))
    {
	my $c = "mkdir -p $loc_refdir";
	qx{$c};
	if (!(-d $loc_refdir))
	{
	    my $str = "You specified $loc_refdir as directory for ref objects.\n";
	    $str .= "This did not exist, so we tried to create it, and that failed\n";
	    die($str);
	}
    }
    $loc_refdir = abs_path($loc_refdir);
    $loc_refdir = BasicUtils::add_slash($loc_refdir);
    BasicUtils::create_dir_if_does_not_exist($loc_refdir, "check_args of prepare.pl");
    my $st_dir = $loc_refdir."stampy/";
    my $ctxdir = $loc_refdir."ctx_bins/";
    BasicUtils::create_dir_if_does_not_exist($st_dir,  "check_args of prepare.pl");
    BasicUtils::create_dir_if_does_not_exist($ctxdir, "check_args of prepare.pl");
    

### is the VCFtools directory actualy a directory with the right stuff?
    if (!(-d $loc_vcftools_dir))
    {
	die("The vcftools directory specified ($loc_vcftools_dir) doesn't exist\n");
    }
    else
    {
	$loc_vcftools_dir = BasicUtils::add_slash($loc_vcftools_dir);

	my $d1 = $loc_vcftools_dir."perl";
	if (!(-e $d1))
	{
	    my $str = "The specified vcftools_dir is not actually the directory ";
	    $str .= "of a full VCFtools install/download. At the least it should ";
	    $str .= " contain the perl/ directory that there when you download a ";
	    $str .= " .tgz file of VCFtools\n";
	    die($str);
	}
    }
    
    
## Sort out the output directory
    if ($outdir ne "no")
    {
	BasicUtils::create_dir_if_does_not_exist($loc_outdir, "check_args of prepare.pl");
	$loc_outdir = abs_path($loc_outdir);
	$loc_outdir = BasicUtils::add_slash($loc_outdir);
	my $filelist_dir = $loc_outdir."filelists/";
	BasicUtils::create_dir_if_does_not_exist($filelist_dir,  "check_args of prepare.pl");
    

## Create the reference fasta filelist
	my $ref_fa_list = $filelist_dir."list_ref_fa";
	my $cmd_create_reffalist = "ls $loc_ref_fa > $ref_fa_list";
	qx{$cmd_create_reffalist};
    }
    
## Check the stampy binary
    if (!(-e $loc_stampy_bin))
    {
	die("You must specify the full path to stampy.py using --stampy_bin\n");
    }
    elsif ($loc_stampy_bin !~ /stampy.py$/)
    {
	die("--stampy_bin should give full path to stampy.py\n")
    }
    
    my $stcmd = "$stampy_bin --help";
    my $stret = qx{$stcmd};
    
    if ($stret !~ /Usage/)
    {
        my $str = "The stampy.py specified does not seem to run (tried "; 
        $str .= " calling --help. Maybe Python 2.6 or 2.7 not present?\n";
        die($str);
    }
    $loc_stampy_bin = abs_path($loc_stampy_bin);
    if ($loc_stampy_hash eq "")
    {
	$loc_stampy_hash=$loc_refdir."stampy/REF";
    }
    $loc_stampy_hash = abs_path($loc_stampy_hash);

    if ($loc_kmer % 2==0)
    {
	die("Kmer must be an odd number\n");
    }	
    
    return ($num_samples, $mem_h, $mem_w, $genome_len,
	    $loc_ref_fa, $loc_stampy_bin, $loc_stampy_hash, $loc_outdir )    
}

sub compile_cortex
{
    my ($ctx_dir, $n_samples, $k) = @_;

    $ctx_dir = BasicUtils::add_slash($ctx_dir);
    my $maxk = 32 * (int($k / 32) +1) -1;

    my @nums = (1);
    if ($n_samples>0)
    {
	push @nums,2;
	push @nums, $n_samples;
	push @nums, $n_samples+1;
    }
    foreach my $n (@nums)
    {
	my $c = "export LD_LIBRARY_PATH=$ctx_dir"."/libs/gsl-1.15/:$ctx_dir"."/libs/gsl-1.15/blas/:$ctx_dir"."/libs/gsl-1.15/.libs/:$ctx_dir"."/libs/gsl-1.15/blas/.libs/:$ctx_dir"."/libs/gsl-1.15/cblas/.libs/:\$"."LD_LIBRARY_PATH; cd $ctx_dir; make NUM_COLS=$n MAXK=$maxk cortex_var 2>&1;";
	my $r = qx{$c};
	my $bin = $ctx_dir."bin/cortex_var_".$maxk."_c".$n;
	if (!(-e $bin))
	{
	    die("Failed to compile Cortex, tried to make $bin and failed\n");
	}

    }
    
    return  $ctx_dir."bin/cortex_var_".$maxk."_c1";
}


sub add_to_config
{
    my ($fh, $href) = @_;

    foreach my $key (keys %$href)
    {
	print $fh $key."\t".$href->{$key};
	print $fh "\n";
    }
}
