#!/usr/bin/perl -w
use strict;

use File::Basename;
use File::Spec;
use Getopt::Long;
use Benchmark;
use Cwd    qw( abs_path );
use FindBin qw($Bin);
use BasicUtils qw ( add_slash, is_fastq, count_bases_in_fasta );

my $cortex_dir = abs_path($0);
$cortex_dir =~ s/scripts\/calling\/prepare.pl//;

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

&GetOptions(
    'index:s'                     =>\$all_samples_index,
    'ref_fa:s'                    =>\$ref_fa,
    'dir_for_ref_objects:s'       =>\$refdir, 
    'vcftools_dir:s'              =>\$vcftools_dir,
    'outdir:s'                    =>\$outdir,
    'stampy_bin :s'               =>\$stampy_bin,
    'stampy_hash:s'               =>\$stampy_hash,
    'ref_id:s'                    =>\$ref_id,
    'kmer:i'                      =>\$kmer,
    );


$outdir = add_slash($outdir);



## will make outdir if it does not exist, as well as refdir
## it returns absolute paths 
my $num_samples, $mem_height, $mem_width;
    ($num_samples, $mem_height, $mem_width, $genome_size,
    $ref_fa, $stampy_bin, $stampy_hash, $outdir )  

        ## includes creating missing dirs and checking Stampy, VCFtools
        = check_args($index, $ref_fa, $refdir, $vcftools_dir, $outdir, 
	   	       		 $stampy_bin, $stampy_hash, $kmer, $genome_size);



##  Will save useful info ain a config file
my $config = $outdir."config.prep.txt";
open ($c_fh, ">".$config)||die("Cannot create $config file - file permissions 
    issue?\n");

##Variables which I will write to the config file
my @arr_config = ($ref_fa, $refdir, $vcftools_dir,
                $stampy_bin, $stampy_hash, $genome_size);
add_to_config($c_fh, \@arr_config);


### Compile Cortex
my $cortex_bin = compile_cortex($cortex_dir, $num_samples);
### Check can run Cortex
check_cortex_runnable($cortex_bin);

### Build reference binary
my $ref_falist =$outdir."filelists/"."ref_list";
my $cmd_list =  "ls $ref_fa > $ref_falist";
qx{$cmd_list};
my $refbin = $refdir."ctx_bins/".$ref_id.".k".$kmer.".ctx";
my $ref_bin_log = $refbin.".log";
my $cmd_build = "$cortex_bin --se_list $ref_falist --kmer $kmer ";
    $cmd_build .= "--mem_height $mem_height --mem_width mem_width ";
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
        print "Created reference genome Cortex binary graph file $refbin\n";
    }
}
else
{
    print "reference binary $refbin already exists\n";
}
@arr_config = ($ref_falist, $ref_fa, $refbin);
add_to_config($c_fh, \@arr_config);


### Build stampy hash
my $stampy_stub = $refdir."stampy/".$ref_id;
my $cmd_stampy1 = "$stampy_bin -G $stampy_stub $ref_fa";
my $ret_stampy1 = qx{$cmd_stampy1};
my $cmd_stampy2 = "$stampy_bin -g $stampy_stub -H $stampy_stub";
my $ret_stampy2 = qx{$cmd_stampy2};
print "$cmd_stampy1\n$ret_stampy1\n$cmd_stampy2\n$ret_stampy2\n";



###################################################################

sub check_args

{
    my ($loc_index, $loc_ref_fa, $loc_refdir, $loc_vcftools_dir, 
        $loc_outdir, $loc_stampy_bin, $loc_stampy_hash, 
        $loc_kmer) = @_;

if (($loc_index eq "") || (!(-e $loc_index)))
{
    my $str = "You must specify an index file with --index.\n";
    $str .= "Each sample corresponds to one line of this file.\n";
    $str .= "Each line is tab separated with 4 fields."
    $str .= "Field1 = sample identifier\n";
    $str .= "Field2 = full path to a filelLIST, listing all fastq or bams";
    $str .= "         for this sample. Do NOT put the name of an actual";
    $str .= "         fastq/bam in this field - it expects a LIST.\n";
    $str .= "         These sequence will be treated as single-ended, but\n";
    $str .= "         Cortex only uses paired-end info to remove PCR \n";
    $str .= "         duplicates.\n"
    $str .= "Field3 = As field 2, for left-hand reads from paired ends\n";
    $str .= "         You can just put a dot \".\" to ignore this.\n";
    $str .= "Field3 = As field 3, for right-hand reads from paired ends\n";
    $str .= "         You can just put a dot \".\" to ignore this.\n";
    die("Please create this file and specify it using --index\n");
}
my $num_samples_cmd = "wc -l $loc_index";
my $num_samples = qx{$num_samples_cmd};
chomp $num_samples;
if ($num_samples =~ /^(\d+)/)
{
    $num_samples = $1;
}
else
{
    die("Failed to count lines in index $loc_index - permissions issue?\n");
}

if (!(-e $loc_ref_fa))
{
    die("Unable to open your specified reference fasta $loc_ref_fa\n")
}
elsif (is_fasta($loc_ref_fa) ne "EFasta")
{
    my $str = "Your specified reference fasta file $loc_ref_fa appears ";
    $str .= "not to be in fasta format\n";
    die($str);
}



### Calculate genome size and memory parameters
my $genome_len = count_bases_in_fasta($loc_ref_fa);

my $mem_w=100;
## 2^height * width = 3g
## => height = log_2(3g/width)
my $mem_h = int(log(1.5*$genome_len/100)/log(2) +0.5);


## Set up the reference directory
create_dir_if_does_not_exist($loc_refdir, "check_args of prepare.pl");
$loc_refir = add_slash($loc_refdir);
my $st_dir = $loc_refdir."stampy/";
my $ctxdir = $loc_refdir."ctx_bins/";
create_dir_if_does_not_exist($st_dir,  "check_args of prepare.pl");
create_dir_if_does_not_exist($ctxdir, "check_args of prepare.pl");


### is the VCFtools directory actualy a directory with the right stuff?
if (!(-d $loc_vcftools_dir))
{
    die("The vcftools directory specified ($loc_vcftools_dir) doesn't exist\n");
}
else
{
    $loc_vcftools_dir=add_slash($loc_vcftools_dir);
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
create_dir_if_does_not_exist($loc_outdir, "check_args of prepare.pl");
$loc_outdir = add_slash($loc_outdir);
my $filelist_dir = $loc_outdir."filelists/";
create_dir_if_does_not_exist($filelist_dir,  "check_args of prepare.pl");

## Create the reference fasta filelist
my $ref_fa_list = $filelist_dir."list_ref_fa";
my $cmd_create_reffalist = "ls $loc_ref_fa > $ref_fa_list";
qx{$cmd_create_reffalist};


## Check the stampy binary
if (!(-e $loc_stampy_bin))
{
    die("You must specify the full path to stampy.py using --stampy_bin\n");
}
elsif ($loc_stampy_bin !~ /stampy.py$/)
{
    die("--stampy_bin shoould give full path to stampy.py\n")
}
my $stcmd = "$stampy_bin --help";
my $stret = qx{$stcmd};

if ($stret !~ /Usage/)
{
        my $str = "The stampy.py specified does not seem to run (tried "; 
        $str .= " calling --help. Maybe Python 2.6 or 2.7 not present?\n";
        die($str);
}

my $stf1 = $stampy_hash.".stidx";
my $stf2 = $stampy_hash.".sthash";
if ( (!(-e $stf1)) || (!(-e $stf2)) )
{
    ##They do not both exist

    if (-e $stf1) 
    {
        my $c1 ="rm $stf1";
        qx{$c1};
    }
    if (!( -e $stf2))
    {
        my $c2 ="rm $stf2";
        qx{$c2};
    }
}

if ($loc_kmer % 2==0)
{
    die("Kmer must be an odd number\n");
}

