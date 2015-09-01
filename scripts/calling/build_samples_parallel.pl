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
    $cortex_dir =~ s/scripts\/calling\/build_samples_parallel.pl//;
    push ( @INC, $cortex_dir."scripts/calling/");
}

use ConfigMgmt qw( get_from_config_if_undef print_to_config get_all_info_from_config_if_undef);
use BasicUtils qw ( add_slash );


my %vars = ( "index" => "",
	     "kmer" => "",
	     "qthresh" => 10,
	     "bc" => "yes",
	     "pd" => "no",
	     "mem_height"=>"",
	     "mem_width"=>"",
	     "vcftools_dir"=>"",
	     "stampy_bin"=>"",
	     "stampy_hash"=>"",
	     "list_ref"=>"",
	     "refbindir"=>"",
	     "ref_fa"=>"",
	     "num"=>-1,
	     "genome_size"=>"",
	     "outdir" => "",
	     "ref_binary" => "",
	     "index_dir"=> "");

#my $all_samples_index = "";
#my $kmer = 31;
#my $qthresh = 10;
#my $bc="yes";
#my $pd = "no";
#my $mem_height = 0;
#my $mem_width=0;
#my $vcftools_dir = "";
#my $stampy_bin="";
#my $stampy_hash="";
#my $list_ref="";
#my $refbindir="";
#my $ref_fa="";
#my $num=-1;
#my $genome_size=3000000;


#my $outdir="";

&GetOptions(
    'num:i'                       =>\$vars{"num"},,
    'index:s'                     =>\$vars{"index"},
    'list_ref:s'                  =>\$vars{"list_ref"},
    'refbindir:s'                 =>\$vars{"refbindir"},
    'vcftools_dir:s'              =>\$vars{"vcftools_dir"},
    'outdir:s'                   =>\$vars{"outdir"},
    'stampy_bin:s'                =>\$vars{"stampy_bin"},
    'stampy_hash:s'               =>\$vars{"stampy_hash"},
    'bc:s'                        =>\$vars{"bc"},
    'pd:s'                        =>\$vars{"pd"},
    'kmer:i'                      =>\$vars{"kmer"},
    'mem_height:i'                =>\$vars{"mem_height"},
    'mem_width:i'                 =>\$vars{"mem_width"},
    'genome_size:i'               =>\$vars{"genome_size"},
    'qthresh:i'                   =>\$vars{"qthresh"},
    );


check_outdir(\%vars);
$vars{"index_dir"} =$vars{"outdir"}."indexes/";
my $prev_config = $vars{"outdir"}."config.prep.txt";
my $this_config = $vars{"outdir"}."config.par.txt";


##anything specified on commandline here gets top priority.
##If not set by that, then config file from running these samples gets priority
#get_all_info_from_config_if_undef(\%vars, $this_config);

## anything not specified in on command line or in par config, we get from prep.config
get_all_info_from_config_if_undef(\%vars, $prev_config);

check_args(\%vars);

##remember many par.pl processes run in parallel, can't all write to same file
if ($vars{"num"}==1)
{
    my $fh_config;
    open($fh_config, ">".$this_config)||die("Unable to create $this_config\n");
    print_to_config("vcftools_dir", \%vars, $fh_config);
    print_to_config("kmer", \%vars, $fh_config);
    print_to_config("ref_binary", \%vars, $fh_config);
    print_to_config("genome_size", \%vars, $fh_config);
    print_to_config("mem_height", \%vars, $fh_config);
    print_to_config("mem_width", \%vars, $fh_config);
    close($fh_config);
}

my $index = $vars{"index_dir"}."index_".$vars{"num"};
my $c1 = "head -n ".$vars{"num"}." ".$vars{"index"}." | tail -n 1 > $index";
qx{$c1};
open(F, $index)||die("Unable to create an index file for job ".$vars{"num"});
my $line= <F>;
chomp $line;
my @sp = split(/\t/, $line);
my $sample = $sp[0];
my $odir = $vars{"outdir"}.$sample.'/';
my $c2 = "mkdir -p $odir";
qx{$c2};
my $log = $odir."log.".$sample;

my $cmd ="perl $cortex_dir"."scripts/calling/run_calls.pl --fastaq_index ".$vars{"index"};
$cmd .= " --first_kmer ".$vars{"kmer"};
$cmd .= " --auto_cleaning yes --bc ".$vars{"bc"};
$cmd .= " --pd ".$vars{"pd"};
$cmd .= " --outdir $odir --ploidy 2 --genome_size ".$vars{"genome_size"};
$cmd .= " --mem_height ".$vars{"mem_height"};
$cmd .= " --mem_width ".$vars{"mem_width"};
$cmd .= " --qthresh ".$vars{"qthresh"};
$cmd .= " --vcftools_dir ".$vars{"vcftools_dir"};
$cmd .= " --do_union yes --logfile $log,f --workflow independent --ref CoordinatesAndInCalling";
$cmd .= " --list_ref ".$vars{"list_ref"};
$cmd .= " --refbindir ".$vars{"refbindir"};
$cmd .= " --stampy_bin ".$vars{"stampy_bin"};
$cmd .= " --stampy_hash ".$vars{"stampy_hash"};
$cmd .= " --outvcf $sample ";

my $ret = qx{$cmd};


sub check_outdir
{
    my ($hashref) = @_;
    if ($hashref->{"outdir"} eq "")
    {
	my $errstr = "You must specify with --outdir the output directory, which \n";
	$errstr .= "must be the same directory used as output dir in prepare.pl\n";
	die($errstr);
    }
    $hashref->{"outdir"} = BasicUtils::add_slash($hashref->{"outdir"});
    if (!(-d $hashref->{"outdir"} ))
    {
	my $errstr = "Output directory does not exist:".$hashref->{"outdir"};
	$errstr .= "\n This should have been created previously by prepare.pl\n";
	die($errstr);
    }
}

sub get_args_from_config_if_undefined
{
    my ($vcft, $st_bin, $st_hash, 
	$list_r, $r_fa, $rbindir, $g_size, $conf) = @_;


    ## My current system does allow 
    ## you to have no config file, am ignorning
    ## the error/return values
    my $err;
    ($err, $vcft)   = get_from_config_if_undef($vcft, "vcftools_dir", $conf); 
    ($err, $st_bin) = get_from_config_if_undef($st_bin, "stampy_bin", $conf);
    ($err, $st_hash)= get_from_config_if_undef($st_hash, "stampy_hash", $conf);
    ($err, $list_r) = get_from_config_if_undef($list_r, "list_ref", $conf);
    ($err, $r_fa)   = get_from_config_if_undef( $r_fa, "ref_fa", $conf);
    ($err, $rbindir)= get_from_config_if_undef( $rbindir, "refbindir", $conf);
    ($err, $g_size) = get_from_config_if_undef( $g_size, "genome_size", $conf);

    return ($vcft, $st_bin, $st_hash, 
	    $list_r, $r_fa, $rbindir, $g_size);

}

sub check_args
{
    my ($hashref)=@_; #$n, $i_dir, $mh, $mw, $g, $km) = @_;

    if ($hashref->{"kmer"} eq "")
    {
	my $errstr = "Kmer size not specified. Normally I expect you to use prepare.pl,\n";
	$errstr .= " and that would enter kmer values into a config file, but that did not happen here.\n";
	$errstr .= " Failing that you should use --kmer in build_samples_parallel.pl \n";
	die($errstr);
    }
    else
    {
	if ($hashref->{"kmer"} % 2==0)
	{
	    die("Kmer must be an odd number\n");
	}
	elsif ($hashref->{"kmer"}<10) 
	{
	    die("I'm just going to stop this being run for k<10, can't believe it will be useful - typo?\n");
	}
    }
    if ($hashref->{"num"} == -1)
    {
	die("--num is a mandatory argument, specifies which sample from the INDEX file to run\n");
    }
    if (!(-d $hashref->{"index_dir"}))
    {
	my $c1 = "mkdir ".$hashref->{"index_dir"};
	qx{$c1};
    }	
    if ($hashref->{"genome_size"} eq "")
    {
	my $errstr = "Genome size is missing from config files - should have beens et by prepare,pl\n";
	$errstr .= "You can enter a genome size (in bp) using --genome_size.\n";
	$errstr .= "An estimate is fine. This will be used to get memory use parameters mem height/width\n";
	$errstr .=" (if you do not enter them) and later on will be used for likelihood calculations\n";
	die($errstr);
    }

    if ( ($hashref->{"mem_height"} eq "") 
	 ||  
	 ($hashref->{"mem_width"} eq "") )
    {
	$hashref->{"mem_width"} =100;

	## 2^height * width = 4g
	## => height = log_2(4g/width)
	$hashref->{"mem_height"} = int(log(4*$hashref->{"genome_size"}/100)/log(2) +0.5);
    }


    if ($hashref->{"vcftools_dir"} eq "")
    {
	my $errstr = "Either you should use prepare.pl, which would leave a config file\n";
	$errstr .= "specifying vcftools_dir (which it appears you have not),\n";
	$errstr .= "or you must specify --vcftools_dir.\n";
	die($errstr);
    }
    if ($hashref->{"stampy_bin"} eq "")
    {
	my $errstr = "Either you should use prepare.pl, which would leave a config file\n";
	$errstr .= "specifying stampy_bin - the path to stampy.py -  (which it appears you have not),\n";
	$errstr .= "or you must specify --stampy_bin.\n";
	die($errstr);
    }
    if ($hashref->{"stampy_hash"} eq "")
    {
	my $errstr = "Either you should use prepare.pl, which would leave a config file\n";
	$errstr .= "specifying stampy_hash (which it appears you have not),\n";
	$errstr .= "or you must specify --stampy_hash.\n";
	die($errstr);
    }
    if ($hashref->{"list_ref"} eq "")
    {
	my $errstr = "Either you should use prepare.pl, which would leave a config file\n";
	$errstr .= "specifying a filelist, listing reference genome fasta(s), (which it appears you have not),\n";
	$errstr .= "or you must specify --list_ref.\n";
	die($errstr);

    }
    if ($hashref->{"ref_fa"} eq "")
    {
	my $errstr = "Either you should use prepare.pl, which would leave a config file\n";
	$errstr .= "specifying a reference fasta file, (which it appears you have not),\n";
	$errstr .= "or you must specify --ref_fa.\n";
	die($errstr);
    }
    my $c = "cat ".$hashref->{"list_ref"};
    my $cr = qx{$c};
    chomp $cr;
    if (abs_path($cr) ne abs_path($hashref->{"ref_fa"}))
    {
	die("The reference fasta list does not correspond to (refer to) the ref fasta - inconsistent config file(s)\n");
    }

    if (!(-d $hashref->{"refbindir"}))
    {
	my $errstr ="The specified directory ".$hashref->{"refbindir"};
	$errstr .= " does not exist\n";
	die($errstr);
    }

    my @files = glob($hashref->{"refbindir"}."/*k".$hashref->{"kmer"}.".ctx");
    if (scalar @files==0)
    {
	my $errstr = "The reference binary directory ".$hashref->{"refbindir"};
	$errstr .= "does not contain any files with names ending k".$hashref->{"kmer"}.".ctx;\n";
	$errstr .= "either you have not followed the naming convention we are asking for, \n";
	$errstr .= "or there is no Cortex binary graph file of the reference genome in that directory\n";
	$errstr .= "If you use prepare.pl, all of this is done for you, so I suggest you do that\n";
	die($errstr);
    }
    elsif (scalar @files>1)
    {
	my $errstr = "There is more than one file in the reference binary directory ".$hashref->{"refbindir"};
	$errstr .= "\n with name ending k".$hashref->{"kmer"}.".ctx, so this script can't work out which is the reference binary\n";
	die($errstr);
    }
    $hashref->{"ref_binary"}=$files[0];
}
