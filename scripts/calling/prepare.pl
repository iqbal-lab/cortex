#!/usr/bin/perl -w
use strict;

use File::Basename;
use File::Spec;
use Getopt::Long;
use Benchmark;
use Cwd    qw( abs_path );
use FindBin qw($Bin);

my $cortex_dir =~ s/scripts\/calling\/prepare.pl//;

my $vcftools_dir = "";
my $stampy_bin="";
my $stampy_hash="";

&GetOptions(
    'index:s'                     =>\$all_samples_index,
    'ref_fa:s'                    =>\$ref_fa,
    'dir_for_ref_objects:s'       =>\$refdir, 
    'vcftools_dir:s'              =>\$vcftools_dir,
    'outdir:s' §                  =>\$outdir,
    'stampy_bin:s'                =>\$stampy_bin,
    'stampy_hash:s'               =>\$stampy_hash,
    'kmer:i'                      =>\$kmer,
    'genome_size:i'                 =>\$genome_size,
    );


$outdir = add_slash($outdir);



### Build reference binary
my $cmd_list =  

