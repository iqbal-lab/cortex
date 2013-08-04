#!/usr/bin/perl -w
use strict;

my $callingscript_dir;
my $analyse_variants_dir;
my $cortex_dir;
my $isaac_bioinf_dir;

BEGIN
{
	use FindBin;
	$callingscript_dir = $FindBin::Bin;
	$callingscript_dir = $callingscript_dir.'/';
	$cortex_dir = $callingscript_dir;
	$cortex_dir =~ s/scripts\/calling//;
	$analyse_variants_dir = $cortex_dir."/scripts/analyse_variants/";
	$isaac_bioinf_dir = $analyse_variants_dir."bioinf-perl/";
	
	push( @INC,
		$cortex_dir
		  . "/scripts/analyse_variants/perl_modules/Statistics-Descriptive-2.6",
	      $isaac_bioinf_dir."lib/",
	      $callingscript_dir
	    );
}

my $check_perl5 = "echo \$PERL5LIB";
my $check_perl5_ret = qx{$check_perl5};
my $isaac_libdir = $isaac_bioinf_dir."lib";
if ($check_perl5_ret !~ /$isaac_libdir/) 
{
    $ENV{PERL5LIB} .= ":$isaac_libdir";
}
if ($check_perl5_ret !~ /$callingscript_dir/)
{
    $ENV{PERL5LIB} .= ":$callingscript_dir";
}
use RunCallsLib;

my $file = "union_calls.vcf.uncleaned";
my @arr=($file);

my $outfilename_joint = "tmp/union_calls.vcf.output.joint";
my $outfilename_indep = "tmp/union_calls.vcf.output.indep";
my $vcftools_dir="/home/zam/installed_apps/vcftools_0.1.9/";


combine_and_filter_dups_overlaps_and_ref_mismatches_from_vcf($outfilename_indep, "/home/zam/dev/hg/bitbucket/CORTEX_mainline/scripts/calling/tmp/", \@arr, $vcftools_dir, "independent",2);


combine_and_filter_dups_overlaps_and_ref_mismatches_from_vcf($outfilename_joint, "/home/zam/dev/hg/bitbucket/CORTEX_mainline/scripts/calling/tmp/", \@arr, $vcftools_dir, "joint",2);
