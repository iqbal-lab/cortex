#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $script_dir;
my $isaac_dir;
my $lib_dir;
BEGIN
{
	use FindBin;
	$script_dir = $FindBin::Bin;
	$isaac_dir = $script_dir. "/bioinf-perl/fastn_scripts/";
	$lib_dir   = $script_dir. "/bioinf-perl/lib/";
	push(@INC, $lib_dir);
}

use lib "/Net/fs1/home/zam/dev/hg/CORTEX_release/scripts/analyse_variants/bioinf-perl/lib";
use FASTNFile;

my ($vcf, $sample, $out_fasta, $num_sites, $kmer, $chrom, $ref_fa, $help);
#default;
$num_sites=10000;
$help='';
&GetOptions(
    'vcf|v:s'      => \$vcf,
    'sample|s:s'   => \$sample,
    'out_fasta|f:s'=> \$out_fasta,
    'num_sites|n:i'=> \$num_sites,
    'kmer|k:s'     => \$kmer,
    'chrom|c:s'    => \$chrom,
    'ref_fa|r:s'   => \$ref_fa,
    'help'         => \$help,
    );

if ($help)
{
    print "Script to generate a fasta with pairs of reads. \n";
    print "First read is ref allele of a site known to be hom-alt in the sample of interest\n";
    print "Second read is the alt allele of that same site\n";
    print "\n\n";
    print "--vcf\t\t\t\tVCF file of SNP genotyping results\n";
    print "--sample\t\t\tName of sample of interest, must be in the VCF\n";
    print "--out_fasta\t\t\tName of output file you want to generate\n";
    print "--num_sites\t\t\tNumber of sites you want to generate\n";
    print "--kmer\t\t\t\tKmer size - need this to get k-1 bases before and after the SNP position\n";
    print "--chrom\t\t\tSpecify from which chromosome you want to find these SNP sites\n";
    print "--ref_fa\t\t\tFasta of the reference genome (or just the chromosome you care about)\n";
    print "\n\n";
    exit();
}


printf("open vcf\n");

open(FILE, $vcf)||die("Cannot open $vcf");

my $tmpfile = "tmp_chrpos";
my $tmpfile2 = "tmp_alleles";
open(OUT, ">".$tmpfile)||die("Cannot open intermediate/working file with the exciting name tmp_chrpos - directory permission issue?");
open(OUT2, ">".$tmpfile2)||die("Cannot open second intermed file, tmp_alleles\n");
my $sam_column=-1;
my $sites_so_far=0;
while (<FILE>)
{
    my $line = $_;
    chomp $line;
    if ($line =~ /^#/)
    {
	if ($line =~ /CHROM/)
	{
	    my @sp = split(/\s+/, $line);
	    my $i;
	    for ($i=0; $i<scalar(@sp); $i++)
	    {
		if ($sp[$i] eq $sample)
		{
		    $sam_column=$i;
		}
	    }
	}
    }
    else
    {
	## first - by the time we get here we should know the column for our sample
	if ($sam_column==-1)
	{
	    die("Unable to find sample $sample in the header of this vcf $vcf\n");
	}
	my @sp = split(/\s+/, $line);
	my $chr = $sp[0];
	my $pos = $sp[1];
	my $name = $sp[2];
	my $ref = $sp[3];
	my $alt = $sp[4];
	my $gt_field = $sp[$sam_column];
	if ( ($chr eq $chrom) && ($gt_field eq "1/1") )
	{
	    #this is a site where our sample is hom non-ref. Good site for us
	    my $start = $pos-$kmer+1;
	    my $len = 2*$kmer-1;
	    print OUT "$chr:$start:$len\n";
	    print OUT2 "$name\t$ref\t$alt\n";

	    $sites_so_far++;
	    if ($sites_so_far>$num_sites)
	    {
		last;
	    }
	}
    }
    
}
close(FILE);
close(OUT);
close(OUT2);


# Now call Isaac's script
my $tmp = "tmp_fasta";
my $cmd  = "perl $isaac_dir/fastn_substr.pl  --file tmp_chrpos $ref_fa > $tmp";
#print "$cmd\n";
my $ret  = qx{$cmd};
#print "$ret\n";

## Now make the final fasta. Reads come in pairs. First one is the alt allele (sites where sample is hom-alt).
##                                                Second one is ref allele
open(OUT, ">".$out_fasta)||die("Cannot open $out_fasta");
open(TMP, $tmp)||die("Cannot open $tmp");
open(TMP2, $tmpfile2)||die("Cannot open $tmpfile2");

while(<TMP2>)
{
    my $tmp2line = $_;
    chomp $tmp2line;
    my @sp2 = split(/\t/, $tmp2line);
    my $name = $sp2[0];
    my $ref  = $sp2[1];
    my $alt  = $sp2[2];

    ## ok, I now have both alleles.
    my $line = <TMP>;
    print OUT ">REF_ALLELE_".$name."\n"; 
    $line = <TMP>;
    print OUT $line;
    chomp $line;
    my $alt_allele = $line;
    $alt_allele = substr($alt_allele, 0, $kmer-1).$alt.substr($alt_allele, $kmer, $kmer-1);
    print OUT ">ALT_ALLELE_".$name." - expected to be homozygous in $sample\n";
    print OUT "$alt_allele\n";
    
} 
close(TMP2);
close(TMP);
close(OUT);
