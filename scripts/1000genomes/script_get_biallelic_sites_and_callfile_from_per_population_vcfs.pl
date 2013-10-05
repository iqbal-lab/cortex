#!/usr/bin/perl -w
use strict;

my $cortex_root = "/home/zam/dev/CORTEX_release_v1.0.5.21/";  ## set this!

use lib $cortex_root."/scripts/analyse_variants/bioinf-perl/lib/";

my $sorted_vcf = shift;
my $list_orig_vcfs = shift; ##list of the raw VCFs generated per population  - these must contain the SITE_CONF value
my $ref = shift; ##full path to human_g1k_v37.fasta.proper_chroms_only


my $scripts_dir = $cortex_root."scripts/analyse_variants/bioinf-perl/vcf_scripts/";
my $g1k_dir = $cortex_root."scripts/1000genomes/";

my $intermediate_vcf = "combine_all_raw.sorted.perl_filtered.tmp";
my $cmd =  $scripts_dir."vcf_remove_dupes.pl $sorted_vcf  | $scripts_dir"."vcf_remove_overlaps.pl --padding 31 --filter_txt OVERLAP > $intermediate_vcf";
print "$cmd\n";
my $ret = qx{$cmd};
print "$ret\n";

## nearly finished. But now, we have a SITE_CONF for each site from each population. Take the MAX of these, and annotate the VCF with this.
##Also move SITE_CONF from the per-sample fields, where it does not belong (same for all samples) to the INFO field


my %sites = ();
print "Get sites..";
get_sites(\%sites, $intermediate_vcf);
print "..done\n";


my %pos_to_siteconf=();
print "Get site confs..";
get_siteconfs_from_list(\%sites, \%pos_to_siteconf, $list_orig_vcfs);
print "..done\n";

my $out = "cortex_phase3_biallelic_sitelist";

open(OUT, ">".$out)||die();
open(IN, $intermediate_vcf)||die();

my $varnum=0;

while (<IN>)
{
    my $line = $_;
    chomp $line;
    if ($line =~ /^\#/)
    {
	if ($line !~ /CHROM/)
	{
	    print OUT "$line\n";
	}
	else
	{
	    my @head = split(/\s+/, $line);
	    print OUT join("\t", @head[0..7]);
	    print OUT "\n";
	}
    }
    else
    {
	my @sp = split(/\t/, $line);
	my $filter = $sp[6];
	if ($filter eq "PASS")
	{
	    $varnum++;

	    my $chr = $sp[0];
	    my $pos = $sp[1];
	    my $name = "var_".$varnum;
	    $sp[2]=$name;
	    print OUT join("\t", @sp[0..6]);
	    my $ref = $sp[3];
	    my $alt = $sp[4];
	    my $info = $sp[7];
	    print OUT "\t$info:SITE_CONF=";
	    print OUT $pos_to_siteconf{$chr."_".$pos."_".$ref."_".$alt};
	    print OUT "\n";
	}
    }
}
close(IN);
close(OUT);


my $out2 = $out.".annot_flanks";
my $cmd2 = "perl $scripts_dir"."vcf_add_flanks.pl 35 $out g1k_grc37 $ref  > $out2";
print "$cmd2\n$ret2\n";

my $callfile_out = "cortex_phase3_biallelic_pseudo_callfile";
my $cmd3 = "perl $g1k_dir"."build_pseudo_callfile_from_vcf.pl $out2 > $final_out";
my $ret3 = qx{$cmd3};


sub get_siteconfs_from_list
{
    my ($href_sites, $href_pos_to_conf, $list) = @_;
    open(LIST, $list)||die();
    while (<LIST>)
    {
	my $line = $_;
	print "Parse file $line to get site confs..";
	chomp $line;
	get_siteconfs_from_file($href_sites, $href_pos_to_conf, $line);
	print "..done\n";
    }
    close(LIST);
}


sub get_siteconfs_from_file
{
    my ($href_sites, $href_confs, $fyle) = @_;
    open(FYLE, $fyle)||die();
    while (<FYLE>)
    {
	my $line = $_;
	chomp $line;
	if ($line !~ /^\#/)
	{
	    my @sp = split(/\t/, $line);
	    my $chr = $sp[0];
	    my $pos = $sp[1];
	    my $ref = $sp[3];
	    my $alt = $sp[4];
	    my $id = $chr."_".$pos."_".$ref."_".$alt;
	    if (exists $href_sites->{$id})
	    {
		my $field = $sp[9];
		my @sp2 = split(/:/, $field);
		my $conf = $sp2[3];  ##  GT:COV:GT_CONF:SITE_CONF
		if (!exists $href_confs->{$id})
		{
		    $href_confs->{$id}=$conf;
		}
		else
		{
		    if ($conf>$href_confs->{$id})
		    {
			$href_confs->{$id}=$conf;
		    }
		}
	    }
	}
    }
    close(FYLE);
}


sub get_sites
{
    my ($href, $vcf) = @_;

    open(FILE, $vcf)||die();
    
    while (<FILE>)
    {
	my $line = $_;
	chomp $line;
	if ($line =~ /^\#/)
	{
	}
	else
	{
	    my @sp = split(/\t/, $line);
	    my $filter = $sp[6];
	    if ($filter eq "PASS")
	    {
		my $chr = $sp[0];
		my $pos = $sp[1];
		my $ref = $sp[3];
		my $alt = $sp[4];
		$href->{$chr."_".$pos."_".$ref."_".$alt}=1;
	    }
	}
    }
    close(FILE);

}

