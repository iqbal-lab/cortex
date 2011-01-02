#!/usr/bin/perl -w
use strict;
use File::Basename;
use Descriptive;


# Takes as input a set of variant calls made by Cortex, maps the flanks, aligns
# branches against each other to determine variant type,
# splits out SNPs from clusters so we have precise loci,
# applies any of the following filters
#  - demand 5prime flank maps with quality> some threshold
#  - for ref-assisted calls, insist that the normalised covg on the ref-allele is close to zero.
# dumps to VCF4.0

my $callfile = shift;
my $algo = shift; #"bubble", or "pd" (for path-divergence)
my $type = shift; # het or hom 
#my $d_eff = shift;#effective diploid covg, can enter -1 if looking at bubble calls, as do not use it
my $prefix = shift; #give prefix to make varnames globally unique
my $outdir = shift;
my $outvcf_filename_stub = shift;
my $mercurial_changeset = shift;#used for making the calls --> vcf header
my $indiv = shift;
my $reference = shift; #text to go into VCF header to describe this callset

### initialisation
my $stampy_bin = "stampy.py";
my $covg_filter_type = "median"; #median covg on both branches must be >= some threshold
my $covg_filter_thresh = 2; 
my $stampy_hash_stub = "stampy_reference"; #name for reference hash file created by Stampy
my $flank_bin = "make_5p_flank_file.pl";
my $process_bubbles_bin = "process_bubbles.pl";

########
#### checks
########
if ( ($type ne "het") && ($type ne "hom") )
{
    die("Type must be het or hom");
}


if (!(-e $flank_bin) )
{
    die("Cannot find $flank_bin");
}

if (!(-e $process_bubbles_bin))
{
    die("Cannot find $process_bubbles_bin");
}


if (!(-e $callfile))
{
    die("Cannot find $callfile");
}



if ($outdir !~ /\/$/)
{
    $outdir = $outdir.'/';
}








## 1. Map 5p flanks
my $bname = basename($callfile);
my $flankfile = $outdir.$bname.".5pflanks";

my $flank_cmd = "perl $flank_bin $callfile > $flankfile";
print "$flank_cmd\n";
my $flank_ret = qx{$flank_cmd};
print "$flank_ret\n";


if ( (!(-e $flankfile)) || (-z $flankfile))
{
    die("$flankfile either failed to be created, or was created with size zero");
}

my $mapped_flanks = $flankfile.".sam";
if (-e $mapped_flanks)
{
    print ("$mapped_flanks already exists, so wont re-do it");
}
else
{
    my $map_flanks_cmd = "$stampy_bin -g $stampy_hash_stub -h $stampy_hash_stub --norefoutput --inputformat=fasta -M $flankfile -o $mapped_flanks";
    print "$map_flanks_cmd\n";
    my $map_flanks_ret = qx{$map_flanks_cmd};
    print "$map_flanks_ret\n";
}

my %var_name_to_cut_flank=();
get_list_vars_with_cut_flanks($mapped_flanks, \%var_name_to_cut_flank);


## 2. Align branches against each other
my $proc_bub_output = $outdir.$bname.".aligned_branches";

if (!(-e $proc_bub_output) )
{
    my $proc_bub_cmd = "perl $process_bubbles_bin $callfile $prefix > $proc_bub_output 2>&1";
    print "$proc_bub_cmd\n";
    my $proc_bub_ret = qx{$proc_bub_cmd};
    print "$proc_bub_ret\n";
}
else
{
    print "Not aligning branches against each other, as $proc_bub_output already exists\n";
}


##  3. Apply filters, and collect a list of good calls.

my %var_name_to_flank_mq_filter=();
my %var_name_to_PD_normalisation_filter=();
my %var_name_to_covg_and_branch_filter=();
my %var_name_to_combined_filtering_result=();


my $filter_check_one_br_in_ref = 0;
if ($type eq "het")
{
    $filter_check_one_br_in_ref  = 1;
}



my $apply_pd_filters=0;
if ($algo eq "pd")
{
    $apply_pd_filters=1;
}

filter_by_coverage_and_compare_branches_with_ref($callfile, $filter_check_one_br_in_ref, $covg_filter_type, $covg_filter_thresh, $type, \%var_name_to_covg_and_branch_filter,
						 $apply_pd_filters);
filter_by_flank_mapqual($mapped_flanks, \%var_name_to_flank_mq_filter);



combine_all_filters(\%var_name_to_covg_and_branch_filter, \%var_name_to_flank_mq_filter, \%var_name_to_combined_filtering_result);


my $fh_calls;
my $fh_map_flanks;
my $fh_proc_bub;
my $fh_simple_vcf;
my $fh_decomp_vcf;

open($fh_calls, $callfile)||die("Cannot open $callfile");
open($fh_map_flanks, $mapped_flanks)||die("Cannot open $mapped_flanks");
open($fh_proc_bub, $proc_bub_output)||die("Cannot open $proc_bub_output");
my $simple_vcf_name = $outdir.$outvcf_filename_stub.".raw.vcf";
open($fh_simple_vcf, "> $simple_vcf_name")||die("Cannot open ");
my $decomp_vcf_name = $outdir.$outvcf_filename_stub.".decomp.vcf";
open($fh_decomp_vcf, "> $decomp_vcf_name")||die("Cannot open $decomp_vcf_name");

##print vcf header
my $header = get_vcf_header($mercurial_changeset, $reference, $indiv);
print $fh_simple_vcf $header;
print $fh_decomp_vcf $header;


my $ret=1;
while ($ret==1)
{

    $ret = print_next_vcf_entry_for_easy_and_decomposed_vcfs($type, $fh_calls, $fh_map_flanks, $fh_proc_bub,
							     \%var_name_to_combined_filtering_result, \%var_name_to_cut_flank,
							     1,1,
							     $fh_simple_vcf, $fh_decomp_vcf);
}
close($fh_calls);
close($fh_map_flanks);
close($fh_proc_bub);



sub get_vcf_header
{
    my ($changeset, $reference, $header) = @_;

    my $date_cmd = "date \'+\%d\/\%m\/\%y\'";
    my $date = qx{$date_cmd};
    chomp $date;

    my $head = "";
    
    $head = $head. "##fileformat=VCFv4.0\n";
    $head = $head. "##fileDate=$date\n";
    $head = $head. "##source=Cortex changeset $changeset\n";
    $head = $head. "##reference=$reference\n";
    $head = $head. "##phasing=none, though some calls involve phasing clustered variants\n";
    $head = $head. "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    $head = $head. "##FORMAT=<ID=COV,Number=2,Type=Integer,Description=\"Median coverage on ref and alt alleles\">\n";
    $head = $head. "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n";
    $head = $head. "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of variant\">\n";
    $head = $head. "##ALT=<ID=SNP,Description=\"SNP\">\n";
    $head = $head. "##ALT=<ID=SNP_FROM_COMPLEX,Description=\"SNP called from a cluster of phased SNPs or complex SNP/indel , split out for easier comparison with other SNP call sets\">\n";
    $head = $head. "##ALT=<ID=DEL,Description=\"Deletion\">\n";
    $head = $head. "##ALT=<ID=INS,Description=\"Insertion of novel sequence\">\n";
    $head = $head. "##ALT=<ID=INDEL,Description=\"Insertion-deletion\">\n";
    $head = $head. "##ALT=<ID=INV,Description=\"Inversion\">\n";
    $head = $head. "##ALT=<ID=INV_INDEL,Description=\"Inversion+indel\">\n";
    $head = $head. "##ALT=<ID=DEL_INV,Description=\"Deletion + Inversion\">\n";
    $head = $head. "##ALT=<ID=INS_INV,Description=\"Insertion + Inversion\">\n";
    $head = $head. "##ALT=<ID=PH_SNPS,Description=\"Phased SNPs\">\n";
    $head = $head. "##ALT=<ID=COMPLEX, Description=\"Complex variant, collection of SNPs and indels\">\n";
    $head = $head. "##FILTER=<ID=flank_mq,Description=\"5prime flank maps to reference with mapping quality below 30\">\n";
    $head = $head. "##FILTER=<ID=allele_cov,Description=\"Median coverage on >=1 allele below 2\">\n";
    $head = $head. "##FILTER=<ID=path_div,Description=\"Nodes in ref allele of Path-divergence hom-non-ref call which are on the ref allele but not alt allele have median covg>0\">\n";
    $head = $head. "##FILTER=<ID=path_div,Description=\"All nodes in ref allele of Path-divergence call occur more than once in reference\">\n";
    $head =  $head."#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$indiv\n";
    return $head;
}




#if het, require both branches have median covg >= filter threshold
#if hom non-ref, require the non-ref branch has median covg >=filter threwshold
#if apply_pd_filt==1, then must be a hom nonref call, and add the following filters
## 1. Covg on nodes that are on the trusted allele but not the var allele, AND are unique in the reference,  should be zero
##    We just take the median covg on trusted-non-var nodes. If there are no trusted non-var nodes, reject the call.
## 2. Fraction of trusted nodes  unique in ref >0?? - make this a parameter of the function 
sub filter_by_coverage_and_compare_branches_with_ref
{
    my ($file, $should_we_check_one_br_in_ref, $filter_type, $filter_threshold, $hom_or_het, $href, $apply_pd_filt)=@_;

    if ($filter_type ne "median")
    {
	die("Have only implemented median filter");
    }

    if ( ($apply_pd_filt==1) && ($hom_or_het ne "hom") )
    {
	die("Cannot call filter_by_coverage_and_compare_branches_with_ref asking to apply pd filters if also say is not a hom call dataset");
    }


    open(FILE, $file)||die("Cannot open $file");

    while (<FILE>)
    {
	my $line = $_;
	chomp $line;
	if ($line =~ /(var_\d+)_5p_flank/) 
	{
	    my $pass_filter="true";
	    my $varname = $prefix."_".$1;
	    my @trusted_non_var=();
	    my @var_non_trusted=();
	    $line = <FILE>;# 5p sequence - need to check is not too small
	    chomp $line;
	    if (length($line)<25)
	    {
		print "$varname fails filter as 5p flank <25bp - will not be able to map it uniquely. Length is ";
		print length($line);
		print "\n";
		$pass_filter="false";
	    }
	    $line = <FILE>;##br1 read_id
	    if ($apply_pd_filt==1)
	    {
		if ($line =~ /covgs of trusted not variant nodes:\s+([\d\s]+)number_of_such_nodes/)
		{
		    my $str = $1;
		    if ($str =~ /(.+)\s+$/)
		    {
			$str=$1;
		    }
		    if ($str =~ /^\s+(.+)/)
		    {
			$str=$1;
		    }

		    @trusted_non_var = split(/\s/, $str);##this gives us an array containing the coverages of trusted non-var nodes
		}
	    }
	    <FILE>;#ignore br1 seq
	    $line = <FILE>;# br2 read id
	    if ($apply_pd_filt==1)
	    {
		if ($line =~ /covgs of variant not trusted nodes:\s+([\d\s]+)number_of_such_nodes/)
		{
		    my $str=$1;
		    if ($str =~ /(.+)\s+$/)
                    {
                        $str=$1;
                    }
                    if ($str =~ /^\s+(.+)/)
                    {
                        $str=$1;
                    }
		    @var_non_trusted = split(/\s/, $str);
		}
	    }
	    <FILE>;#ignore br1 seq
	    <FILE>;#ignore 3p flank read id
	    $line = <FILE>;# 3p flank seq
	    chomp $line;
	    if (length($line)<3)
	    {
		print "$varname fails filter, as 3p flank <3bp\n";
		$pass_filter="false";
	    }

	    $line = <FILE>;
	    if ($line =~ /extra information/)
	    {
		$line = <FILE>;
	    }
	    $line = <FILE>;

	    if ($line !~ /branch1 coverages/)
	    {
		die("Zam Expected to see \"branch 1 coverages\" but instead saw $line");
	    }
	    $line = <FILE>;
	    if ($line !~ /Mult in  hum ref/)
	    {
		die("Expected to see \"Mult in  hum ref\" but instead saw $line");
	    }
	    $line = <FILE>;
	    chomp $line;
	    my @br1_refmult=split(/\s+/, $line);


	    $line = <FILE>;
	    if ($line !~ /Covg in indiv/)
	    {
		die("Expected to see \"Covg in indiv\" but instead saw $line");
	    }
	    $line = <FILE>;
	    chomp $line;
	    my @br1_cov = split(/\s+/, $line);

	    <FILE>;
	    $line = <FILE>;
	    if ($line !~ /Mult in  hum ref/)
	    {
		die("Second time - Expected to see \"Mult in  hum ref\" but instead saw $line");
	    }
	    $line = <FILE>;
	    chomp $line;
	    
	    my @br2_refmult = split(/\s+/, $line);

	    $line = <FILE>;
	    if ($line !~ /Covg in indiv/)
	    {
		die("Second branch - Expected to see \"Covg in indiv\" but instead saw $line");
	    }
	    $line = <FILE>;
	    chomp $line;
	    my @br2_cov = split(/\s+/, $line);


	    my $stat_br1_humref = Statistics::Descriptive::Full->new();
	    $stat_br1_humref->add_data(@br1_refmult);
	    my $stat_br2_humref = Statistics::Descriptive::Full->new();
	    $stat_br2_humref->add_data(@br2_refmult);
	    my $stat_br1_cov =  Statistics::Descriptive::Full->new();
	    $stat_br1_cov -> add_data(@br1_cov);
	    my $stat_br2_cov =  Statistics::Descriptive::Full->new();
	    $stat_br2_cov -> add_data(@br2_cov);
	    
	    my $min_refmult_br1 = $stat_br1_humref->min();
	    my $min_refmult_br2 = $stat_br2_humref->min();

	    my $median_br1_cov = $stat_br1_cov->median();
	    my $median_br2_cov = $stat_br2_cov->median();



	    ## check one branch is in reference 
	    if ($should_we_check_one_br_in_ref==1)
	    {
		##one of the two branches must have min multiplicity in the ref of  >=1
		if ( ($min_refmult_br1==0) && ($min_refmult_br2==0) )
		{
		    print "$varname Fails at this point - het call, but not true that one branch is in reference\n";
		    $pass_filter="false";
		}
	    }
	    if ($pass_filter eq "true")
	    {
		

		## for a het, apply condition to both branches
		if ($hom_or_het eq "het")
		{
		    if ( ($median_br1_cov>=$filter_threshold) && ($median_br2_cov>=$filter_threshold) )
		    {
			$pass_filter = "true";
		    }
		    else
		    {
			$pass_filter="false";
			print "$varname Fails branch filter, median covgs on br1 and br2 are $median_br1_cov and $median_br2_cov and threshold is $filter_threshold\n";
		    }
		}
		else
		{
		    #for a hom, determine the nonref branch, and apply to that
		    my $which_br_is_ref;

		    if ( ($min_refmult_br1>=1) && ($min_refmult_br2==0) )
		    {
			$which_br_is_ref=1;
		    }
		    elsif ( ($min_refmult_br2>=1) && ($min_refmult_br1==0) )
		    {
			$which_br_is_ref=2;
		    }
		    elsif ( ($min_refmult_br1>=1) && ($min_refmult_br2>=1) )
		    {
			$pass_filter="false";
			print("Hom call $varname has both branches in reference - fails filter");
			$which_br_is_ref="b";
		    }
		    else
		    {
			$which_br_is_ref="neither";
			$pass_filter="false";
                        print("Hom call $varname has neither branch in reference - fails filter");

		    }

		    if ($pass_filter eq "true")##in which case which_br_is_ref is numeric
		    {
			if ( ($which_br_is_ref==1) && ($median_br2_cov < $filter_threshold) )
			{
			    $pass_filter="false";
			    print "Hom call Fails branch filter, median covgs on alt allele, br2 is $median_br2_cov and threshold is $filter_threshold\n";
			    
			}
			elsif ( ($which_br_is_ref==2) && ($median_br1_cov < $filter_threshold) )
			{
			    $pass_filter="false";
			    print "Hom call Fails branch filter, median covgs on br1 is $median_br1_cov and threshold is $filter_threshold\n";
			    
			}
			
			
		    }
		}


		if ($apply_pd_filt==1)
		{
		    if (scalar(@trusted_non_var)==0)
		    {
			$pass_filter="false";
			print "$varname fails PD filter; all trusted nodes poccur >1 time in reference\n";
		    }
		    else
		    {

			my $stat_trusted_non_var = Statistics::Descriptive::Full->new();
			$stat_trusted_non_var->add_data(@trusted_non_var);
			my $tr_non_var_med = $stat_trusted_non_var->median();

			if ($tr_non_var_med >0 )
			{
			    $pass_filter="false";
			    print "$varname fails PD filter: median covg on trusted_non_var nodes is $tr_non_var_med, which is >0";
			}

		    }

		}




	    }



	    if ($pass_filter eq "true")
	    {
		$href->{$varname} = "PASS";
	    }
	    else
	    {
		$href->{$varname} = "FAIL";
	    }
	}
	else
	{

	}
    }
    close(FILE);


}




sub filter_by_flank_mapqual
{
    my ($file, $href) = @_;

    open(FILE, $file)||die();
    
    while(<FILE>)
    {
	my $line = $_;
	
	if ($line =~ /^@/)
	{
	}
	else
	{
	    my @sp = split(/\t/, $line);
	    
	    if (scalar @sp < 4)
	    {
		die("problem parsing $line");
	    }
	    
	    my $query = $sp[0];
	    my $varname;
	    if ($query =~ /(var_\d+)_5p_flank/)
	    {
		$varname = $prefix."_".$1;
		
	    }
	    else
	    {
		die("Expected query name in sam file would be of form var_1_5p_flank, but is $query");
	    }	
	    if ($sp[4]>=30)
	    {
		$href->{$varname}="PASS";
	    }
	    else
	    {
		$href->{$varname} = "FAIL";
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

    my ($type, $file_handle_calls, $file_handle_map_flanks, $file_handle_proc_bubble_output,
	$href_varname_to_passing_all_filters, $href_var_name_to_cut_flank,
	$print_easy_vcf, $print_decomp_vcf, $fh_easy, $fh_decomp) = @_;



    ##1. Get next var info from all three files
    my ($eof, $var_name, $flank5p, $br1_seq, $br2_seq, $flank3p, $br1_med_cov, $br2_med_cov, $which_is_ref) = get_next_var_from_callfile($file_handle_calls);
    
    my ($eof2, $var_name2, $strand, $chr, $coord) = get_next_var_from_flank_mapfile($file_handle_map_flanks);

    my ($eof3, $var_name3, $align_num_bases_agreement_at_start, $align_num_bases_agreement_at_end, $align_direction,
	$num_snps_align, $num_indels_align,$alignment_br1, $alignment_middle, $alignment_br2, $possible_inversion, 
	$clean_indel, $aref_snp_coords, $aref_snp_alleles,
        $aref_coords_of_indels, $aref_alleles_of_indels, $aref_indel_needs_extra_base)

	= get_next_alignment_from_procfile($file_handle_proc_bubble_output, $strand, $which_is_ref);


    ## if niether allele is the ref allele - ignore this
    if ($which_is_ref eq "neither")
    {
	print "Ignore var $var_name, both alleles non ref!\n";
	return 1;
    }
    elsif ($which_is_ref eq "b")
    {
	print "Ignore var $var_name, both alleles REF!\n";
	return 1;
    }


    ## if one but not all reach the end of file
    if (  ( ($eof eq"EOF") || ($eof2 eq "EOF") || ($eof3 eq "EOF") )  &&    (!( ($eof eq "EOF") && ($eof2 eq "EOF") && ($eof3 eq "EOF") )) )
    {
	die("SOme of these files end before the others");
    }
    ##if all hit the end of file
    if ( ($eof eq "EOF") && ($eof2 eq "EOF") && ($eof3 eq "EOF") )
    {
	return 0;
    }


	
    if ( ($var_name ne $var_name2) || ($var_name ne $var_name3) )
    {
	die("Some ordering problem, next var in the three files are $var_name, $var_name2, $var_name3");
    }


    if ( (exists $href_varname_to_passing_all_filters->{$var_name}) 
	&&
	($href_varname_to_passing_all_filters->{$var_name} eq "PASS")
	)
    {
	if ($print_easy_vcf==1)
	{
	    my $split_phased_snps=0;
	    my @empty=();
	    print_vcf_entry( $type,$fh_easy, $var_name, $flank5p, $br1_seq, $br2_seq, $flank3p,
			     $strand, $chr, $coord, $which_is_ref, 
			     $align_num_bases_agreement_at_start, $align_num_bases_agreement_at_end, $align_direction,
			     $num_snps_align, $num_indels_align,$alignment_br1, $alignment_middle, $alignment_br2,
			     $br1_med_cov, $br2_med_cov, $href_var_name_to_cut_flank, $possible_inversion, $clean_indel, $split_phased_snps,
			     0,0,0,0,\@empty) ;
	}
	if ($print_decomp_vcf==1)
	{
	    my $split_phased_snps=1;
	    print_vcf_entry( $type,$fh_decomp, $var_name, $flank5p, $br1_seq, $br2_seq, $flank3p,
			     $strand, $chr, $coord, $which_is_ref, 
			     $align_num_bases_agreement_at_start, $align_num_bases_agreement_at_end, $align_direction,
			     $num_snps_align, $num_indels_align,$alignment_br1, $alignment_middle, $alignment_br2,
			     $br1_med_cov, $br2_med_cov, $href_var_name_to_cut_flank, $possible_inversion, $clean_indel, $split_phased_snps, 
			     $aref_snp_coords, $aref_snp_alleles,
			     $aref_coords_of_indels, $aref_alleles_of_indels, $aref_indel_needs_extra_base) ;


	}
    }
    else
    {
	#print "Fails filter: $var_name\n";
    }
    return 1;

}
		


sub print_vcf_entry
{
    my ( $hom_or_het, $fh_output_vcf, $var_name, $flank5p, $br1_seq, $br2_seq, $flank3p,
	 $strand, $chr, $coord, $which_is_ref,
	 $align_num_bases_agreement_at_start, $align_num_bases_agreement_at_end, $align_direction,
	 $num_snps_align, $num_indels_align,$alignment_br1, $alignment_middle, $alignment_br2,
	 $br1_med_cov, $br2_med_cov, $href_var_name_to_cut_flank, $poss_inversion, $clean_indel,
	 $split_phased_snps, $aref_snp_coords, $aref_snp_alleles,
	 $aref_indel_coords, $aref_indel_alleles, $aref_indel_needs_extra_base) = @_;

    my $vcf_entry_chr = $chr;
    my $vcf_entry_pos;
    my $vcf_entry_ref_allele;
    my $vcf_entry_alt_allele;
    my $error;

    ($vcf_entry_pos, $vcf_entry_ref_allele, $vcf_entry_alt_allele, $error)
	 = get_simple_vcf_entry_pos_and_alleles($var_name, $strand, $coord, $flank5p, $br1_seq, $br2_seq, $flank3p,
						$which_is_ref,
						$align_num_bases_agreement_at_start, $align_num_bases_agreement_at_end, $align_direction,
						$num_snps_align, $num_indels_align, $href_var_name_to_cut_flank);

    if ($error ne "0")
    {
	print ("Ignore this $var_name - due to this error $error\n");
	return;
    }

    ## quick check

    my $count_num_indels_needing_extra_base=0;
    foreach my $in (@$aref_indel_needs_extra_base)
    {
	if ($in==1)
	{
	    $count_num_indels_needing_extra_base++;
	}
    }
    if ($count_num_indels_needing_extra_base>1)
    {
	die("More than one indel in $var_name nees an extra base - namely $count_num_indels_needing_extra_base");
    }


    my $genotype;
    if ($hom_or_het eq "het")
    {
	$genotype="0/1";
    }
    else
    {
	$genotype="1/1";
    }
    my $cov;
    if ($which_is_ref==1)
    {
	$cov=int($br1_med_cov).",".int($br2_med_cov);
    }
    else
    {
	$cov=int($br2_med_cov).",".int($br1_med_cov);
    }
    my $svlen = length($vcf_entry_ref_allele)-length($vcf_entry_alt_allele);
    my $svtype;

    if ( ($num_snps_align==1) && ($num_indels_align==0) )
    {
	$svtype="SNP";
    }
    elsif ($poss_inversion==1)
    {
	$svtype="INV";
    }
    elsif ($poss_inversion==2)
    {
	$svtype="INV_INDEL";
    }
    elsif ( ($num_snps_align>1) && ($num_indels_align==0) ) ## this can only happen if the two branches are the same length :-)
    {
	$svtype="PH_SNPS";
    }
    elsif ( ($clean_indel==1) && ($svlen<0) )
    {
	$svtype="INS";
    }
    elsif ( ($clean_indel==1) && ($svlen>0) )
    {
	$svtype="DEL";
    }
    elsif ( ($num_snps_align==0) && ($num_indels_align==1)  )
    {
	$svtype="INDEL";
    }
    else
    {
	#$svtype = "PH_SNPS_INDELS";
	$svtype = "COMPLEX";
    }

    my $info = "SVTYPE=$svtype;SVLEN=$svlen";

    if ( ( ($svtype  !~ /COMPLEX/) && ($svtype !~ /PH_SNPS/)  )|| ($split_phased_snps==0) )
    {
	print $fh_output_vcf "$vcf_entry_chr\t$vcf_entry_pos\t$var_name\t$vcf_entry_ref_allele\t$vcf_entry_alt_allele\t.\tPASS\t$info\tGT:COV\t$genotype:$cov\n";
    }
    else
    {

	## double check
	if (scalar(@$aref_snp_coords) != scalar(@$aref_snp_alleles) )
	{
	    print("SNP coord and Allele arrays of different lengths: ");
	    print scalar(@$aref_snp_coords) ;
	    print " ";
	    print scalar(@$aref_snp_alleles);
	    die("\n");
	}


	## we have coordinates of phased SNPs with respect to fw_middle.
	## we will inherit the $genotype and $cov data from the overall variant. 
	## So all that changes are: POS, NAME=var_name_snp_NUMBER REF, ALT, sv_type=SNP_FROM_PHASED_SNP_CALL, svlen=0,genotype,
	my $cnt=0;
	foreach my $c (@$aref_snp_coords)
	{
	    $cnt++;
	    my $two_alleles = $aref_snp_alleles->[$cnt-1];
	    my $this_snp_ref_allele;
	    my $this_snp_alt_allele;
	    if ($two_alleles =~/^([ACGT])_([ACGT])$/)
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

	    my $this_snp_name = $var_name."_sub_snp_".$cnt;
	    my $this_snp_info = "SVTYPE=SNP_FROM_COMPLEX;SVLEN=0";
	    print $fh_output_vcf "$this_snp_chr\t$this_snp_pos\t$this_snp_name\t$this_snp_ref_allele\t$this_snp_alt_allele\t.\tPASS\t$this_snp_info\tGT:COV\t$genotype:$cov\n";
	}
	$cnt=0;
	my $index_of_indel=0;
	foreach my $c (@$aref_indel_coords)
	{
	    $cnt++;
	    my $needs_extra_base = $aref_indel_needs_extra_base->[$index_of_indel];
            my $two_alleles = $aref_indel_alleles->[$cnt-1];
            my $this_indel_ref_allele;
            my $this_indel_alt_allele;
            if ($two_alleles =~/^([ACGT]+)_([ACGT]+)$/)
            {
                $this_indel_ref_allele = $1;
                $this_indel_alt_allele = $2;
            }
	    elsif ($two_alleles =~/^_([ACGT]+)$/)
            {
                $this_indel_ref_allele = "";
                $this_indel_alt_allele = $1;
            }
	    elsif ($two_alleles =~/^([ACGT]+)_$/)
            {
                $this_indel_ref_allele = $1;
                $this_indel_alt_allele = "";
            }
            else
            {
                die("Bad formatting for indel alleles: $two_alleles");
            }

	    if ($needs_extra_base==1)
	    {
		$this_indel_ref_allele=substr($vcf_entry_ref_allele, 0,1).$this_indel_ref_allele;
		$this_indel_alt_allele=substr($vcf_entry_alt_allele, 0,1).$this_indel_alt_allele;
	    }

            #we have already got the snp coords with respect to the start of the overall-ref allele, AND the ref/alt alleles for each snp
            #in the right direction

            my $this_indel_chr = $vcf_entry_chr;

            ## take care here. $vcf_entry_pos is the right coordinate if we are going to call an indel for the whole shebang - ie one base BEFORE the first variant.
	    ## hence add 1, as $c has already also subtracted 1
            my $this_indel_pos= $vcf_entry_pos + $c +1; 

	     
            my $this_indel_name = $var_name."_sub_indel_".$cnt;
	    my $svlen = length($this_indel_ref_allele)-length($this_indel_alt_allele);
            my $this_indel_info = "SVTYPE=INDEL_FROM_COMPLEX;SVLEN=$svlen";
            print $fh_output_vcf "$this_indel_chr\t$this_indel_pos\t$this_indel_name\t$this_indel_ref_allele\t$this_indel_alt_allele\t.\tPASS\t$this_indel_info\tGT:COV\t$genotype:$cov\n";


	    $index_of_indel++;
	}
    }
					 
}					 
	





sub get_simple_vcf_entry_pos_and_alleles

{
 
    my ($name, $str, $pos, $flank5p, $br1_seq, $br2_seq, $flank3p, $which_br_ref, 
	$align_num_bp_agreement_at_start, $align_num_bp_agreement_at_end, $align_dir,
	$num_snps_align, $num_indels_align, $href_var_name_to_cut_flank) = @_;

    ## if the 5p flank is >1000bp, then stampy fails, so we cut the 5p flank and only take the last 1000bp.
    if (exists $href_var_name_to_cut_flank->{$name})
    {
	##the flank5p passedin is directly taken from the callfile
	$flank5p = substr($flank5p, -1000);
	if (length($flank5p)!=1000)
	{
	    die("perl issue with $flank5p");
	}
    }




    ## snp for single SNP, else not_snp
    my $what_type = determine_type($num_snps_align, $num_indels_align, length($br1_seq), length($br2_seq) );

    my $vcf_entry_ref_allele;
    my $vcf_entry_alt_allele;
    my $vcf_entry_pos;
    ## if 5p flank maps in + direction
    if ($str==0)
    {
	##bloody vcf format wants the ref allele to start at the first variant base for a SNP, and the base BEFORE that for a non-snp
	if ($what_type eq "snp")
	{
	    $vcf_entry_pos        =  $pos+length($flank5p);
	}
	else
	{
	    $vcf_entry_pos        =  $pos+length($flank5p)-1;
	}
	
	if ($which_br_ref==1)
	{
	    ## if br1 and br2 are the same at the end, cut off that bit

	    if ($what_type eq "snp")
	    {
		$vcf_entry_ref_allele =  substr($br1_seq, 0, length($br1_seq) - $align_num_bp_agreement_at_end);
		$vcf_entry_alt_allele =  substr($br2_seq, 0, length($br2_seq) - $align_num_bp_agreement_at_end);
	    }
	    else
	    {
		$vcf_entry_ref_allele =  substr($flank5p, -1).substr($br1_seq, 0, length($br1_seq) - $align_num_bp_agreement_at_end);
		$vcf_entry_alt_allele =  substr($flank5p, -1).substr($br2_seq, 0, length($br2_seq) - $align_num_bp_agreement_at_end);
	    }
	}
	elsif ($which_br_ref==2)
	{
	    ## if br1 and br2 are the same at the end, cut off that bit

	    if ($what_type eq "snp")
	    {
		$vcf_entry_ref_allele =  substr($br2_seq, 0, length($br2_seq) - $align_num_bp_agreement_at_end);
		$vcf_entry_alt_allele =  substr($br1_seq, 0, length($br1_seq) - $align_num_bp_agreement_at_end);
	    }
	    else
	    {
		$vcf_entry_ref_allele =  substr($flank5p, -1).substr($br2_seq, 0, length($br2_seq) - $align_num_bp_agreement_at_end);
		$vcf_entry_alt_allele = substr($flank5p, -1).substr($br1_seq, 0, length($br1_seq) - $align_num_bp_agreement_at_end);
	    }
	}
	else
	{
	    return (0,0,0,"This var has both alleles entirely in the reference: $which_br_ref");
	}
    }
    ##else, 5prime flank mapped in the reverse direction
    elsif ($str==16)
    {
	my $br1_excepting_agreement_at_end = substr($br1_seq, 0,length($br1_seq)-$align_num_bp_agreement_at_end);
	my $br2_excepting_agreement_at_end = substr($br2_seq, 0,length($br2_seq)-$align_num_bp_agreement_at_end);

	my $last_base_in_branch_before_variant;
	if ($align_num_bp_agreement_at_end>0)
	{
	    $last_base_in_branch_before_variant = substr($br1_seq, length($br1_seq)-$align_num_bp_agreement_at_end,1);
	}



	if ($which_br_ref==1)
	{
	    ## if br1 and br2 are the same at the end, cut off that bit

	    if ($what_type eq "snp")
	    {
		$vcf_entry_pos        =  $pos-length($br1_seq) + $align_num_bp_agreement_at_end ;
	    }
	    else
	    {
		$vcf_entry_pos        =  $pos-length($br1_seq)-1 + $align_num_bp_agreement_at_end;
	    }


	    if ($what_type eq "snp")
	    {
		$vcf_entry_ref_allele =  rev_comp($br1_excepting_agreement_at_end);
		$vcf_entry_alt_allele =  rev_comp($br2_excepting_agreement_at_end);
		#$vcf_entry_ref_allele =  rev_comp(substr($br1_seq, - (length($br1_seq) - $align_num_bp_agreement_at_start) ) );
		#$vcf_entry_alt_allele =  rev_comp(substr($br2_seq, - (length($br2_seq) - $align_num_bp_agreement_at_start) ) );
	    }
	    else
	    {
		if ($align_num_bp_agreement_at_end==0)
		{
		    $vcf_entry_ref_allele = substr( rev_comp($flank3p), -1)  .rev_comp($br1_excepting_agreement_at_end);
		    $vcf_entry_alt_allele = substr( rev_comp($flank3p), -1)  .rev_comp($br2_excepting_agreement_at_end);
		}
		else
		{
		    $vcf_entry_ref_allele = rev_comp($last_base_in_branch_before_variant).rev_comp($br1_excepting_agreement_at_end);
		    $vcf_entry_alt_allele = rev_comp($last_base_in_branch_before_variant).rev_comp($br2_excepting_agreement_at_end);
		}
		#$vcf_entry_ref_allele =  substr( rev_comp($flank3p), -1)  .rev_comp(substr($br1_seq, - (length($br1_seq) - $align_num_bp_agreement_at_start) ) );
		#$vcf_entry_alt_allele =  substr( rev_comp($flank3p), -1)  .rev_comp(substr($br2_seq, - (length($br2_seq) - $align_num_bp_agreement_at_start) ) );
	    }
	}
	elsif ($which_br_ref==2)
	{
	    ## if br1 and br2 are the same at the end, cut off that bit

	    if ($what_type eq "snp")
	    {
		$vcf_entry_pos        =  $pos-length($br2_seq) + $align_num_bp_agreement_at_end;
	    }
	    else
	    {
		$vcf_entry_pos        =  $pos-length($br2_seq)-1 + $align_num_bp_agreement_at_end;
	    }


	    if ($what_type eq "snp")
	    {
		$vcf_entry_ref_allele =  rev_comp($br2_excepting_agreement_at_end);
		$vcf_entry_alt_allele =  rev_comp($br1_excepting_agreement_at_end);
	    }
	    else
	    {
		if ($align_num_bp_agreement_at_end==0)
                {
		    $vcf_entry_ref_allele = substr( rev_comp($flank3p), -1)  .rev_comp($br2_excepting_agreement_at_end);
		    $vcf_entry_alt_allele = substr( rev_comp($flank3p), -1)  .rev_comp($br1_excepting_agreement_at_end);
		}
		else
		{
		    $vcf_entry_ref_allele = rev_comp($last_base_in_branch_before_variant).rev_comp($br2_excepting_agreement_at_end);
		    $vcf_entry_alt_allele = rev_comp($last_base_in_branch_before_variant).rev_comp($br1_excepting_agreement_at_end);
		}

	    }
	}
	else
	{
	    return (0,0,0,"This var has both alleles entirely in the reference: $which_br_ref");	    
	}

	
    }
    else
    {
	die("Given str $str\n");
    }

    return ($vcf_entry_pos, $vcf_entry_ref_allele, $vcf_entry_alt_allele, "0");
    
}				 
						     

sub rev_comp
{
  my ($seq) = @_;

  my $r_seq = reverse($seq);
  $r_seq =~ tr/acgtACGT/tgcaTGCA/;
 # print join(" ",$seq,$r_seq),"\n";
  return $r_seq;

}

sub determine_type
{
    my ($num_snps_align, $num_indels_align, $lenbr1, $lenbr2) = @_;

    my $what_type;
    if ( ($num_snps_align==1) && ($num_indels_align==0) && ($lenbr1==$lenbr2) )
    {
	$what_type="snp";
    }
    else
    {
	$what_type="not_snp"; 
    }

    return $what_type;
}




sub get_next_var_from_callfile
{
    my ($fh) = @_;

    my $line = "";
    my $varname;
    my $flank5p;
    my $br1;
    my $br2;
    my $flank3p;
    
    while ($line !~  /(var_\d+)_5p_flank/) 
    {	
	if (eof($fh))
	{
	    return ("EOF", 0,0,0,0,0,0,0,0);
	}
	else
	{
	    $line = <$fh>;
	}
	
    }
    if ($line =~ /(var_\d+)_5p_flank/)
    {
	$varname = $prefix."_".$1;
	
	$flank5p=<$fh>;
	chomp $flank5p;
	<$fh>;#ignore br1 read id
	$br1 = <$fh>;
	chomp $br1;
	<$fh>;#ignore br2 read id
	$br2=<$fh>;
	chomp $br2;
	<$fh>;#ignore 3p flank read id
	$flank3p= <$fh>;
	$line = <$fh>;
	if ($line =~ /extra information/)
	{
	    <$fh>;
	}
	$line = <$fh>;
	if ($line !~ /branch1 coverages/)
	{
	    die("Expected to see \"branch 1 coverages\" but instead saw $line");
	}
	$line = <$fh>;
	if ($line !~ /Mult in  hum ref/)
	{
	    die("Expected to see \"Mult in  hum ref\" but instead saw $line");
	}
	$line = <$fh>;
	chomp $line;
	my @br1_refmult=split(/\s+/, $line);
	
	
	$line = <$fh>;
	if ($line !~ /Covg in indiv/)
	{
	    die("Expected to see \"Covg in indiv\" but instead saw $line");
	}
	$line = <$fh>;
	chomp $line;
	my @br1_cov = split(/\s+/, $line);
	
	<$fh>;
	$line = <$fh>;
	if ($line !~ /Mult in  hum ref/)
	{
	    die("Second time - Expected to see \"Mult in  hum ref\" but instead saw $line");
	}
	$line = <$fh>;
	chomp $line;
	
	my @br2_refmult = split(/\s+/, $line);
	
	$line = <$fh>;
	if ($line !~ /Covg in indiv/)
	{
	    die("Second branch - Expected to see \"Covg in indiv\" but instead saw $line");
	}
	$line = <$fh>;
	chomp $line;
	my @br2_cov = split(/\s+/, $line);
	
	
	my $stat_br1_humref = Statistics::Descriptive::Full->new();
	$stat_br1_humref->add_data(@br1_refmult);
	my $stat_br2_humref = Statistics::Descriptive::Full->new();
	$stat_br2_humref->add_data(@br2_refmult);
	my $stat_br1_cov =  Statistics::Descriptive::Full->new();
	$stat_br1_cov -> add_data(@br1_cov);
	my $stat_br2_cov =  Statistics::Descriptive::Full->new();
	$stat_br2_cov -> add_data(@br2_cov);
	
	my $median_br1_cov = $stat_br1_cov->median();
	my $median_br2_cov = $stat_br2_cov->median();
	my $min_refmult_br1 = $stat_br1_humref->min();
	my $min_refmult_br2 = $stat_br2_humref->min();
	
	
	##determine which is ref allele
	my $which_is_ref;
	if ( ($min_refmult_br1>=1) && ($min_refmult_br2==0) )
	{
	    $which_is_ref=1;
	}
	elsif ( ($min_refmult_br2>=1) && ($min_refmult_br1==0) )
	{
	    $which_is_ref=2;
	}
	elsif ( ($min_refmult_br2>=1) && ($min_refmult_br1>=1) )
	{
	    $which_is_ref="b";#both
	}
	else
	{
	    $which_is_ref="neither"
	}
	my $eof="";
	return  ($eof, $varname, $flank5p, $br1, $br2, $flank3p, $median_br1_cov, $median_br2_cov, $which_is_ref);
	
    }
}



sub get_next_var_from_flank_mapfile
{
    my ($fh) = @_;


    my $line="@";
    while ($line =~ /^\@/)
    {

	if (eof($fh))
	{
	    return ("EOF", 0,0,0,0);
	}
	$line = <$fh>;
    }
    chomp $line;
    my @sp = split(/\t/, $line);
    if (scalar(@sp)<10)
    {
	die("Unexpected num <10 of fields on $line");
    }
    my $name ;
    if ($sp[0] =~ /(var_\d+)/)
    {
	$name = $prefix."_".$1;
    }
    else
    {
	die("Unexpected format of first field in $line");
    }
    my $strand = $sp[1];
    my $chr = $sp[2];
    my $coord = $sp[3];
    my $eof="";
    return ($eof, $name, $strand, $chr, $coord);
}



## return info about how the alignment of the two branches goes.
## if it turns out that an alignment of the rev comp of the branches looks good,
## then return that info as a binary datum, so can flag the variant as possible inversion
sub get_next_alignment_from_procfile
{

    #we almost do not need $dir_of_alignment_of_5pflank, except at the end when we get coordinates of SNPs in a set pf phased SNPS
    #similarly for $which_is_ref. Must handle gracefully if this has value b (both) or neither. In both these cases we will ignore the variant
    #and we only call this function to make sure we be at the right point in the file to go to the NEXT one, next time we call this. ie just read past this var.
    my ($fh, $dir_of_alignment_of_5pflank, $which_is_ref) = @_;

    my $line="";
    while ($line !~ /START NEW VAR/)
    {
	if (eof($fh))
	{
	    return ("EOF", 0,0,0,0,0,0,0,0,0);
	}
	$line = <$fh>;
    }

    $line = <$fh>;
    my $name;
    if ($line =~ /^(\S+)/)
    {
	$name = $1;
    }
    else
    {
	die("Expected name on $line");
    }
    <$fh>;

    my $fw_br1="";
    my $rev_br1="";
    my $fw_middle="";
    my $rev_middle="";
    my $fw_br2="";
    my $rev_br2="";
    my $fw_num_snps=0;
    my $rev_num_snps=0;
    my $fw_num_indels=0;
    my $rev_num_indels=0;

    $line = <$fh>;
    if ($line !~ /FORWARD ALIGNMENT/)
    {
	die("format issue - this line $line should have said FORWARD ALIGNMENT. Current var is $name\n");
    }
    $fw_br1=<$fh>;
    chomp $fw_br1;
    $fw_middle=<$fh>;
    chomp $fw_middle;
    $fw_br2=<$fh>;
    chomp $fw_br2;
    $fw_br1 =~ s/Br1://;
    $fw_br2 =~ s/Br2://;
    if ($fw_middle =~ /^\s{4}(.+)/)
    {
	$fw_middle = $1;
    }
    else
    {
	die("bad format of fw middle $fw_middle. which does not have 4 space sat the start - put there so it all aligns visually when printed, but nothing to do with the alignment");
    }

   $line = <$fh>;
    chomp $line;

    if ($line =~ /(\d+)\s+(\d+)/)
    {
	$fw_num_snps = $1;
	$fw_num_indels=$2;
    }
    else
    {
	die("Bad format of $line - expect num snps and indels");
    }
    <$fh>;


    my $there_is_rev_alignment=0;
    $line = <$fh>;
    if ($line =~ /NO REVERSE ALIGNMENT/)
    {
	
    }
    elsif ($line =~ /REVERSE ALIGNMENT/)
    {
	$there_is_rev_alignment=1;
	$rev_br1=<$fh>;
	chomp $rev_br1;
	$rev_middle=<$fh>;
	chomp $rev_middle;
	$rev_br2=<$fh>;
	chomp $rev_br2;
	$rev_br1 =~ s/Br1://;
	$rev_br2 =~ s/Br2://;
	if ($rev_middle =~ /^\s{4}(.+)/)
	{
	    $rev_middle = $1;
	}
	else
	{
	    die("bad format of rev middle $rev_middle. which does not have 4 space sat the start - put there so it all aligns visually when printed, but nothing to do with the alignment");
	}
	
	my $line = <$fh>;
	chomp $line;

	if ($line =~ /(\d+)\s+(\d+)/)
	{
	    $rev_num_snps = $1;
	    $rev_num_indels=$2;
	}
	else
	{
	    die("Bad format of $line - expect num snps and indels");
	}
    }


    if ( ($which_is_ref eq "b") || ($which_is_ref eq "neither") )
    {
	## No need for the rest of this function, we will ignore this variant.
	return ("", 0,0,0,0,0,0,0,0,0,0,0,0,0);
    }





    ##  if in doubt, say there are 0 bases of agreement at start/end, so that no clever trimming is done

    ## 0 =no, 1=inversion (pure), 2=deletion+inversion or insertion + inversion
    my $possible_inversion=0;
    if ($there_is_rev_alignment==1)
    {
	##then we have an inversion + maybe deletion
	if ($rev_middle =~ /^\s+[\|]+$/)
	{
	    $possible_inversion=2;
	}
	elsif ($rev_middle =~ /^[\|]+\s+$/)
	{
	    $possible_inversion=2;
	}
	elsif ($rev_middle =~ /^[\|]+$/) ##perfect inversion
	{
	    $possible_inversion=1;
	}
	else
	{
	    my $max_consec_matches_fw = get_max_consecutive_matches($fw_middle);
	    my $max_consec_matches_rev = get_max_consecutive_matches($rev_middle);

	    if ($max_consec_matches_fw<$max_consec_matches_rev )
	    {
		    $possible_inversion=2;
	    }


	}
    }


    my $clean_indel=0;
    if ($fw_middle =~ /^\s+[\|]+$/)
    {
	$clean_indel=1;
    }
    elsif ($fw_middle =~ /^[\|]+\s+$/)
    {
	$clean_indel=1;
    }




    my $num_bases_agreement_at_start=0;
    my $num_bases_agreement_at_end=0;
    my $align_direction;
    my $num_snps;
    my $num_indels;

    $align_direction="+";
    if ( (length($fw_br1)==length($fw_br2) ) 
	 &&
	 ($fw_middle =~ /^([\|]+)/)
	)
    {
	if ($fw_middle =~ /^([\|]+)/)
	{
	    $num_bases_agreement_at_start=length($1);
	}
    }
    if ( (length($fw_br1)==length($fw_br2) ) 
	 &&
	 ($fw_middle =~ /([\|]+)$/)
	)
    {
	if ($fw_middle =~ /([\|]+)$/)
	{
	    $num_bases_agreement_at_end=length($1);
	}
    }
    $num_snps = $fw_num_snps;
    $num_indels = $fw_num_indels;
    


    #if this is a bunch of phased SNPs, then get their positions
    my @coords_of_snps=();
    my @indices_of_snps=();
    my @alleles_of_snps=();


    my @indices_of_indels=();
    my @indices_of_indel_ends=();
    my @coords_of_indels=();
    my @alleles_of_indels=();
    my @indel_flag_needs_extra_base=(); #will be 1 if that indel needs an extra base on the front. Should only be possible for ONE of th indels

    my $warning_need_to_add_flank_base=0;
    
    if ( ($num_snps>0) && ($num_indels==0) ) ## this can only happen if the two branches are the same length
    {

	if ($dir_of_alignment_of_5pflank==0)
	{
	    ### first get the SNP coords
	    my $j;
	    my @sp = split(//, $fw_middle);
	    for ($j=0; $j<length($fw_middle); $j++)
	    {
		if ($sp[$j] eq "\*")
		{
		    push @coords_of_snps, $j; 
		}
	    }
	    ##then get the SNP alleles, always in form REF_allele,ALT_allele. Last argument is ignored as in fw direction
	    get_snp_alleles(\@coords_of_snps, \@alleles_of_snps, $fw_br1, $fw_middle, $fw_br2, "fw",$which_is_ref, -1);

	}
	elsif ($dir_of_alignment_of_5pflank==16)
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
	    my @sp = split(//, $fw_middle);
	    my $leng = scalar(@sp);
	    ##first find first base of disagreement starting from the end
	    my $fwd_pos_of_last_snp=$leng-1;

            for ($j=$leng-1 ; $j>=0; $j--)
            {
                if ($sp[$j] ne "\|")
                {
                    $fwd_pos_of_last_snp=$j;
		    last;
		}
            }
	    for ($j=0 ;$j<=$fwd_pos_of_last_snp; $j++)
	    {
		if ($sp[$fwd_pos_of_last_snp-$j] eq "\*")
		{
		    push @coords_of_snps, $j;
		}
	    }

	    ##then get the SNP alleles, always in form REF_allele,ALT_allele.  Note we DELIBERATELY pass in fw_br, fw_middle etc AND "rev".
	    get_snp_alleles(\@coords_of_snps, \@alleles_of_snps, $fw_br1, $fw_middle, $fw_br2, "rev",$which_is_ref, $fwd_pos_of_last_snp);

	    
	}
    }


    elsif ( ($possible_inversion==0)  && ($num_snps>0) && ($num_indels>0) )
    {
	
	if ($dir_of_alignment_of_5pflank==0)
	{
	    ### first get the SNP coords and indices
	    my $j;
	    my @sp = split(//, $fw_middle);
	    my @split_br1 = split(//, $fw_br1);
	    my @split_br2 = split(//, $fw_br2);

	    my $underscore_on_ref_allele_so_far=0;

	    for ($j=0; $j<length($fw_middle); $j++)
	    {


		if (($which_is_ref==1)&&($split_br1[$j] eq "_") )
		{
		    $underscore_on_ref_allele_so_far++;
		}
		elsif (($which_is_ref==2)&&($split_br2[$j] eq "_") )
		{
		    $underscore_on_ref_allele_so_far++;
		}

		if ($sp[$j] eq "\*")
		{
		    #check base before/after is not an indel
		    if  (   ( ($j>0) && ($j<length($fw_middle)-1) && ($sp[$j-1] ne " ") && ($sp[$j+1] ne " ") )
			    ||
			    ( ($j==0) && ($sp[$j+1] ne " ") )
			    || 
			    ( ($j==length($fw_middle)-1) && ($sp[$j-1] ne " ") )
			)
		    {
			push @indices_of_snps, $j;
			push @coords_of_snps, $j - $underscore_on_ref_allele_so_far; 
		    }
		}


	    }
	    ##then get the SNP alleles, always in form REF_allele,ALT_allele. Last argument is ignored as in fw direction
	    get_snp_alleles(\@indices_of_snps, \@alleles_of_snps, $fw_br1, $fw_middle, $fw_br2, "fw",$which_is_ref, -1);


	    ## now get the indel INDICES of where the indels are in these strings, and hence will get the COORDINATES
	    ##  eg ig we have AAA_T__C on the ref allele, then the index of C is 7, and the coordinate is 5
	    ## will get the coord just before the indel in the forward direction along the ref, and for ref/alt alleles give that prior base also.
	    ## look for gaps in ||| ||||| that are preceded by two || and succeeded also by two.


	    $underscore_on_ref_allele_so_far=0;

	    for ($j=0; $j<length($fw_middle); $j++)
            {

                if ($sp[$j] eq " ")
                {
		    if (
			( ($j==0) )
			||
			( ($j==1) && ($sp[$j-1] ne " ")  )
			||
			( ($j>1) && ($sp[$j-1] ne " ") && ($sp[$j-2] ne " ") )
			)
		    {

			my $indel_start_coord=$j - $underscore_on_ref_allele_so_far - 1;##base before indel
			my $indel_start_index = $j-1 ;## base before first base of space 

			##we have found an indel. Now find the end: - keep going to you see two |, or reach end
			while (
			    (     ! ( ($j<scalar(@sp)-2) && ($sp[$j+1]  eq "\|") && ($sp[$j+2] eq "\|") )
				  &&
				  !( ($j==scalar(@sp)-2) && ($sp[$j+1] eq "\|") )
				  &&
				  !( ($j==scalar(@sp)-1) )
			    )
			    &&
			    ($j<scalar(@sp))
			    )
			{

			    if (($which_is_ref==1)&&($split_br1[$j] eq "_") )
			    {
				$underscore_on_ref_allele_so_far++;
			    }
			    elsif (($which_is_ref==2)&&($split_br2[$j] eq "_") )
			    {
				$underscore_on_ref_allele_so_far++;
			    }

			    $j++;
			}
			
			push @indices_of_indels, $indel_start_index;
			push @coords_of_indels, $indel_start_coord;
			push @indices_of_indel_ends, $j; 

		    }
		}

		if (($which_is_ref==1)&&($split_br1[$j] eq "_") )
		{
		    $underscore_on_ref_allele_so_far++;
		}
                elsif (($which_is_ref==2)&&($split_br2[$j] eq "_") )
                {
                    $underscore_on_ref_allele_so_far++;
                }

	    }

	    ## now get the indel alleles
	    get_indel_alleles(\@indices_of_indels, \@indices_of_indel_ends, \@alleles_of_indels, \@indel_flag_needs_extra_base,
			      $fw_br1, $fw_middle, $fw_br2, "fw",$which_is_ref, -1);



	}
	
	## get forward INDEX of base just after last space in the indel -  ie when we reverse, will be the base just befroe the indel
	## also, get the coord of that same base.
	elsif ($dir_of_alignment_of_5pflank==16)
	{
	    ## e.g.
	    ##Br1:CGCCGTTGTTGAGTGTTCTATGG__TTGTCGTTTATTGAGCACAACTACAGCATTT
            ##    *|||*||||||||||||||*|||  |||||||||||||||||||||||||||||||
            ##Br2:TGCCCTTGTTGAGTGTTCTTTGGAATTGTCGTTTATTGAGCACAACTACAGCATTT
	    ##                             ^find this base


	    my $j;
	    my @sp = split(//, $fw_middle);
	    my @split_br1 = split(//, $fw_br1);
	    my @split_br2 = split(//, $fw_br2);
	    my $underscore_on_ref_allele_so_far=0;

	    my $leng = scalar(@sp);
	    ##first find first base of disagreement starting from the end
	    my $fwd_pos_of_last_non_match=$leng-1;

            for ($j=$leng-1 ; $j>=0; $j--)
            {
                if ($sp[$j] ne "\|")
                {
                    $fwd_pos_of_last_non_match=$j;
		    last;
		}
            }
	    for ($j=0 ;$j<=$fwd_pos_of_last_non_match; $j++)
	    {

		if ($sp[$fwd_pos_of_last_non_match-$j] eq "\*")
		{
		    push @indices_of_snps, $j;
		    push @coords_of_snps, $j -$underscore_on_ref_allele_so_far;
		}

		if (($which_is_ref==1)&&($split_br1[$fwd_pos_of_last_non_match-$j] eq "_") )
		{
		    $underscore_on_ref_allele_so_far++;
		}
		elsif (($which_is_ref==2)&&($split_br2[$fwd_pos_of_last_non_match-$j] eq "_") )
		{
		    $underscore_on_ref_allele_so_far++;
		}

	    }

	    ##then get the SNP alleles, always in form REF_allele,ALT_allele.  Note we DELIBERATELY pass in fw_br, fw_middle etc AND "rev".
	    get_snp_alleles(\@indices_of_snps, \@alleles_of_snps, $fw_br1, $fw_middle, $fw_br2, "rev",$which_is_ref, $fwd_pos_of_last_non_match);
	 

	    ##now do indels.

	    $underscore_on_ref_allele_so_far=0;

	    for ($j=0; $j<= $fwd_pos_of_last_non_match ; $j++)
            {

                if ($sp[$fwd_pos_of_last_non_match - $j] eq " ")
                {
		    if (##there are two | on the RHS
			( ($j==0) )
			||
			( ($j==1) && ($sp[$fwd_pos_of_last_non_match-$j+1] eq "\|")  )
			||
			( ($j>1) && ($sp[$fwd_pos_of_last_non_match-$j+1] eq "\|") && ($sp[$fwd_pos_of_last_non_match-$j+2] eq "\|") )
			)

		    {

			#                 v one after fwd pos of last non-match - in fact might actually be AFTER the whole string.
			#  A___CCCCCCCCCC_C
			#      ^- 10 is the index_start_coord for this deletion ___
			##indel_start_coord is the absolute size of the difference between the fwd_pos_of_last_non_match, and the base after (going forwards) the start
			my $indel_start_coord= $j - $underscore_on_ref_allele_so_far-1;##base before indel
			my $indel_start_index = $fwd_pos_of_last_non_match -$j+1 ;## base before first base of space 

                        ##we have found an indel. Now find the end: - keep going to you see two |, or reach end                                                                                 
                        while (    !(
					( ($j<= $fwd_pos_of_last_non_match-2) && ($sp[$fwd_pos_of_last_non_match -$j-1]  eq "\|") && ($sp[$fwd_pos_of_last_non_match-$j-2] eq "\|") )
					||
					( ($j==$fwd_pos_of_last_non_match-1) && ($sp[$fwd_pos_of_last_non_match-$j-1] eq "\|") )
					||
					( ($j==$fwd_pos_of_last_non_match) )
				   )
				   &&
                                   ($j<=$fwd_pos_of_last_non_match)
                            )
                        {

			    if (($which_is_ref==1)&&($split_br1[$fwd_pos_of_last_non_match-$j] eq "_") )
			    {
				$underscore_on_ref_allele_so_far++;
			    }
			    elsif (($which_is_ref==2)&&($split_br2[$fwd_pos_of_last_non_match-$j] eq "_") )
			    {
				$underscore_on_ref_allele_so_far++;
			    }
			    

			    $j++;
			}
			
			push @indices_of_indels, $indel_start_index;
			push @coords_of_indels, $indel_start_coord;
			push @indices_of_indel_ends, $fwd_pos_of_last_non_match -$j;

		    }
		}

		if (($which_is_ref==1)&&($split_br1[$fwd_pos_of_last_non_match - $j] eq "_") )
		{
		    $underscore_on_ref_allele_so_far++;
		}
                elsif (($which_is_ref==2)&&($split_br2[$fwd_pos_of_last_non_match - $j] eq "_") )
                {
                    $underscore_on_ref_allele_so_far++;
                }

	    }

	    ## now get the indel alleles
            get_indel_alleles(\@indices_of_indels, \@indices_of_indel_ends, \@alleles_of_indels, \@indel_flag_needs_extra_base,
			      $fw_br1, $fw_middle, $fw_br2, "rev",$which_is_ref, $fwd_pos_of_last_non_match);



   
	}

	




    }


    
    my $eof="";
    

    return  ($eof, $name, $num_bases_agreement_at_start, $num_bases_agreement_at_end, $align_direction,
	     $num_snps, $num_indels,$fw_br1, $fw_middle, $fw_br2, $possible_inversion, 
	     $clean_indel, \@coords_of_snps, \@alleles_of_snps,
	     \@coords_of_indels, \@alleles_of_indels, \@indel_flag_needs_extra_base);

    

}

## string of form |||||| |||| | | || |          || || | || 
sub get_max_consecutive_matches
{
    my ($str) = @_;

    my $max=0;

    my @sp = split(/\s+/, $str);
    foreach my $word (@sp)
    {
	if (length($word)>$max)
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


    my ($aref_snp_coor, $aref_snp_allel, $b1, $middle, $b2, $dir, $which_is_ref, $index_of_last_non_match) = @_;


    ##double check
    if (  (length($b1) != length($b2) )
	  ||
	  (length($b1) != length($middle) )
	)
    {
	die("This alignment has the 3 lines of different lengths.\n$b1\n$middle\n$b2\n");
    }
    
    if ($dir eq "fw")
    {
	my @sp1 = split(//, $b1);
	my @sp_middle = split(//, $middle);
	my @sp2 = split(//, $b2);
	my $i;
	for ($i=0; $i<scalar(@$aref_snp_coor); $i++)
	{
	    if ($which_is_ref==1)
	    {
		push @$aref_snp_allel, substr($b1, $aref_snp_coor->[$i], 1)."_".substr($b2, $aref_snp_coor->[$i], 1);
	    }
	    elsif ($which_is_ref==2)
	    {
		push @$aref_snp_allel, substr($b2, $aref_snp_coor->[$i], 1)."_".substr($b1, $aref_snp_coor->[$i], 1);
	    }
	}
    }
    else##dir is rev
    {
	my @sp1 = split(//, $b1);
	my @sp_middle = split(//, $middle);
	my @sp2 = split(//, $b2);
	my $i;
	for ($i=0; $i<scalar(@$aref_snp_coor) ; $i++)
	{
	    if ($which_is_ref==1)
	    {
		push @$aref_snp_allel, 
		    rev_comp($sp1[$index_of_last_non_match - ($aref_snp_coor->[$i])] )
		    ."_".
		    rev_comp($sp2[$index_of_last_non_match - ($aref_snp_coor->[$i])] );
	    }
	    elsif ($which_is_ref==2)
	    {
		push @$aref_snp_allel, 
		    rev_comp($sp2[$index_of_last_non_match - ($aref_snp_coor->[$i])] )
		    ."_".
		    rev_comp($sp1[$index_of_last_non_match - ($aref_snp_coor->[$i])] );
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


    my ($aref_indel_indices, $aref_indel_end_indices, $aref_indel_allel, $aref_indel_need_base_adding, $b1, $middle, $b2, $dir, $which_is_ref, $index_of_last_non_match) = @_;


    ##double check
    if (  (length($b1) != length($b2) )
	  ||
	  (length($b1) != length($middle) )
	)
    {
	die("This alignment has the 3 lines of different lengths.\n$b1\n$middle\n$b2\n");
    }
    

    if ($dir eq "fw")
    {

	my @sp1 = split(//, $b1);
	my @sp_middle = split(//, $middle);
	my @sp2 = split(//, $b2);
	my $i;
	for ($i=0; $i<scalar(@$aref_indel_indices); $i++)
	{
	    my $flag_needs_extra_base=0;
	    my $ref="";
	    my $alt="";
	    my $k  =$aref_indel_indices->[$i];##base before indel's first space

	    if ($k<-1)
	    {
		die("k should not be $k\n");
	    }

	    if ($k==-1 )#removed constraint that index of last nonmatch==0
	    {
		$flag_needs_extra_base=1;
		$k++;
	    }

	    if ($which_is_ref==1)
	    {
		$ref = $ref.$sp1[$k];
		$alt = $alt.$sp2[$k];
	    }
	    else
	    {
		$ref = $ref.$sp2[$k];
		$alt = $alt.$sp1[$k];		
	    }
	    $k++;
	    
	    while ($k<= $aref_indel_end_indices->[$i]  )
	    {
		if ($which_is_ref==1)
		{
		    $ref = $ref.$sp1[$k];
		    $alt = $alt.$sp2[$k];
		}
		else
		{
		    $ref = $ref.$sp2[$k];
		    $alt = $alt.$sp1[$k];
		}
		
		$k++;
	    }
	    $ref =~ s/_//g;
	    $alt =~ s/_//g;
	    push @$aref_indel_allel, $ref."_".$alt;
	    push @$aref_indel_need_base_adding, $flag_needs_extra_base;
	}
    }
    else##dir is rev
    {

	my @sp1 = split(//, $b1);
	my @sp_middle = split(//, $middle);
	my @sp2 = split(//, $b2);
	my $i;
	for ($i=0; $i<scalar(@$aref_indel_indices) ; $i++)
	{
	    my $flag_needs_extra_base=0;
	    my $ref="";
	    my $alt="";
	    my $k  =$aref_indel_indices->[$i];

	    
	    if ($which_is_ref==1)
	    {

		if ($k>scalar(@sp1)-1 )
		{
		    $flag_needs_extra_base=1;
		    $k--;
		}

		$ref = $ref.$sp1[$k];
		$alt = $alt.$sp2[$k];
	    }
	    else
	    {

		if ($k>scalar(@sp2)-1)
		{
		    $flag_needs_extra_base=1;
		    $k--;
		}

		$ref = $ref.$sp2[$k];
		$alt = $alt.$sp1[$k];		


	    }

	    $k--;
	    
	    while ( $k>=$aref_indel_end_indices->[$i])
	    {
		if ($which_is_ref==1)
		{
		    $ref = $ref.$sp1[$k];
		    $alt = $alt.$sp2[$k];
		}
		else
		{
		    $ref = $ref.$sp2[$k];
		    $alt = $alt.$sp1[$k];
		}
		
		$k--;
	    }
	    $ref =~ s/_//g;
	    $alt =~ s/_//g;
	    $ref = reverse($ref);
	    $alt = reverse($alt);
	    $ref = rev_comp($ref);
	    $alt = rev_comp($alt);

	    push @$aref_indel_allel, $ref."_".$alt;
	    push @$aref_indel_need_base_adding, $flag_needs_extra_base;
	    
	}

    }


}




sub combine_all_filters
{
    my ($href_var_name_to_covg_and_branch_filter, 
	$href_var_name_to_flank_mq_filter, 
	$href_var_name_to_combined_filtering_result) = @_;


    foreach my $key (keys %$href_var_name_to_flank_mq_filter)
    {
	if ( ($href_var_name_to_flank_mq_filter->{$key} eq "PASS")
	     &&
	     (exists $href_var_name_to_covg_and_branch_filter->{$key})
	     &&
	     ($href_var_name_to_covg_and_branch_filter->{$key} eq "PASS")
	    )
	{
		$href_var_name_to_combined_filtering_result->{$key}="PASS";
	}
	elsif (!exists $href_var_name_to_covg_and_branch_filter->{$key})
	{
	    die("$key is not in the covg and branch filter hash\n");
	}
	elsif ($href_var_name_to_covg_and_branch_filter->{$key} ne "PASS")
	{
	    print "$key has failed the covg & branch filter  \n";
	}
	elsif ($href_var_name_to_flank_mq_filter->{$key} ne "PASS")
	{
	    print "$key fails mapqual filter\n";
	}
    }
    
    
}



sub get_list_vars_with_cut_flanks
{
    my ($file, $href) = @_;

    open(FILE, $file)||die();

    while (<FILE>)
    {
	my $line = $_;
	if ( ($line =~ /var_\d+_5p_flank/) && ($line =~ /cut_at_1000/) )
	{
	    if ($line =~ /(var_\d+)_5p_flank/)
	    {
		my $name = $prefix."_".$1;
		$href->{$name}=1;
	    }
	    else
	    {
		die("programming error on $line");
	    }
	}
    }
    close(FILE);
}
