package RunCallsLib;

use strict;
use warnings;
use File::Basename;
use List::Util qw(min max sum);


my $callingscript_dir;
my $analyse_variants_dir;
my $cortex_dir;
my $isaac_bioinf_dir;

BEGIN
{
	use FindBin;
	$callingscript_dir = $FindBin::Bin;
	$cortex_dir = $callingscript_dir;
	$cortex_dir =~ s/scripts\/calling//;
	$analyse_variants_dir = $cortex_dir."/scripts/analyse_variants/";
	$isaac_bioinf_dir = $analyse_variants_dir."bioinf-perl/";
}

my $make_union = $analyse_variants_dir."make_union_varset_at_single_kmer.pl";
if (!(-e $make_union))
{
    die("Cannot find make_union_varset script: $make_union\n");
}


# Use current directory to find modules
use base 'Exporter';
our @EXPORT = qw(make_a_union_callset_for_each_kmer genotype_union make_sure_dirs_exist_and_create_if_necessary get_max_cleaning_for_given_sample get_min_cleaning_for_given_sample apply_pop_classifier delete_and_add_header combine_and_filter_dups_overlaps_and_ref_mismatches_from_vcf);

sub combine_and_filter_dups_overlaps_and_ref_mismatches_from_vcf
{
    my ($outfilename, $tmpdir, $aref, $vcftools_dir, $workflow, $num_samples) = @_;
    
    my $tmp = $tmpdir.basename($outfilename).".temp_cated_vcf";
    my $catout_fh;
    open ($catout_fh, ">".$tmp)||die("Unable to open temp file $tmp");
    print_header_from_vcf($aref->[0], $catout_fh);
	
    my $i;
    ## concatenate all the vcfs we have so far
    for ($i=0; $i<scalar @$aref; $i++)
    {
	print_vcf_contents($aref->[$i], $catout_fh);
    }
    close($catout_fh);

    
    ## Now cleanup - remove dups etc
    my $tmp_midpoint_fails = $tmpdir.basename($outfilename).".temp_vcf_of_fails_midway_through_process";
    delete_and_add_header($tmp, $tmp_midpoint_fails);
    my $tmp_midpoint_passes = $tmpdir.basename($outfilename).".temp_vcf_of_passes_midway_through_process";
    delete_and_add_header($tmp, $tmp_midpoint_passes);
    my $tmp_nearly_there = $tmpdir.basename($outfilename).".next_temp_vcf";
    if (-e $tmp_nearly_there)
    {
	my $c = "rm $tmp_nearly_there";
	qx{$c};
    }
    #delete_and_add_header($tmp, $tmp_nearly_there);
    my $tmp_final_unsorted = $tmpdir.basename($outfilename).".final_unsorted";

    ##separate passes and fails
    my $cmd1 = "($vcftools_dir/perl/vcf-sort $tmp | " .
               $isaac_bioinf_dir."vcf_scripts/vcf_remove_dupes.pl --take_first  | " .
               $isaac_bioinf_dir."vcf_scripts/vcf_combine_alleles.pl --pass  | grep PASS | grep -v \"\#\"  >> $tmp_midpoint_passes) 2>&1";
    my $cmd2 = "($vcftools_dir/perl/vcf-sort $tmp | " .
               $isaac_bioinf_dir."vcf_scripts/vcf_remove_dupes.pl --take_first  | " .
               $isaac_bioinf_dir."vcf_scripts/vcf_combine_alleles.pl --pass  | grep -v PASS | grep -v \"\#\" >> $tmp_midpoint_fails) 2>&1";
    
    my $cmd3;
    if ( ($workflow eq "independent") && ($num_samples>1) )
    {
	    $cmd3 =  "( ($vcftools_dir/perl/vcf-sort $tmp_midpoint_passes | ".
		$isaac_bioinf_dir."vcf_scripts/vcf_remove_overlaps.pl  --filter_txt OVERLAPPING_SITE  | $vcftools_dir/perl/vcf-sort | ". 
		$isaac_bioinf_dir."vcf_scripts/vcf_revert_consistent_overlaps.pl) >> $tmp_nearly_there)  2>&1 ";
    }
    else #no overlap removal for joint, or for independent with one sample only
    {
	$cmd3 =  "($vcftools_dir/perl/vcf-sort $tmp_midpoint_passes  >> $tmp_nearly_there)  2>&1 ";
    }
    my $cmd4 = "(cat $tmp_nearly_there | head -100 | grep ^#; cat $tmp_nearly_there | grep -v \"\#\" ; cat $tmp_midpoint_fails | grep -v \"\#\" ;) > $tmp_final_unsorted ";
    my $cmd5  = $vcftools_dir."/perl/vcf-sort $tmp_final_unsorted >  $outfilename ";
    
    print "$cmd1\n";
    my $ret1 = qx{$cmd1};
    print "$ret1\n";
    print "$cmd2\n";
    my $ret2 = qx{$cmd2};
    print "$ret2\n";
    print "$cmd3\n";
    my $ret3 = qx{$cmd3};
    print "$ret3\n";
    print "$cmd4\n";
    my $ret4 = qx{$cmd4};
    print "$ret4\n";
    print "$cmd5\n";
    my $ret5 = qx{$cmd5};
    print "$ret5\n";
    
}

sub delete_and_add_header
{
    my ($file_with_header, $file_missing_header) = @_;
    my $fh;
    if (-e $file_missing_header)
    {
	print "Delete pre-existing file $file_missing_header ready to make clean start\n";
	my $cmd = "rm $file_missing_header";
	qx{$cmd};
    }
    open($fh, ">".$file_missing_header)||die("Cannot open $file_missing_header - out of disk space? permissions issue?\n");
    print_header_from_vcf($file_with_header, $fh);
    close($fh);
}


##for now assume there is a ref - will need to fix this
sub apply_pop_classifier
{
    my ($classifier, $make_covg_script, $make_table_script,
	$genotyping_output, $genotyping_output_log,
	$num_cols, $use_reference, $kmer, $genome_size, $ploidy) = @_;

    my $is_ref;
    my $refcol=-1;
    if ($use_reference eq "Absent")
    {
	$is_ref=0;
    }
    else
    {
	$is_ref=1;
	$refcol=0;
    }

    ## first - make covg file
    my $covgfile = $genotyping_output.".covg_for_classifier";

    if (!(-e $covgfile))
    {
	my $cmd1 = "perl $make_covg_script $genotyping_output $num_cols $refcol";
	print "$cmd1\n";
	my $ret1 = qx{$cmd1};
	print "$ret1\n";
	if (!(-e $covgfile))
	{
	    die("Unable to run $make_covg_script successfully on $genotyping_output\n");
	}
    }
    else
    {
	print "$covgfile already exists, no need to regenerate\n";
    }
    
    ## second make a table of covgs etc
    my $table = $genotyping_output_log.".table";
    my $cmd2 = "perl $make_table_script $genotyping_output_log > $table 2>&1";
    print "$cmd2\n";
    my $ret2 = qx{$cmd2};
    print "$ret2\n";

    ## count the calls
    my $num_calls = count_calls_in_callfile($genotyping_output);
    
    
    my $out = $covgfile.".classified.split_start_1";
    if ($num_calls == 0) { system "touch $out" }
#    my $log = $out.".log";

    if (!(-e $out))
    {
	my $cmd3 = "cat $classifier  | R --vanilla --args 1 $num_calls $covgfile $num_calls $num_cols $is_ref $table $genome_size $kmer $ploidy $out ";
	print "$cmd3\n";
	my $ret3 = qx{$cmd3};
	print "$ret3\n";

	if (!(-e $out))
	{
	    die("Ran classifier but could not find output file $out\n");
	}
    }

    return $out;
}


sub count_calls_in_callfile
{
    my ($file) = @_;
    
    my $cmd = "cat $file | grep 5p_flank | tail -n 1";
    print "Count how many calls in $file\n";
    print "$cmd\n";
    my $ret = qx{$cmd};
    print "$ret\n";
    my $num=-1;
    if ($ret =~ /var_(\d+)_5p_flank/)
    {
	$num = $1;
    }
    else
    {
	#die("Cannot parse $ret to find number of calls, in count_calls_in_callfile\n");
	$num = 0;
    }
    if ($num==-1)
    {
	die("num is -1 in count_calls_in_callfile");
    }
    return $num;
}


sub make_sure_dirs_exist_and_create_if_necessary
{
    my ($aref) = @_;
    foreach my $dir (@$aref)
    {
	if (! (-d $dir))
	{
	    my $mkdir_cmd = "mkdir -p $dir";
      print "$mkdir_cmd\n";
	    my $mkdir_ret = qx{$mkdir_cmd};
      print "$mkdir_ret";
	}
    }

}

sub get_max_cleaning_for_given_sample
{
    my ($href_sample_to_cleaned_bin, $sample, $kmer) = @_;
    my $min =999999999;
    my $max=0;
    
    foreach my $c (keys %{$href_sample_to_cleaned_bin->{$sample}{$kmer}})
    {
	if ($c>$max)
	{
	    $max=$c;
	}
    }	
    return $max;
}
    


sub get_min_cleaning_for_given_sample
{
    my ($href_sample_to_cleaned_bin, $sample, $kmer) = @_;
    my $min =999999999;
    
    foreach my $c (keys %{$href_sample_to_cleaned_bin->{$sample}{$kmer}})
    {
	if ($c<$min)
	{
	    $min=$c;
	}
    }
    return $min;

   
}


sub genotype_union
{
    my ($ctx_bin, $colour_list, $kmer_size, $mem_height, $mem_width, 
	$ref_colour, $file_to_genotype,$outfile, $logfile, 
	$which_caller, $max_read_len, $expt_type, $genome_size, $gt_assemblies) = @_;

    	    
    if (!(-e $outfile))
	    {
		my $gt_cmd = $ctx_bin." --colour_list $colour_list  --kmer_size $kmer_size --mem_height $mem_height --mem_width $mem_width --ref_colour $ref_colour --gt $file_to_genotype,$outfile,$which_caller --max_read_len $max_read_len  --print_colour_coverages --experiment_type $expt_type --genome_size $genome_size ";
		if ($gt_assemblies eq "yes")
		{
		    printf("Note - since you used --gt_assemblies yes, we assume these are reference genomes, with a very low per-base error rate, of 0.0000001, in order to genotype\n");
		    $gt_cmd = $gt_cmd." --estimated_error_rate 0.000001 ";
		}
		$gt_cmd = $gt_cmd." > $logfile 2>&1";
		print "$gt_cmd\n";
		my $gt_ret = qx{$gt_cmd};
		print "$gt_ret\n";
	    }
	    else
	    {
		print "$outfile already exists so no need to genotype the union of $which_caller calls at kmer $kmer_size\n";
	    }
}



## returns a hashref with kmer-->filename of union file as $href_results
## and a hashref with kmer -> BC/PD -> max read length

sub make_a_union_callset_for_each_kmer
{
    my ($aref_kmers, $aref_samples, $href_sample_to_kmer_to_cleanings_to_binary, ##inputs
	$href_sample_to_kmer_to_cleaning_to_bc_callfile, $href_sample_to_kmer_to_cleaning_to_pd_callfile,
	$href_joint_bc_callfiles, 
	$href_results, $href_max_readlens, $do_bc, $do_pd, $outdir_calls, $outdir_jointcalls, $tmpdir,	##outputs
	$workflow, $use_ref, $num_cleaning_levels,
	$href_kmer_to_bc_stub, $href_kmer_to_pd_stub)
	= @_;

    print "#********* Make a union callset at each k ************* \n";

    foreach my $k (@$aref_kmers)
    {
	## you are either working with the JOINT or INDEPENDENT workflow, but not both. If joint, you only do BC for the moment
	my $union_of_bc_callsets     = $outdir_calls."union_all_".$workflow."_k".$k."_bc_callsets";
	my $union_of_pd_callsets     = $outdir_calls."union_all_".$workflow."_k".$k."_pd_callsets";
	if ($workflow eq "joint" )
	{
	    $union_of_bc_callsets =~ s/$outdir_calls/$outdir_jointcalls/;
	    $union_of_pd_callsets =~ s/$outdir_calls/$outdir_jointcalls/;### right now this is not used anyway
	}

	my $union_of_bc_callsets_log = $union_of_bc_callsets.".log";
	my $union_of_pd_callsets_log = $union_of_pd_callsets.".log";

	if ( ($do_bc eq "yes") && ($workflow eq "independent" ))
	{
	    my $bc_call_index  = $tmpdir."/index_k".$k."_bc_callfiles";
	    $href_kmer_to_bc_stub->{$k} = make_index_and_union_callset_for_specific_k($bc_call_index, $aref_samples, $href_sample_to_kmer_to_cleanings_to_binary, 
										      $href_sample_to_kmer_to_cleaning_to_bc_callfile, $href_sample_to_kmer_to_cleaning_to_pd_callfile,
										      $href_joint_bc_callfiles,
										      $k, $make_union, $union_of_bc_callsets, $union_of_bc_callsets_log, 
										      $href_max_readlens, $href_results, "BC");
	}
	
	if ( ($do_pd eq "yes") && ($workflow eq "independent") )
	{
	    my $pd_call_index  = $tmpdir."/index_k".$k."_pd_callfiles";
 	    $href_kmer_to_pd_stub->{$k} = make_index_and_union_callset_for_specific_k($pd_call_index, $aref_samples, $href_sample_to_kmer_to_cleanings_to_binary, 
										      $href_sample_to_kmer_to_cleaning_to_bc_callfile, $href_sample_to_kmer_to_cleaning_to_pd_callfile,
										      $href_joint_bc_callfiles,
										      $k, $make_union, $union_of_pd_callsets, $union_of_pd_callsets_log, 
										      $href_max_readlens, $href_results, "PD");
	}


	## if calling jointly, might as well genotype then and there, why print out and then regenotype on the same graph.

	#if ( ($call_jointly_noref eq "yes") || ($call_jointly_incref eq "yes") )
	#{	
	#    my $str="inc_ref";
	#    if ($call_jointly_noref eq "yes") 
	#    {
	#	$str = "exc_ref";
	#    }
	#    my $bc_call_index  = $tmpdir."/list_bc_joint_".$str."_callfiles_k".$k;
	#}## end of joint stuff.
	
	
    }## end for each k
    
    
}




##returns the stub used for naming variants
sub make_index_and_union_callset_for_specific_k
{
    my ($index_name, $aref_samples, $href_sample_to_kmer_to_cleanings_to_binary, 
	$href_sample_to_kmer_to_cleaning_to_bc_callfile, $href_sample_to_kmer_to_cleaning_to_pd_callfile,
	$href_joint_bc_callfiles, 
	$k, $make_union, $outfilename, $outfile_log, 
	$href_max_readlens, $href_results, $which_caller) = @_;


    
    open(CALL, ">".$index_name)||die("Cannot open call_index $index_name in RunCallsLib make_a_union_callset_for_specific_k");
    
    foreach my $sam (@$aref_samples)
    {
	foreach my $cleaning (keys %{$href_sample_to_kmer_to_cleanings_to_binary->{$sam}->{$k}})
	{
	    if ($which_caller eq "BC")
	    {
		print CALL $href_sample_to_kmer_to_cleaning_to_bc_callfile->{$sam}->{$k}->{$cleaning};
	    }
	    else
	    {
		print CALL $href_sample_to_kmer_to_cleaning_to_pd_callfile->{$sam}->{$k}->{$cleaning};
	    }
	    print CALL "\t$k\t$cleaning";
	    print CALL "\n";
	}
    }
    
    close(CALL);
    
    my $stub = "UNION_".$which_caller."_k".$k;

    my $cmd = "perl $make_union --kmer $k --index $index_name --varname_stub $stub  --outfile $outfilename > $outfile_log 2>&1";
    if (! -e($outfilename))
    {
	print "$cmd\n";
	my $ret = qx{$cmd};
  print "$ret\n";
    }
    else
    {
	print "Union of $which_caller callfiles for k=$k already exists: $outfilename,  so will not regenerate it\n";
    }
    $href_max_readlens->{$k}->{$which_caller} = get_max_read_len_of_fasta($outfilename)+$k+10;
    $href_results->{$k}->{$which_caller} = $outfilename;


    return $stub;

    #my $cl=0;
    
    #for ($cl=0; $cl < $num_cleaning_levels; $cl++)
    #{		
    #	print BCCALL $joint_callfiles{$k}{$cl};
#	print BCCALL "\t$k\t$cl\n";
 #   }
  #  close(BCCALL);
#	    
#	    my $bc_cmd = "perl $make_union --filelist $bc_call_index --varname_stub UNION_BC --outfile $union_of_bc_callsets > $union_of_bc_callsets_log 2>&1";
#	    if (! -e($union_of_bc_callsets))
#	    {
#		print "$bc_cmd\n";
#		my $bc_ret = qx{$bc_cmd};
#	    }
#	    else
#	    {
#		print "Union of all joint calls, $union_of_bc_callsets, already exists so will not regenerate\n";
#	    }
#	    $href_max_readlens->{$k}->{"BC"} = get_max_read_len_of_fasta($union_of_bc_callsets)+$last_kmer+10;
#	    $href_results->{$k}->{"BC"} = $union_of_bc_callsets;
#
    
}



sub get_max_read_len_of_fasta
{
    my ($fasta) = @_;
    if (!(-e $fasta))
    {
	die("Cannot open nonexistent file $fasta\n");
    }
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

sub print_header_from_vcf
{
    my ($file, $fh) = @_;
    open(FILE, $file)||die("Cannot open $file");
    
    while (<FILE>)
    {
	my $line = $_;
	if ($line =~ /^\#/)
	{
	    print $fh $line;
	}
	else
	{
	    last;
	}
    }
    close(FILE);
    return;
}


sub print_vcf_contents
{
    my ($file, $fh) = @_;
    open(FILE, $file)||die("Cannot open $file");
    
    while (<FILE>)
    {
	my $line = $_;
	if ($line !~ /^\#/)
	{
	    print $fh $line;
	}
    }
    close(FILE);
    return;
}


1;
