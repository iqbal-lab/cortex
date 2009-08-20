#!/usr/bin/perl -w
use strict;
use lib "/homes/marioc/perl-modules/ensembl/modules";
use lib "/homes/marioc/perl-modules/ensembl-variation/modules";
use lib "/homes/marioc/perl-modules/bioperl-live";
use Bio::EnsEMBL::Registry;


sub atoi;



## Takes as input one of the output files of the trusted-path SV caller, for one chromosome + the number/letter of that chromosome


my $file = shift;
my $chrom = shift;
my $dbsnp="SNPChrPosAllele_b129.txt";
my $must_do_filtering = shift; ##"yes" or "no"
my $dir_with_mapping_data = shift; ##for this chromosome only.
my $registry = shift;


##ensembl stuff  for annotating variants further down
Bio::EnsEMBL::Registry->load_all($registry);
my $ensembl_core_db = Bio::EnsEMBL::Registry->get_DBAdaptor('human','core');
my $ensembl_var_db = Bio::EnsEMBL::Registry->get_DBAdaptor('human','variation');
my $sa = $ensembl_core_db->get_SliceAdaptor;
my $va = $ensembl_var_db->get_VariationAdaptor;
my $ra = $ensembl_var_db->get_RepeatFeatureAdaptor;




## If appropriate: Get hash of variants whose 5prime anchors map uniquely
my %fivep_hash=();
my %threep_hash=();
if ($must_do_filtering eq "yes")
{
    if ($dir_with_mapping_data !~ /\/$/)
    {
	$dir_with_mapping_data = $dir_with_mapping_data."/";
    }

    my $get_filenames_cmd = "ls ".$dir_with_mapping_data."\*.out";
    my $files = qx{$get_filenames_cmd};
    my @files = split(/\n/,$files);


    foreach my $f (@files)
    {
	
	open(FILE,$f)||die("cannot open $f");
	
	while(<FILE>)
	{
	    my $line =$_;
	    chomp $line;
	    
	    ## Example line you're looking for
	    ## vulgar: chromosome_1_VARIATION_100_3p_flank_length_208 0 208 + 1 938315 938523 + 1040 M 208 208
	    ## Will be followed by a line of this form
	    ## position 1:938315-938523 orientation:+ pid:100.00 length:208 start:0 end:208
	    
	    my $varid;
	    my $which_flank;
	    my $query_length;
	    my $pid;
	    my $alignment_length;
	    
	    if ($line =~ /^vulgar: (\S+)_([35]p_flank)_length_(\d+)/) 
	    {
		$varid=$1;
		$which_flank=$2;
		$query_length=$3;
		
		$line=<FILE>;
		chomp $line;
		if ($line =~ /pid:(\d+).*length:(\d+) start/)
		{
		    $pid=$1;
		    $alignment_length=$2;
		}
		else
		{
		    die("Problem with this like $line which should follow a vulgar line and be the pid line\n");
		}
		
		
		if (($pid==100) && ($query_length==$alignment_length))
		{
		    ## We have found a perfect alignment. If it is first, add it to the appropriate hash, with hash_value 1.
		    ## If there is already a hash_value for it, this shows the match is not unique, so set the value to zero.
		    
		    if ($which_flank eq "3p_flank")
		    {
			if (exists $threep_hash{$varid})
			{
			    $threep_hash{$varid}=0;
			}
			else
			{
			    $threep_hash{$varid}=1;
			}
		    }
		    elsif ($which_flank eq "5p_flank")
		    {
			if (exists $fivep_hash{$varid})
			{
			    $fivep_hash{$varid}=0;
			}
			else
			{
			    $fivep_hash{$varid}=1;
			}
		    }
		    else
		    {
			die("Flank is not 3p or 5p - $which_flank");
		    }
		}
		
		
	    }
	}
	
    }

}
elsif ($must_do_filtering eq "no")
{
    die("Hey! You should be filtering");
}
else
{
    die("3rd argument must be yes/no, and you have entered $must_do_filtering");
}





open(FILE, $file)||die("Cannot open $file");
open(DBSNP, $dbsnp)||die("Cannot open $dbsnp");
##ignore first 4 lines
<DBSNP>;
<DBSNP>;
<DBSNP>;
<DBSNP>;


##Now get the file pointer to the start of the DBSNP filesection abut the chromosome we're interested in
my $dbsnp_fptr_at_start_of_right_chrom = "false";
my $next_dbsnpline;
while ($dbsnp_fptr_at_start_of_right_chrom eq "false")
{
    $next_dbsnpline = <DBSNP>;
    if ($next_dbsnpline =~ /^$chrom\|/)
    {
	$dbsnp_fptr_at_start_of_right_chrom="true";
    }
}



while(<FILE>)
{
    my $line = $_;
    chomp $line;
    
    if ($line =~ /^VARIATION: (\d+), start coordinate (\d+)/)
    {
	my $varno = $1;
	my $varname = "chromosome_".$chrom."_VARIATION_".$varno;
	my $start = $2;
	
	##We must collect the following info: varname (done), type, size, filter, chrom (done), start (done), end, genotype, genotype confidence
	my $type="Q";
	my $size="Q";
	my $end="Q";
	my $genotype="Q";
	my $genotype_confidence="Q";
	my $rsid="";


	my $filter="";
	my $filter_reason="";
	if ($must_do_filtering eq "yes")
	{
	    if (exists $fivep_hash{$varname})
	    {
		if ($fivep_hash{$varname}==1)
		{
		    $filter="0";
		    $filter_reason="-";
		}
		else
		{

		    ## OK - the 5prime flank does not map uniquely. This qualifies the variant for filtering out UNLESS the other place where the 5prime flank maps to is in the trusted branch!
		    open(CHECK,"sv_called_in_chrom_".$chrom)||die("Cannot open sv_called_in_chrom $chrom");

		    while (<CHECK>)
		    {
			my $ln = $_;
			chomp $ln;
			my $tmp1="var_".$varno."_5p_flank";
			my $tmp2="var_".$varno."_trusted_branch";
			my $fivepseq;

			if ($ln =~ /$tmp1/)
			{
			    $fivepseq = <CHECK>;
			    chomp $fivepseq;
			    <CHECK>;
			    $ln = <CHECK>;##the trusted branch sequence
			    if ($ln =~ /$fivepseq/)
			    {
				## OK - we expect that sometimes the 5p flank to map also into the trusted branch - no surprise. Do not filter
				$filter="0";
				$filter_reason="warning_5pflank_maps_into_deletion"; #print warning as need to double check
			    }
			    else
			    {
				$filter="1"; ##does not map uniquely
				$filter_reason="5prime_flank_not_map_uniquely";				
			    }
			    last;
			}
		    }
		    close(CHECK);
		    


		}
	    }
	    else
	    {
		##5p flank does not map at all - this should be impossible
		print("This variant $varname, the 5prime flank appears not to map at all to the reference!!! Investigate\n");
		$filter="1";
		$filter_reason="5prime_flank_not_map_AT_ALL";
	    }
	}
	else
	{
	    $filter="-";
	    $filter_reason="-";
	}


	my $flank5p_info = <FILE>;
	chomp $flank5p_info;
	my $flank5p_seq= <FILE>;
	chomp $flank5p_seq;

	my $trusted_info=<FILE>;
	chomp $trusted_info;
	my $trusted_seq=<FILE>;
	chomp $trusted_seq;

	my $variant_info=<FILE>;
	chomp $variant_info;
	my $variant_seq=<FILE>;
	chomp $variant_seq;

	my $flank3p_info=<FILE>;
	chomp $flank3p_info;
	my $flank3p_seq=<FILE>;
	chomp $flank3p_seq;

	my $trusted_len;
	my $variant_len;

	if ($trusted_info =~ /length:(\d+)/)
	{
	    $trusted_len=$1;
	}
	else
	{
	    die("Problem with this line: $trusted_info");
	}

	if ($variant_info =~ /length:(\d+)/)
	{
	    $variant_len=$1;
	}
	else
	{
	    die("Problem with this line: $variant_info");
	}


	## ***********    classifying type - the only way to do this properly is with Needleman-Wunsch. The case we handle worst is very long balanced things, that may have only a couple of SNPs in them:

	if  ($trusted_len-$variant_len>32) 
	{
	    $type="deletion";
	    $size = $trusted_len - $variant_len;
	    $end =  $start+ $trusted_len - $variant_len;
	}

	elsif ($variant_len-$trusted_len>32)
	{

	    $type="insertion";
	    $size = $variant_len-$trusted_len;
	    $end =  $start;
	    $start--;
	}

	elsif (($trusted_len == $variant_len) && ($trusted_len == 32))
	{
	    $type="snp";
	    $size=1;
	    $end=$start;
	}
	elsif (($trusted_len == $variant_len) && ($trusted_len > 32))
	{
	    $type = "balanced_unresolved";
	    $end=$start+$trusted_len;
	    $size="<".$trusted_len;
	}
	elsif ( ($variant_len == $trusted_len) && ($trusted_len<32) )
	{
	    ##don't think this is possible
	    $type="snp or indel";
	    $end="$start+1";
            $size="?";
	}
	else
	{
	    $type="indel";
	    $end=$start+$trusted_len;
            $size=abs($trusted_len-$variant_len);

	}


	## First attempt at genotyping
	

	#printf("Looking at $trusted_info\n");
                               
	if ($trusted_info =~ /covgs of trusted not variant (.+)number_of_such_nodes: (\d+)/)
	{
	    my $trusted_nonvar_covgs = $1;
	    my $num = $2;
	    my $count_zeroes_in_trusted_nonvar=0;
	    my $total_trusted_non_var_covg=0;
	    my $total_var_non_trusted_covg=0;
	    my $var_non_trusted_covgs;

	    my $i=0;


	    my @split1 = split(/ /, $trusted_nonvar_covgs);
	    foreach my  $bit1 (@split1)
	    {
		if ($bit1 =~ /\d+/)
		{
		    if ($bit1==0)
		    {
			$count_zeroes_in_trusted_nonvar++;
		    }
		    $total_trusted_non_var_covg=$total_trusted_non_var_covg+$bit1;
		}
	    }


	    if ($variant_info =~  /covgs of variant not trusted(.+)number_of_such_nodes: (\d+)/)
	    {
		$var_non_trusted_covgs=$1;
	    }

	    if ($var_non_trusted_covgs)
	    {
		my @split = split(/ /, $var_non_trusted_covgs);
		foreach my  $bit (@split)
		{
		    if ($bit =~ /\d+/)
		    {
			$total_var_non_trusted_covg=$total_var_non_trusted_covg+$bit;		    
		    }
		}
	    }


	    ## All the nodes in the trusted path that are not in the variant path have coverage zero.
	    if ($count_zeroes_in_trusted_nonvar == $num)
	    {
		$genotype="hom";
		$genotype_confidence="1";
	    }
	    else
	    {##this is rubbish - improve
		if ($total_trusted_non_var_covg>$total_var_non_trusted_covg)
		{
		    $genotype="het";
		    $genotype_confidence="0.5";
		}
		elsif ($total_trusted_non_var_covg<$total_var_non_trusted_covg)
		{
		    $genotype="hom";
		    $genotype_confidence="0.5";
		}
		else
		{
		    $genotype="";
		    $genotype_confidence="";
		}
	    }
	}
	else
	{
	    #print "Cannot match $trusted_info\n";
	}




	## Some simplifications - if two branches are identical except at first base, then use that info to clean up output:
	
	my $init1 = substr($trusted_seq,0,1);
	my $init2 = substr($variant_seq,0,1);
	
	my $tmp1 = substr($trusted_seq,1,$trusted_len-1);
	my $tmp2 = substr($variant_seq,1,$variant_len-1);
	

	##The only simplification we can easily make is this. If we have a SNP, only the first character of variant_seq/trusted_seq matters
	##check this is true, then use it
	if ($type eq "snp")
	{

	    if ($tmp1 eq $tmp2)
	    {
		$trusted_seq = $init1;
		$variant_seq = $init2;
	    }
	    else
	    {
		print "Hold on. This isa a SNP, but the two sequences differ in more than just the frst character\n";
		print "$trusted_info\t$trusted_seq\n$variant_info\t$variant_seq\n";
		die();
	    }
	}



	## Finally - is it in dbSNP.
	# format is chr,chrPosFrom, chrPosTo, rs,ChrAllele,variantAllele,snpAlleleChrOrien, snp2chrOrien, snpClassAbbrev
	# eg 1|885|885|rs62636509|G|A|A/G|0|single base|1|2|2|3|reference||
	

	if ($type eq "snp")
	{

	    my $st = 0;
	    my $en = 0;
	    my $id;
	    my $ref_al;
	    my $var_al;
	    my $or;
	    my $dbsnp_type;

	    while ($st<$start)
	    {
		$next_dbsnpline = <DBSNP>;
		chomp $next_dbsnpline;

		my @cols = split(/\|/, $next_dbsnpline);
		
		#f ($next_dbsnpline=~ /$chrom\|(\d+)\|(\d+)\|([^\|]*)\|([^\|]*)\|([^\|]*)\|[^\|]*\|[^\|]*\|([^\|]*)\|/)
		if (scalar(@cols) >8)
		{
		    $st = $cols[1]; ##$1;
		    $en = $cols[2]; ##$2;
		    $id=  $cols[3]; ##$3;
		    $ref_al = $cols[4]; ##$4;
		    $var_al = $cols[5]; ##$5;
		    $dbsnp_type=$cols[8]; ##$6;

		}
		else
		{
		    die("Failed to match $next_dbsnpline");
		}

	    }

	    if ( ($st==$start) && ($en==$start) )
	    {
		#print "This dbsnp line seems to match start $start and end $end and variant $varname\n$next_dbsnpline\n";
		
		if ($dbsnp_type eq "single base")
		{
		    if (($ref_al =~ /$trusted_seq/)&& ($var_al =~ /$variant_seq/) )
		    {
			$rsid = $id;
		    }
		}
		else
		{
		    ##print STDERR "We this ref and var are $trusted_seq and $variant_seq, but they think they are $ref_al and $var_al\n$next_dbsnpline";
		    
		}
	    }
	}


	## What genes are in this variant
	my $slice = $sa->fetch_by_region("chromosome",$chrom,$start,$end);

	#my $var_feats = $slice->get_all_VariationFeatures;
	my $genes = $slice->get_all_Genes;
	my $rfeats = $ra->fetch_all_by_Slice($slice, undef);

	#my $v_features = (map {$_->display_id} @$var_feats);
	my $g = (map {$_->external_name} @$genes);
	my $r = (map {$_->display_id} @$rfeats);

		
	#print "#$varname\t$type\t$size\t$chrom\t$start\t$end\t$genotype\t$genotype_confidence\t$rsid\t$filter\t$filter_reason\t$final_cols\n";
	print join("\t", $varname,$type,$size,$chrom,$start,$end,$genotype,$genotype_confidence,$rsid,$filter,$filter_reason, join(",", (map {$_->external_name} @$genes)), join(",", (map {$_->display_id} @$rfeats))  );
	print "\n";
	print "Ref:$trusted_seq\n";
	print "Var:$variant_seq\n";

    }
}




close(FILE);
close(DBSNP);















