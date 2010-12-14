#This script takes a chromosome (at the moment from ensembl) to generate 3 modified chromosomes called A, B and C.This is useful to benchamark software with simulated data
#The modifications are: 
# - random SNPs
# - remove ALUs (at the moment taken from the ensembl annotation)  -- later to be replaced by a list of custom regions
# - generate deletions
# - complex regions - these are SNPs in a close range followed by a deletion


use strict;
# these lines point to Mario's perl modules for ensembl/bioperl 
use lib "/homes/marioc/perl-modules/ensembl/modules";
use lib "/homes/marioc/perl-modules/ensembl-variation/modules";
use lib "/homes/marioc/perl-modules/bioperl-live";
use Getopt::Long;

use Bio::EnsEMBL::Registry;



my %revcomp_hash=( "A" => "T",
                   "a" => "T",
                   "T" => "A",
                   "t" => "A",
                   "G" => "C",
                   "g" => "C",
                   "C" => "G",
                   "c" => "G",
		   "N" => "N");




my ($registry,$chromosome,$snp_ratio,@delspec,$repeats,@complex_regions, @invspec, $space, $kmer_size, $hom_file, $het_file);

&GetOptions(
            'registry|r:s'           => \$registry,
            'chr|c:s'                => \$chromosome,
	    'snps|m:i'               => \$snp_ratio,      # frequency: SNPs every $snp_ratio bases (positions for SNPs are sampled uniformly at this rate)
	    'del|d:s'                => \@delspec,        # specification of deletions. Format size_min-size_max:number, eg 2-100:20 means 20 insertions ranging from 2 to 100 in size
 	    'rep|s:i'                => \$repeats,        # number of ALU repeats to be removed
	    'complex|t:s'            => \@complex_regions,# number of complex regions (same format as deletions)
	    'inv|v:s'                => \@invspec,        # specification of inversions. Format as for deletions
	    'clear_space|l:i'        => \$space,          # how much space to ensure there is around variants       
	    'kmer_size|k:i'          => \$kmer_size,       # this is only used when printing the files of variants (flanks + branches) 
	    'hom_file|m:s'           => \$hom_file,       # output filename for homozygous (non-ref) variants
	    'het_file|h:s'           => \$het_file,       # output filename for heterozygous variants
           );

Bio::EnsEMBL::Registry->load_all($registry);
my $ensembl_core_db = Bio::EnsEMBL::Registry->get_DBAdaptor('human','core');
die "database not defined\n" if not defined $ensembl_core_db;
my $sa = $ensembl_core_db->get_SliceAdaptor;

## check args
if ( ($space < 2*$kmer_size) ) ##|| ($space < ($kmer_size + 100)) )
{
    die("Must have space variable greater than 2*kmer_size"); # and also > kmer_size+100. Just leaving space for flanks. Space is $space, kmer is $kmer_size")
}

#open files for chromosomes A, B and C
my $fa;
my $fb;
my $fc;

open($fa,">chromosome_${chromosome}_A") || die "cannot open file A\n";
open($fb,">chromosome_${chromosome}_B") || die "cannot open file B\n";
open($fc,">chromosome_${chromosome}_C") || die "cannot open file C\n";
print $fa ">chromosome_${chromosome}_A\n";
print $fb ">chromosome_${chromosome}_B\n";
print $fc ">chromosome_${chromosome}_C\n";
my $streamA;
my $streamB;
my $streamC;


my $slice = $sa->fetch_by_region("chromosome",$chromosome);
die "chromosome $chromosome not defined\n" if not defined $slice;
my $chromosome_length = $slice->length;

print "chromosome $chromosome length: $chromosome_length\n";


# The algorithm works as follows:
# the array @events has one entry per position in the chromosome with information about the event (SNP,REP,DEL,etc) affecting the position - this is to avoid overlaping events.
# - first "pass" generates the events
# - second "pass" prints the output for the 3 chromosomes - every event is considered in turn randomly selecting the chromosomes affectd (hom/het).  


my @event;

## Mark positions in the ref chromosome that have Ns as we dont want variants there
my %Ns=();



## Later on we will need to know the precise start/end positions of variants in order to get the flanking regions.
## These are coordinates in the reference - i.e. the original chromosome. 
my %first_base_of_event=();#in case of a complex event, this is the first SNP before the deletion
my %last_base_of_event=();
my %start_of_5p_flank_of_event=();
my %end_of_5p_flank_of_event=();
my %start_of_3p_flank_of_event=();#<<<<<<<< this depends on choice of kmer. So when we are asked to produce output for multiple kmer sizes (so we can test the same
                                  #         simulated varianst at multiple kmer_sizes, we will have to loop through adding different kmer_sizes
my %end_of_3p_flank_of_event=();
#my %end_of_variation_in_event=();


#######################################################################################
# FIRST PASS:  define the events, and ensure they do not overlap



#defined complex regions
print "define complex regions (SNPs+indels)\n";

my $which_complex_set="";
foreach my $complex (@complex_regions)
{
    my ($st,$ed,$number) = ($complex =~ /(\d+)\-(\d+)\:(\d+)/);
    $which_complex_set="COMPLEX_SET_".$st."_".$ed."_";
    
    if (defined $st and $st<=$ed){
	my $comp_count=0;
	print "generate $number complex regions each with a deletion sized somewhere between $st and $ed \n";
	
	while ($comp_count<$number){
	    
	    #get size of deletion
	    my $size = int(rand($ed-$st))+$st;
	    
	    #get position
	    my $pos = int(rand($chromosome_length-1))+1;
	    
	    my $overlaps_event = 0;
	    
	    my $i=$pos;
	    
	    #in a complex variant, there is a region directly before which is  prone to SNPs
	    my $length_snp_prone_reg = 20;
	    
	    # we need to ensure compex region + space on either side for flanking region does not overlap another variant
	    # so have <---space bases---> POS <--- $length snp prone reg ---> START DELETION....END DELETION<-----space bases----->
	    my $ctr;
	    #check $space bases before $pos.
	    for ($ctr = $pos-$space-1; ($ctr < $pos) &&(not $overlaps_event) ; $ctr++)
	    {
		if ( (defined $event[$ctr]) || (exists $Ns{$ctr}) )
		{
		    $overlaps_event=1;
		}
	    }
	    
	    #check region of variant plus $space bases afterwards. Is separate loop as involving $i here.
	    while( ($i<$pos+$size+$length_snp_prone_reg+$space) && not $overlaps_event){
		if ((defined $event[$i]) || (exists $Ns{$ctr}) )
		    {
			$overlaps_event = 1;
		    }
		$i++;
	    }
	    
	    
	    if (not $overlaps_event){
		
		#introduce SNPs
		#define how many SNPs before deletion
		my $snps = int(rand(3))+2;
		
		my %snps_pos;
		while ((keys %snps_pos)<$snps){
		    $snps_pos{int(rand($length_snp_prone_reg))}=1;
		}
		
		my $j=0;
		my $pos_first_snp=-1;#we will need this to get the 5prime flank
		my $pos_end_deletion = $pos+$size+$length_snp_prone_reg; #we will need this to get the 3prime flank    
		
		for($i=$pos;$i<$pos_end_deletion;$i++)
		{
		    if (defined $snps_pos{$j})
		    {
			#print "TEST $i C_SNP $j\n";
			$event[$i]=$which_complex_set."COM_SNP_".$comp_count;
			
			if ($pos_first_snp == -1)
			{
			    #then this is the first SNP
			    $pos_first_snp=$i;
			}
		    }
		    elsif ($j>=$length_snp_prone_reg)
		    {
			#print "TEST $i C_DEL $j\n";
			$event[$i]=$which_complex_set."COM_DEL_".$comp_count;
		    }
		    else
		    {
			#print "TEST $i C $j\n";
			$event[$i]=$which_complex_set."COM_".$comp_count;
		    }
		    $j++;
		}
		if ($pos_first_snp == -1)
		{
		    die("programmig error. No first snp coord");
		}
		my $name = $which_complex_set.$comp_count;
		$first_base_of_event{$name}       = $pos_first_snp;
		$last_base_of_event{$name}        = $pos_end_deletion ;
		$start_of_5p_flank_of_event{$name}= $pos-$space;
		$end_of_5p_flank_of_event{$name}= $pos_first_snp-1;
		$start_of_3p_flank_of_event{$name}= $pos_end_deletion+$kmer_size+1;
		$end_of_3p_flank_of_event{$name}= $pos_end_deletion+$space;
		$comp_count++;
		#print "Complex event generated: size $size with $snps snps\n";

	    }
	}
	
    }
    
}

#remove ALU repeats
print "select $repeats ALU repeats to be removed\n";
my $all_repeat_feats = $slice->get_all_RepeatFeatures;
my @alu_repeats;

foreach (@$all_repeat_feats){

    my $class = $_->repeat_consensus->repeat_class;
    if ($class eq 'SINE/Alu'){
	push @alu_repeats, $_;
    }
    
}

print scalar(@alu_repeats)," ALU repeats in chromosome $chromosome\n";
if (@alu_repeats>0){
    my $rep_count=0;
    while($rep_count<$repeats){
	my $rep_index = int(rand(scalar(@alu_repeats)));
	my $rep = $alu_repeats[$rep_index];
	my $overlaps_event=0;
	my $i=$rep->start;
	
	
	# we need to ensure repeat + space on either side for flanking region does not overlap another variant
	my $ctr;
	#check $space bases before $rep->start.
	for ($ctr = $rep->start-$space-1; ($ctr < $rep->start) &&(not $overlaps_event) ; $ctr++)
	{
	    if ((defined $event[$ctr]) || (exists $Ns{$ctr}) )
	    {
		$overlaps_event=1;
	    }
	}
	
	#check the bases of the repeat itself and $space bases after
	while($i<=$rep->end+$space && not $overlaps_event){
	    if ( (defined $event[$i]) || (exists $Ns{$ctr}) )
	    {
		$overlaps_event = 1;
	    }
	    $i++;
	}
	
	if (not $overlaps_event){
	    
	    for($i=$rep->start;$i<=$rep->end;$i++){
		#print "TEST add rep event at $i\n";
		$event[$i]="REP_".$rep_count."_".$rep->repeat_consensus->name;
	    }
	    $first_base_of_event{"REP_".$rep_count."_".$rep->repeat_consensus->name}=$rep->start;
	    $last_base_of_event{"REP_".$rep_count."_".$rep->repeat_consensus->name}=$rep->end;
	    $start_of_5p_flank_of_event{"REP_".$rep_count."_".$rep->repeat_consensus->name}=$rep->start-$space;
	    $end_of_5p_flank_of_event{"REP_".$rep_count."_".$rep->repeat_consensus->name}=$rep->start-1;
	    $start_of_3p_flank_of_event{"REP_".$rep_count."_".$rep->repeat_consensus->name}=$rep->end+$kmer_size+1;
	    $end_of_3p_flank_of_event{"REP_".$rep_count."_".$rep->repeat_consensus->name}=$rep->end+$space;
	    $rep_count++;
	}
	
    }
}


#create all the deletions
my $which_del_set="";
foreach my $del (@delspec){
    #format i-j:n eg 1-100:1000 means 1000 deletions ranging from 1 to 100 
    my ($st,$ed,$number) = ($del =~ /(\d+)\-(\d+)\:(\d+)/);

    $which_del_set="DEL_SET_".$st."_".$ed."_";

    if (defined $st and $st<=$ed){
	my $del_count=0;
	print "generate dels for $st-$ed:$number\n";
	
	while ($del_count<$number){
	    
	    #get size
	    my $size = int(rand($ed-$st))+$st;
	    
	    #get position
	    my $pos = int(rand($chromosome_length-1))+1;
	    my $overlaps_event = 0;
	    
	    my $i=$pos;
	    
	    # we need to ensure deletion + space on either side for flanking region does not overlap another variant
	    my $ctr;
	    #check $space bases before start
	    for ($ctr = $pos-$space-1; ($ctr < $pos) &&(not $overlaps_event) ; $ctr++)
	    {
		if ((defined $event[$ctr]) || (exists $Ns{$ctr}) )
		{
		    $overlaps_event=1;
		}
	    }
	    
	    #check deletion + $space bases after
	    while($i<$pos+$size+$space && not $overlaps_event){
		if ( (defined $event[$i]) || (exists $Ns{$ctr}) )
		{
		    $overlaps_event = 1;
		}
		$i++;
	    }
	    
	    if (not $overlaps_event){
		#print "TEST DEL $pos $size\n";
		for($i=$pos;$i<$pos+$size;$i++){
		    #print "TEST add del event at $i\n";
		    $event[$i]=$which_del_set."DEL_".$del_count;
		}
		my $name = $which_del_set."DEL_".$del_count;
		$first_base_of_event{$name}=$pos;
		$last_base_of_event{$name}=$pos+$size-1;
		$start_of_5p_flank_of_event{$name}=$pos-$space;
		$end_of_5p_flank_of_event{$name}=$pos-1;
		#$end_of_variation_in_event{$name}=$pos+$size;
		$start_of_3p_flank_of_event{$name}=$pos+$size+$kmer_size-1;#note subtracted one 
		$end_of_3p_flank_of_event{$name}=$pos+$size+$space;
		$del_count++;
		#print "Generated clean indel size $size\n";
	    }
	}
    }
}


## create the inversions
my $which_inv_set="";
foreach my $inv (@invspec){
    #format i-j:n eg 1-100:1000 means 1000 deletions ranging from 1 to 100 
    my ($st,$ed,$number) = ($inv =~ /(\d+)\-(\d+)\:(\d+)/);
    
    $which_inv_set="INV_SET_".$st."_".$ed."_";

    if (defined $st and $st<=$ed){
	my $inv_count=0;
	print "generate inversions for $st-$ed:$number\n";
	
	while ($inv_count<$number){
	    
	    #get size
	    my $size = int(rand($ed-$st))+$st;
	    
	    #get position
	    my $pos = int(rand($chromosome_length-1))+1;
	    my $overlaps_event = 0;
	    
	    my $i=$pos;
	    
	    # we need to ensure inversion + space on either side for flanking region does not overlap another variant
	    my $ctr;
	    #check $space bases before start
	    for ($ctr = $pos-$space-1; ($ctr < $pos) &&(not $overlaps_event) ; $ctr++)
	    {
		if ((defined $event[$ctr]) || (exists $Ns{$ctr}) )
		{
		    $overlaps_event=1;
		}
	    }
	    
	    #check inversion + $space bases after
	    while($i<$pos+$size+$space && not $overlaps_event){
		if ((defined $event[$i]) || (exists $Ns{$ctr}) )
		{
		    $overlaps_event = 1;
		}
		$i++;
	    }
	    
	    if (not $overlaps_event){
		#print "TEST INV $pos $size\n";
		for($i=$pos;$i<$pos+$size;$i++){
		    #print "TEST add inv event $inv_count at $i\n";
		    $event[$i]=$which_inv_set."INV_".$inv_count;
		}
		my $name = $which_inv_set."INV_".$inv_count;
		
		$first_base_of_event{$name}=$pos;
		$last_base_of_event{$name}=$pos+$size-1;
		$start_of_5p_flank_of_event{$name}=$pos-$space;
		$end_of_5p_flank_of_event{$name}=$pos-1;
		$start_of_3p_flank_of_event{$name}=$pos+$size+$kmer_size-1;##may need to check this -1
		$end_of_3p_flank_of_event{$name}=$pos+$size+$space;
		$inv_count++;
		#print "Generated clean inversion size $size\n";
	    }
	}
    }
}





#create all the SNPs
my $number_SNPs = int($chromosome_length/$snp_ratio);

print "generate $number_SNPs SNPs\n";

my $snp_count=0;
my $limit = 2000000; ## if you ask for more SNPs than we can fit, we try $limit times, but then quit
my $tries_so_far=0;

while($snp_count<$number_SNPs){

    my $pos = int(rand($chromosome_length+1)+1);

    ## ensure does not lie within $space of any other variant. 
    ## Ensures for example that if kmer_size is $space, a clean bubble will close off without touching 
    
    my $locus_is_well_away_from_other_vars="true";
    my $j;

    if ($tries_so_far>$limit)
    {
	die("Unable to fit so many SNPS in - only managed $snp_count out of $number_SNPs\n");
    }

    for ($j=$pos-$space-1; $j<=$pos+$space+1; $j++)
    {
	if (defined $event[$j])
	{
	    $locus_is_well_away_from_other_vars="false";
	    $tries_so_far++;
	}
    }

    if ($locus_is_well_away_from_other_vars eq "true")
    {
	#print "SNP $pos \n";
	$event[$pos]="SNP_".$snp_count; #1 is SNPs
	$first_base_of_event{$event[$pos]}=$pos;
	$last_base_of_event{$event[$pos]}=$pos;
	$start_of_5p_flank_of_event{$event[$pos]}=$pos-$space;
	$end_of_5p_flank_of_event{$event[$pos]}=$pos-1;
	$start_of_3p_flank_of_event{$event[$pos]}=$pos+$kmer_size+1;
	$end_of_3p_flank_of_event{$event[$pos]}=$pos+$space;
	
	$snp_count++;
    }
    
    
}


#******************************************************
# SECOND PASS :   print the chromosomes to file,


my %event_to_affected=();
my %event_to_5p_seq=();
my %event_to_b1_seq=();#b1 will be the reference branch ecept when it is a HOM non-ref variant
my %event_to_b2_seq=();
my %event_to_3p_seq=();
my %event_to_A_coord=();
my %event_to_B_coord=();
my %event_to_C_coord=();



my ($coordA,$coordB,$coordC) = (1,1,1);

my $i;
my $last_pos=1;

print "print chromosomes\n";


for($i=1;$i<=$chromosome_length;$i++)
{
    
    if ($event[$i]=~/SNP|DEL|REP|INV/){
	#print "TEST $last_pos ",$i-1,"\n";
	#print "EVENT $event[$i]\n";
	
	if ($last_pos>$i-1)
	{
	    die("last_pos $last_pos is > i-1, $-1");
	}
	my $subseq = $slice->subseq($last_pos,$i-1);
	$streamA = print_fasta($fa,$subseq,$streamA);
	$streamB = print_fasta($fb,$subseq,$streamB);
	$streamC = print_fasta($fc,$subseq,$streamC);
	$coordA += length($subseq);    
	$coordB += length($subseq);
	$coordC += length($subseq);
	
	
	if ( ($event[$i]=~/SNP/) && ($event[$i] !~ /COM/) ){ #SNP which is not in a complex event - ie isolated
	    
	    my $ref_base   = $slice->subseq($i,$i);
	    my $other_allele  = pick_other_allele($ref_base);

	    #sanity
	    sanity_check_args($event[$i], $i, $start_of_5p_flank_of_event{$event[$i]}, $end_of_5p_flank_of_event{$event[$i]} , $last_base_of_event{$event[$i]}, 
			      $start_of_3p_flank_of_event{$event[$i]}, $end_of_3p_flank_of_event{$event[$i]});

	    $event_to_5p_seq{$event[$i]}  = $slice->subseq($start_of_5p_flank_of_event{$event[$i]}, $end_of_5p_flank_of_event{$event[$i]});
	    $event_to_b1_seq{$event[$i]} = $slice->subseq($i, $start_of_3p_flank_of_event{$event[$i]}-1) ;#ref_branch    *****
	    $event_to_b2_seq{$event[$i]} = $other_allele.($slice->subseq($i+1, $start_of_3p_flank_of_event{$event[$i]}-1)) ;#variant branch  ******
	    $event_to_3p_seq{$event[$i]}= $slice->subseq($start_of_3p_flank_of_event{$event[$i]}, $end_of_3p_flank_of_event{$event[$i]});
	    $event_to_A_coord{$event[$i]}=$coordA;
	    $event_to_B_coord{$event[$i]}=$coordB;
	    $event_to_C_coord{$event[$i]}=$coordC;
	    
	    #print "SNP $i $ref_base $other_allele - A $coordA B $coordB C $coordC\n";
	    
	    $streamA = print_fasta($fa,$ref_base,$streamA);
	    
	    #select which haploid to affect - 0 is homozygous - 1 is het B - 2 is het C
	    my $sel = int(rand(3));
	    
	    my @affected;
	    if ($sel == 0){
		@affected = ('B','C');
		$streamB = print_fasta($fb,$other_allele,$streamB);
		$streamC = print_fasta($fc,$other_allele,$streamC);
		$event_to_affected{$event[$i]}="HOM_snpBC";
	    }
	    elsif ($sel == 1){
		@affected = ('B');
		$streamB = print_fasta($fb,$other_allele,$streamB);
		$streamC = print_fasta($fc,$ref_base,$streamC);
		$event_to_affected{$event[$i]}="HET_snpB";
	    }
	    elsif ($sel == 2){
		@affected = ('C');
		$streamB = print_fasta($fb,$ref_base,$streamB);
		$streamC = print_fasta($fc,$other_allele,$streamC);
		$event_to_affected{$event[$i]}="HET_snpC";
	    }


	    #print join (" ",'SNP',@affected,'-', $i, $ref_base, $other_allele,'-',"A $coordA","B $coordB","C $coordC",$event[$i]),"\n";
	    
	    
	    $coordA++;    
	    $coordB++;
	    $coordC++;
	    $last_pos = $i+1;
	    
	}
	elsif (($event[$i]=~/(DEL|REP)/) && ($event[$i] !~ /COM/) ){  ## repeat deletion or clean deletion  - but not complex event
	    my $event = $event[$i];
	    my $event_name;
	    if ($event =~ /DEL/)
	    {
		$event_name = "DEL";
	    }
	    elsif ($event =~ /REP/)
	    {
		$event_name = "REP";
	    }

	    sanity_check_args($event[$i], $i, $start_of_5p_flank_of_event{$event[$i]}, $end_of_5p_flank_of_event{$event[$i]}, $last_base_of_event{$event[$i]}, 
			      $start_of_3p_flank_of_event{$event[$i]}, $end_of_3p_flank_of_event{$event[$i]} );

	    $event_to_5p_seq{$event[$i]}  = $slice->subseq($start_of_5p_flank_of_event{$event[$i]}, $end_of_5p_flank_of_event{$event[$i]});
	    #$event_to_b1_seq{$event[$i]} = $slice->subseq($i, $start_of_3p_flank_of_event{$event[$i]}-1) ;#ref_branch
            #$event_to_b2_seq{$event[$i]} = $slice->subseq($last_base_of_event{$event[$i]}+1, $start_of_3p_flank_of_event{$event[$i]}-1  ); # variant branch
	    $event_to_3p_seq{$event[$i]}= $slice->subseq($start_of_3p_flank_of_event{$event[$i]}, $end_of_3p_flank_of_event{$event[$i]});


	    #select which haploid to affect - 0 A - 1 B - 2 C - 3 A+B - 4 A+C - 5 B+C
	    my $sel = int(rand(6));
	    my $st = $i;
	    
	    #retrieve region to be deleted
	    while ($event[$i+1] eq $event) {$i++};
	    
	    $last_pos = $i+1;
	    
	    if ($st>$i)
	    {
		die("st $st is > i $i - breaks slice->subseq\n");
		die();
	    }
	    my $subseq = $slice->subseq($st,$i);


	    if ( ($sel==0) || ($sel==3) || ($sel==4) ) #then we are deleting from A  - the ref
	    {
		$event_to_b2_seq{$event[$i]} = $slice->subseq($st, $start_of_3p_flank_of_event{$event[$i]}-1) ;#var_branch
		$event_to_b1_seq{$event[$i]} = $slice->subseq($last_base_of_event{$event[$i]}+1, $start_of_3p_flank_of_event{$event[$i]}-1  ); # ref branch
	    }
	    else
	    {
		$event_to_b1_seq{$event[$i]} = $slice->subseq($st, $start_of_3p_flank_of_event{$event[$i]}-1) ;#ref_branch
		$event_to_b2_seq{$event[$i]} = $slice->subseq($last_base_of_event{$event[$i]}+1, $start_of_3p_flank_of_event{$event[$i]}-1  ); # variant branch
	    }

	    if ($sel == 0){
		$streamB = print_fasta($fb,$subseq,$streamB);
		$streamC = print_fasta($fc,$subseq,$streamC);
		#print "$event_name A - $st $i - A $coordA B $coordB C $coordC $event\n";
		$event_to_A_coord{$event[$i]}=$coordA;
		$event_to_B_coord{$event[$i]}=$coordB;
		$event_to_C_coord{$event[$i]}=$coordC;

		$coordB += length($subseq);
		$coordC += length($subseq);
		$event_to_affected{$event[$i]}="HOM_delA";
		##note - deletion is in reference here:
		#$event_to_b2_seq{$event[$i]} = $slice->subseq($i, $start_of_3p_flank_of_event{$event[$i]}-1) ;#var_branch
		#$event_to_b1_seq{$event[$i]} = $slice->subseq($last_base_of_event{$event[$i]}+1, $start_of_3p_flank_of_event{$event[$i]}-1  ); # ref branch
		#$event_to_b2_seq{$event[$i]} =$event_to_b2_seq{$event[$i]}.$subseq;

	    }
	    elsif ($sel == 1){
		$streamA = print_fasta($fa,$subseq,$streamA);
		$streamC = print_fasta($fc,$subseq,$streamC);
		#print "$event_name B - $st $i - A $coordA B $coordB C $coordC $event\n";
		$event_to_A_coord{$event[$i]}=$coordA;
		$event_to_B_coord{$event[$i]}=$coordB;
		$event_to_C_coord{$event[$i]}=$coordC;
		$coordA += length($subseq);    
		$coordC += length($subseq);
		$event_to_affected{$event[$i]}="HET_delB";
		#$event_to_b1_seq{$event[$i]} = $slice->subseq($i, $start_of_3p_flank_of_event{$event[$i]}-1) ;#ref_branch
		#$event_to_b2_seq{$event[$i]} = $slice->subseq($last_base_of_event{$event[$i]}+1, $start_of_3p_flank_of_event{$event[$i]}-1  ); # variant branch		
		#$event_to_b1_seq{$event[$i]} =$event_to_b1_seq{$event[$i]}.$subseq;

	    }
	    elsif ($sel == 2){
		$streamA = print_fasta($fa,$subseq,$streamA);
		$streamB = print_fasta($fb,$subseq,$streamB);
		#print "$event_name C - $st $i - A $coordA B $coordB C $coordC $event\n";
		$event_to_A_coord{$event[$i]}=$coordA;
		$event_to_B_coord{$event[$i]}=$coordB;
		$event_to_C_coord{$event[$i]}=$coordC;
		$coordA += length($subseq);    
		$coordB += length($subseq);
		$event_to_affected{$event[$i]}="HET_delC";	
		#$event_to_b1_seq{$event[$i]} = $slice->subseq($i, $start_of_3p_flank_of_event{$event[$i]}-1) ;#ref_branch
		#$event_to_b2_seq{$event[$i]} = $slice->subseq($last_base_of_event{$event[$i]}+1, $start_of_3p_flank_of_event{$event[$i]}-1  ); # variant branch		
		#$event_to_b1_seq{$event[$i]} =$event_to_b1_seq{$event[$i]}.$subseq;


	    }   
	    elsif ($sel == 3){
		$streamC = print_fasta($fc,$subseq,$streamC);
		#print "$event_name A B  - $st $i - A $coordA B $coordB C $coordC $event\n";
		$event_to_A_coord{$event[$i]}=$coordA;
		$event_to_B_coord{$event[$i]}=$coordB;
		$event_to_C_coord{$event[$i]}=$coordC;
		$coordC += length($subseq);
		$event_to_affected{$event[$i]}="HET_delAB";
		##note - deletion is in reference aswell as B here:
		#$event_to_b2_seq{$event[$i]} = $slice->subseq($i, $start_of_3p_flank_of_event{$event[$i]}-1) ;#var_branch
		#$event_to_b1_seq{$event[$i]} = $slice->subseq($last_base_of_event{$event[$i]}+1, $start_of_3p_flank_of_event{$event[$i]}-1  ); # ref branch
		#$event_to_b2_seq{$event[$i]} =$event_to_b2_seq{$event[$i]}.$subseq;
		

	    }   
	    elsif ($sel == 4){
		$streamB = print_fasta($fb,$subseq,$streamB);
		#print "$event_name A C - $st $i - A $coordA B $coordB C $coordC $event\n";
		$event_to_A_coord{$event[$i]}=$coordA;
		$event_to_B_coord{$event[$i]}=$coordB;
		$event_to_C_coord{$event[$i]}=$coordC;
		$coordB += length($subseq);
		$event_to_affected{$event[$i]}="HET_delAC";
		##note - deletion is in reference aswell as C here:
		#$event_to_b2_seq{$event[$i]} = $slice->subseq($i, $start_of_3p_flank_of_event{$event[$i]}-1) ;#var_branch
		#$event_to_b1_seq{$event[$i]} = $slice->subseq($last_base_of_event{$event[$i]}+1, $start_of_3p_flank_of_event{$event[$i]}-1  ); # ref branch
		#$event_to_b2_seq{$event[$i]} =$event_to_b2_seq{$event[$i]}.$subseq;

	    } 
	    elsif ($sel == 5){
		$streamA = print_fasta($fa,$subseq,$streamA);
		#print "$event_name B C - $st $i - A $coordA B $coordB C $coordC $event\n";
		$event_to_A_coord{$event[$i]}=$coordA;
		$event_to_B_coord{$event[$i]}=$coordB;
		$event_to_C_coord{$event[$i]}=$coordC;
		$coordA += length($subseq);
		$event_to_affected{$event[$i]}="HOM_delBC";
		#$event_to_b1_seq{$event[$i]} = $slice->subseq($i, $start_of_3p_flank_of_event{$event[$i]}-1) ;#ref_branch
		#$event_to_b2_seq{$event[$i]} = $slice->subseq($last_base_of_event{$event[$i]}+1, $start_of_3p_flank_of_event{$event[$i]}-1  ); # variant branch		
		#$event_to_b2_seq{$event[$i]} =$event_to_b2_seq{$event[$i]}.$subseq;
		
		
	    } 
	    
	}
	elsif (($event[$i]=~/INV/) && ($event[$i] !~ /COM/))##clean inversion
	{

	    my $event = $event[$i];
	    my $event_name = "INV";

	    sanity_check_args($event[$i], $i, $start_of_5p_flank_of_event{$event[$i]}, $end_of_5p_flank_of_event{$event[$i]}, $last_base_of_event{$event[$i]}, 
			      $start_of_3p_flank_of_event{$event[$i]}, $end_of_3p_flank_of_event{$event[$i]} );

	    $event_to_5p_seq{$event[$i]}  = $slice->subseq($start_of_5p_flank_of_event{$event[$i]}, $end_of_5p_flank_of_event{$event[$i]});
	    #$event_to_b1_seq{$event[$i]} = $slice->subseq($i, $start_of_3p_flank_of_event{$event[$i]}-1) ;#ref_branch
            #$event_to_b2_seq{$event[$i]} = $slice->subseq($last_base_of_event{$event[$i]}+1, $start_of_3p_flank_of_event{$event[$i]}-1  ); # variant branch
	    $event_to_3p_seq{$event[$i]}= $slice->subseq($start_of_3p_flank_of_event{$event[$i]}, $end_of_3p_flank_of_event{$event[$i]});


	    #select which haploid to affect - 0 A - 1 B - 2 C - 3 A+B - 4 A+C - 5 B+C
	    my $sel = int(rand(6));
	    my $st = $i;
	    
	    #retrieve region to be inverted
	    while ($event[$i+1] eq $event) {$i++};
	    
	    $last_pos = $i+1;
	    
	    if ($st>$i)
	    {
		die("st $st is > i $i - breaks slice->subseq\n");
		die();
	    }
	    my $subseq = $slice->subseq($st,$i);
	    my $revcomp_subseq = get_reverse_complement_of_clean_sequence($subseq);

	    if ( ($sel==0) || ($sel==3) || ($sel==4) ) #then we are inverting  from A  - the ref
	    {
		$event_to_b2_seq{$event[$i]} = $slice->subseq($st, $start_of_3p_flank_of_event{$event[$i]}-1) ;#var_branch
		$event_to_b1_seq{$event[$i]} = $revcomp_subseq.($slice->subseq($i+1, $start_of_3p_flank_of_event{$event[$i]}-1 )); # ref branch
	    }
	    else
	    {
		$event_to_b1_seq{$event[$i]} = $slice->subseq($st, $start_of_3p_flank_of_event{$event[$i]}-1) ;#ref_branch
		$event_to_b2_seq{$event[$i]} = $revcomp_subseq.($slice->subseq($i+1, $start_of_3p_flank_of_event{$event[$i]}-1)); # variant branch
	    }

	    if ($sel == 0){
		$streamA = print_fasta($fa,$revcomp_subseq,$streamA);
		$streamB = print_fasta($fb,$subseq,$streamB);
		$streamC = print_fasta($fc,$subseq,$streamC);
		#print "$event_name A - $st $i - A $coordA B $coordB C $coordC $event\n";
		$event_to_A_coord{$event[$i]}=$coordA;
		$event_to_B_coord{$event[$i]}=$coordB;
		$event_to_C_coord{$event[$i]}=$coordC;
		$coordA += length($revcomp_subseq);
		$coordB += length($subseq);
		$coordC += length($subseq);
		$event_to_affected{$event[$i]}="HOM_invA";

	    }
	    elsif ($sel == 1){
		$streamA = print_fasta($fa,$subseq,$streamA);
		$streamB = print_fasta($fb,$revcomp_subseq,$streamB);
		$streamC = print_fasta($fc,$subseq,$streamC);
		#print "$event_name B - $st $i - A $coordA B $coordB C $coordC $event\n";
		$event_to_A_coord{$event[$i]}=$coordA;
		$event_to_B_coord{$event[$i]}=$coordB;
		$event_to_C_coord{$event[$i]}=$coordC;
		$coordA += length($subseq);    
		$coordB += length($revcomp_subseq);    
		$coordC += length($subseq);
		$event_to_affected{$event[$i]}="HET_invB";
	    }
	    elsif ($sel == 2){
		$streamA = print_fasta($fa,$subseq,$streamA);
		$streamB = print_fasta($fb,$subseq,$streamB);
		$streamC = print_fasta($fc,$subseq,$streamC);
		#print "$event_name C - $st $i - A $coordA B $coordB C $coordC $event\n";
		$event_to_A_coord{$event[$i]}=$coordA;
		$event_to_B_coord{$event[$i]}=$coordB;
		$event_to_C_coord{$event[$i]}=$coordC;
		$coordA += length($subseq);    
		$coordB += length($subseq);
		$coordC += length($revcomp_subseq);
		$event_to_affected{$event[$i]}="HET_invC";	
	    }   
	    elsif ($sel == 3){
		$streamA = print_fasta($fa,$revcomp_subseq,$streamA);
		$streamB = print_fasta($fb,$revcomp_subseq,$streamB);
		$streamC = print_fasta($fc,$subseq,$streamC);
		#print "$event_name A B  - $st $i - A $coordA B $coordB C $coordC $event\n";
		$event_to_A_coord{$event[$i]}=$coordA;
		$event_to_B_coord{$event[$i]}=$coordB;
		$event_to_C_coord{$event[$i]}=$coordC;
		$coordA += length($revcomp_subseq);
		$coordB += length($revcomp_subseq);
		$coordC += length($subseq);
		$event_to_affected{$event[$i]}="HET_invAB";

	    }   
	    elsif ($sel == 4){
		$streamA = print_fasta($fa,$revcomp_subseq,$streamA);
		$streamB = print_fasta($fb,$subseq,$streamB);
		$streamC = print_fasta($fc,$revcomp_subseq,$streamC);
		#print "$event_name A C - $st $i - A $coordA B $coordB C $coordC $event\n";
		$event_to_A_coord{$event[$i]}=$coordA;
		$event_to_B_coord{$event[$i]}=$coordB;
		$event_to_C_coord{$event[$i]}=$coordC;
		$coordA += length($revcomp_subseq);
		$coordB += length($subseq);
		$coordC += length($revcomp_subseq);
		$event_to_affected{$event[$i]}="HET_delAC";
	    } 
	    elsif ($sel == 5){
		$streamA = print_fasta($fa,$subseq,$streamA);
		$streamB = print_fasta($fb,$revcomp_subseq,$streamB);
		$streamC = print_fasta($fc,$revcomp_subseq,$streamC);
		#print "$event_name B C - $st $i - A $coordA B $coordB C $coordC $event\n";
		$event_to_A_coord{$event[$i]}=$coordA;
		$event_to_B_coord{$event[$i]}=$coordB;
		$event_to_C_coord{$event[$i]}=$coordC;
		$coordA += length($subseq);
		$coordB += length($revcomp_subseq);
		$coordC += length($revcomp_subseq);

		$event_to_affected{$event[$i]}="HOM_delBC";
		
	    } 

	}
	elsif ($event[$i] =~ /COM/) 
	{
	    ## idea is choose for the whole complex event, which haploids to apply it to.
	    ## ie the SNPs and the deletion will have identical zygosity

	    #select which haploid to affect - 0 A - 1 B - 2 C - 3 A+B - 4 A+C - 5 B+C
	    my $sel = int(rand(6));
	    my $st = $i;
	    my $event = $event[$i];

	    #sanity
	    my $overall_event_name = get_main_com_name($event[$i]);
	    sanity_check_args($overall_event_name, $i, $start_of_5p_flank_of_event{$overall_event_name}, $end_of_5p_flank_of_event{$overall_event_name} , 
			      $last_base_of_event{$overall_event_name}, 
			      $start_of_3p_flank_of_event{$overall_event_name}, $end_of_3p_flank_of_event{$overall_event_name});


	    # we can work these out straight away
	    $event_to_5p_seq{$overall_event_name}  = $slice->subseq($start_of_5p_flank_of_event{$overall_event_name}, $end_of_5p_flank_of_event{$overall_event_name});
	    $event_to_3p_seq{$overall_event_name}= $slice->subseq($start_of_3p_flank_of_event{$overall_event_name}, $end_of_3p_flank_of_event{$overall_event_name});

	    #$event_to_b1_seq{$overall_event_name} = $slice->subseq($i, $start_of_3p_flank_of_event{$overall_event_name}-1) ;#ref_branch
	    #$event_to_b1_seq{$overall_event_name} = $slice->subseq($end_of_5p_flank_of_event{$overall_event_name}+1, $start_of_3p_flank_of_event{$overall_event_name}-1) ;#ref_branch

	    #we will keep appending bases on the end of this as we work through the complex region.
	    $event_to_b1_seq{$overall_event_name} = "";
            $event_to_b2_seq{$overall_event_name} = "";


	    my $have_we_hit_first_thing_in_complex_event=0;

	    #move through the complex region
            while (complex_equal($event[$i+1], $event)==1)
	    {
		if ( ($event[$i] !~ /SNP/) && ($event[$i] !~ /DEL/)  )
		{
		    #this is a base within the complex region, but not a SNP or del - ie is same as reference
		    my $base = $slice->subseq($i,$i);
		    $streamA = print_fasta($fa,$base,$streamA);
		    $streamB = print_fasta($fb,$base,$streamB);
		    $streamC = print_fasta($fc,$base,$streamC);
		    $coordA++;    
		    $coordB++;
		    $coordC++;
		    $event_to_b1_seq{$overall_event_name}=$event_to_b1_seq{$overall_event_name}.$base;
		    $event_to_b2_seq{$overall_event_name}=$event_to_b2_seq{$overall_event_name}.$base;
		}
		elsif ($event[$i] =~ /SNP/)
		{
		    
		    my $ref_base   = $slice->subseq($i,$i);
		    my $other_allele  = pick_other_allele($ref_base);
		    #$event_to_b1_seq{$overall_event_name}=$event_to_b1_seq{$overall_event_name}.$ref_base;
		    #$event_to_b2_seq{$overall_event_name}=$event_to_b2_seq{$overall_event_name}.$other_allele;
		    if ($have_we_hit_first_thing_in_complex_event==0)
		    {
			$have_we_hit_first_thing_in_complex_event=1;
			$event_to_A_coord{get_main_com_name($event[$i])}=$coordA;
			$event_to_B_coord{get_main_com_name($event[$i])}=$coordB;
			$event_to_C_coord{get_main_com_name($event[$i])}=$coordC;
		    }
		    my @affected;
		    if ($sel == 0){ ## SNP applies to A only
			@affected = ('A');
			$streamA = print_fasta($fa,$other_allele,$streamA);
			$streamB = print_fasta($fb,$ref_base,$streamB);
			$streamC = print_fasta($fc,$ref_base,$streamC);
			$event_to_affected{$overall_event_name}="HOM_complexA"; 
			#confusingly, the SNP is in chr A, so branch1 is what will later be seen as the ref allele = what is in chrA
			# - so the other_allele goes in br1, as it is what is in chrA
                        $event_to_b1_seq{$overall_event_name}=$event_to_b1_seq{$overall_event_name}.$other_allele;
			$event_to_b2_seq{$overall_event_name}=$event_to_b2_seq{$overall_event_name}.$ref_base;
		    
		    }
		    elsif ($sel == 1){## SNP applies to B only
			@affected = ('B');
			$streamA = print_fasta($fa,$ref_base,$streamA);
			$streamB = print_fasta($fb,$other_allele,$streamB);
			$streamC = print_fasta($fc,$ref_base,$streamC);
			$event_to_affected{$overall_event_name}="HET_complexB";
			$event_to_b1_seq{$overall_event_name}=$event_to_b1_seq{$overall_event_name}.$ref_base;
			$event_to_b2_seq{$overall_event_name}=$event_to_b2_seq{$overall_event_name}.$other_allele;
			
		    }
		    elsif ($sel == 2){
			@affected = ('C');## SNP applies to C only
			$streamA = print_fasta($fa,$ref_base,$streamA);
			$streamB = print_fasta($fb,$ref_base,$streamB);
			$streamC = print_fasta($fc,$other_allele,$streamC);
			$event_to_affected{$overall_event_name}="HET_complexC";
			$event_to_b1_seq{$overall_event_name}=$event_to_b1_seq{$overall_event_name}.$ref_base;
			$event_to_b2_seq{$overall_event_name}=$event_to_b2_seq{$overall_event_name}.$other_allele;

		    }
		    elsif ($sel == 3){
			@affected = ('A', 'B');## SNP applies to A and B
			$streamA = print_fasta($fa,$other_allele,$streamA);
			$streamB = print_fasta($fb,$other_allele,$streamB);
			$streamC = print_fasta($fc,$ref_base,$streamC);
			$event_to_affected{$overall_event_name}="HET_complexAB";
			#confusingly, the SNP is in chr A, so branch1 is what will later be seen as the ref allele = what is in chrA
			# - so the other_allele goes in br1, as it is what is in chrA
                        $event_to_b1_seq{$overall_event_name}=$event_to_b1_seq{$overall_event_name}.$other_allele;
			$event_to_b2_seq{$overall_event_name}=$event_to_b2_seq{$overall_event_name}.$ref_base;

		    }

		    elsif ($sel == 4){
			@affected = ('A', 'C');## SNP applies to A and C
			$streamA = print_fasta($fa,$other_allele,$streamA);
			$streamB = print_fasta($fb,$ref_base,$streamB);
			$streamC = print_fasta($fc,$other_allele,$streamC);
			$event_to_affected{$overall_event_name}="HET_complexAC";
			#confusingly, the SNP is in chr A, so branch1 is what will later be seen as the ref allele = what is in chrA
			# - so the other_allele goes in br1, as it is what is in chrA
                        $event_to_b1_seq{$overall_event_name}=$event_to_b1_seq{$overall_event_name}.$other_allele;
			$event_to_b2_seq{$overall_event_name}=$event_to_b2_seq{$overall_event_name}.$ref_base;

		    }

		    elsif ($sel == 5){
			@affected = ('B', 'C');## SNP applies to B and C
			$streamA = print_fasta($fa,$ref_base,$streamA);
			$streamB = print_fasta($fb,$other_allele,$streamB);
			$streamC = print_fasta($fc,$other_allele,$streamC);
			$event_to_affected{$overall_event_name}="HET_complexBC";
			$event_to_b1_seq{$overall_event_name}=$event_to_b1_seq{$overall_event_name}.$ref_base;
			$event_to_b2_seq{$overall_event_name}=$event_to_b2_seq{$overall_event_name}.$other_allele;

		    }


		    #print join (" ",'COMPLEX SNP',@affected,'-', $i, $ref_base, $other_allele,'-',"A $coordA","B $coordB","C $coordC",$event[$i]),"\n";
		    

		    $coordA++;    
		    $coordB++;
		    $coordC++;
		    
		}
		elsif ($event[$i] =~ /DEL/)#deletion within complex region
		{

		    my $event = $event[$i];
		    
		    #retrieve region to be deleted
		    while ($event[$i+1] eq $event) {$i++};
		    
		    if ($st>$i)
		    {
			die("In complex del. st is $st and is > than i $i, breaks slice -> subseq");
		    }
		    my $subseq = $slice->subseq($st,$i);#remember we set $st=$i at the beginning of this if-statement - ie start of this region

		    # Reminder - we want the sequence that ends up in chrA to be in br1
		    #$event_to_b1_seq{$overall_event_name}=$event_to_b1_seq{$overall_event_name}.$subseq;			
		    # Don't append anything to the variant branch - this bit is deleted :-)

		    if ($have_we_hit_first_thing_in_complex_event==0)
		    {
			$have_we_hit_first_thing_in_complex_event=1;
			$event_to_A_coord{get_main_com_name($event[$i])}=$coordA;
			$event_to_B_coord{get_main_com_name($event[$i])}=$coordB;
			$event_to_C_coord{get_main_com_name($event[$i])}=$coordC;
		    }

		    if ($sel == 0){
			#deleted from A - print nothing to A
			$streamB = print_fasta($fb,$subseq,$streamB);
			$streamC = print_fasta($fc,$subseq,$streamC);
			#print "DEL A - $st $i - A $coordA B $coordB C $coordC $event\n";
			$coordB += length($subseq);
			$coordC += length($subseq);
			$event_to_affected{$overall_event_name}="HOM_complexA";
			$event_to_b2_seq{$overall_event_name}=$event_to_b2_seq{$overall_event_name}.$subseq;#deletion is in chrA = ref			

		    }
		    elsif ($sel == 1){
			$streamA = print_fasta($fa,$subseq,$streamA);
			#deleted from B - print nothing to B
			$streamC = print_fasta($fc,$subseq,$streamC);
			#print "DEL B - $st $i - A $coordA B $coordB C $coordC $event\n";
			$coordA += length($subseq);    
			$coordC += length($subseq);
			$event_to_affected{$overall_event_name}="HET_complexB";
			$event_to_b1_seq{$overall_event_name}=$event_to_b1_seq{$overall_event_name}.$subseq;			
			
		    }
		    elsif ($sel == 2){
			$streamA = print_fasta($fa,$subseq,$streamA);
			$streamB = print_fasta($fb,$subseq,$streamB);
			#print "DEL C - $st $i - A $coordA B $coordB C $coordC $event\n";
			$coordA += length($subseq);    
			$coordB += length($subseq);
			$event_to_affected{$overall_event_name}="HET_complexC";	
			$event_to_b1_seq{$overall_event_name}=$event_to_b1_seq{$overall_event_name}.$subseq;			

		    }   
		    elsif ($sel == 3){#delete from A and B - will look like an insertion in C
			$streamC = print_fasta($fc,$subseq,$streamC);
			#print "DEL A B  - $st $i - A $coordA B $coordB C $coordC $event\n";
			$coordC += length($subseq);
			$event_to_affected{$overall_event_name}="HET_complexAB";
			$event_to_b2_seq{$overall_event_name}=$event_to_b2_seq{$overall_event_name}.$subseq;			

		    }   
		    elsif ($sel == 4){#delete from A and C - will look like an insertion in B
			$streamB = print_fasta($fb,$subseq,$streamB);
			#print "DEL A C - $st $i - A $coordA B $coordB C $coordC $event\n";
			$coordB += length($subseq);
			$event_to_affected{$overall_event_name}="HET_complexAC";
			$event_to_b2_seq{$overall_event_name}=$event_to_b2_seq{$overall_event_name}.$subseq;			

		    } 
		    elsif ($sel == 5){#delete from B and C - will look like a hom non ref deletion in B and C
			$streamA = print_fasta($fa,$subseq,$streamA);
			#print "DEL B C - $st $i - A $coordA B $coordB C $coordC $event\n";
			$coordA += length($subseq);
			$event_to_affected{$overall_event_name}="HOM_complexBC";
			$event_to_b1_seq{$overall_event_name}=$event_to_b1_seq{$overall_event_name}.$subseq;			

		    } 
		    
		 
		    
		}
		
		$i++;
	    }

	    ##Here you need to add bases to b1_seq and b2_seq from the end of the complex region to the star of the 3p flank? 
	    if ($last_base_of_event{$overall_event_name}<$start_of_3p_flank_of_event{$overall_event_name}-1)
	    {
		my $added_seq = $slice->subseq($last_base_of_event{$overall_event_name}+1, $start_of_3p_flank_of_event{$overall_event_name}-1);
		$event_to_b1_seq{$overall_event_name}=$event_to_b1_seq{$overall_event_name}.$added_seq;
		$event_to_b2_seq{$overall_event_name}=$event_to_b2_seq{$overall_event_name}.$added_seq;
	    }

            $last_pos = $i+1;

	}
	
    }
}

#print "TEST $last_pos ",$i-1,"\n";
if ($last_pos>$i-1)
{
    print "Unexpected. last pos $last_pos > i-1, ";
    print $i-1;
    die();
}

my $subseq = $slice->subseq($last_pos,$i-1);
$streamA = print_fasta($fa,$subseq,$streamA);
$streamB = print_fasta($fb,$subseq,$streamB);
$streamC = print_fasta($fc,$subseq,$streamC);

print $fa $streamA;
print $fb $streamB;
print $fc $streamC;

print $fa "\n";
print $fb "\n";
print $fc "\n";






### Third pass, print out a pair of variant files, for hom and het variants


($coordA,$coordB,$coordC) = (1,1,1);
$last_pos=1;
print "print variants\n";


my %vars_printed_already;

open(HOM, "> ".$hom_file)||die("Cannot open $hom_file");
open(HET, "> ".$het_file)||die("Cannot open $het_file");

for($i=1;$i<=$chromosome_length;$i++)
{
    
    #if there is any event
    if ($event[$i]=~/SNP|DEL|REP|INV/)
    {
	#print "TEST $last_pos ",$i-1,"\n";
	#print "EVENT $event[$i]\n";
	

	#get the name for the overall event. 
	#ie a complex deletion will contain a deletion + SNPs labelled (eg) COM_DEL_23, and COM_SNP_23. For simplicity, I'm going to call these all COM_23 for this
	my $curr_event_name=get_main_com_name($event[$i]);
	my $homhetinfo = $event_to_affected{$curr_event_name};
	
	#print 5prime flank
	if (!exists $event_to_affected{$curr_event_name})
	{
	    die("Do not have event to affected for $curr_event_name");
	}

	my $info = $curr_event_name."_A_coord:".$event_to_A_coord{$curr_event_name}."_B_coord:".$event_to_B_coord{$curr_event_name}."_C_coord:".$event_to_C_coord{$curr_event_name}."_type_".$homhetinfo;

	if ($event_to_affected{$curr_event_name} =~ /^HET/)
	{
	    if (!exists $vars_printed_already{$curr_event_name})
	    {
		print HET ">var_5p_flank_".$info."\n";
		print HET $event_to_5p_seq{$curr_event_name};
		print HET "\n";
		print HET ">var_b1_ref_branch_".$info."\n";
		print HET $event_to_b1_seq{$curr_event_name};
		print HET "\n";
		print HET ">var_b2_var_branch_".$info."\n";
		print HET $event_to_b2_seq{$curr_event_name};
		print HET "\n";
		print HET ">var_3p_flank_".$info."\n";
		print HET $event_to_3p_seq{$curr_event_name};
		print HET "\n";
		$vars_printed_already{$curr_event_name}=1;
	    }
	}
	elsif ($event_to_affected{$curr_event_name} =~ /^HOM/)
	{

	    if (!exists $vars_printed_already{$curr_event_name})
	    {
		print HOM ">var_5p_flank_".$info."\n";
		print HOM $event_to_5p_seq{$curr_event_name};
		print HOM "\n";
		print HOM ">var_b1_ref_branch_".$info."\n";
		print HOM $event_to_b1_seq{$curr_event_name};
		print HOM "\n";
		print HOM ">var_b2_var_branch_".$info."\n";
		print HOM $event_to_b2_seq{$curr_event_name};
		print HOM "\n";
		print HOM ">var_3p_flank_".$info."\n";
		print HOM $event_to_3p_seq{$curr_event_name};
		print HOM "\n";
		$vars_printed_already{$curr_event_name}=1;
	    }


	}
	
	
    } 
    
}

close(HOM);
close(HET);
print "Completed\n";




sub pick_other_allele{
   my %num2base = (0 => 'A', 1 => 'C', 2 => 'G', 3 => 'T');

   my ($allele) = @_;
   my $index;
   do{
     $index = int(rand(4));

  } until (uc $allele ne $num2base{$index});

  return $num2base{$index};

}


sub print_fasta{
  my ($f,$string,$stream) = @_;
 
  my $seq =$stream.$string;
  
  my $tail_length = length($seq) % 80;

  my ($print,$return) = ($seq =~ /(\w*)(\w{$tail_length})$/);
  
  $print =~ s/(\w{80})/$1\n/ig;

  print $f $print;
  return $return;
  
}



#returns COM_xxx from COM_xxx, COM_SNP_xxx and COM_DEL_xxx
sub get_main_com_name
{
    my $a = $_[0];
    if ($a !~ /COM/)
    {
	return $a;
    }

    if ($a =~ /^(COMPLEX_SET_\d+_\d+_)COM_(\d+)$/)
    {
	return $1.$2;
    }
    elsif ($a =~ /^(COMPLEX_SET_\d+_\d+_)COM_SNP_(\d+)$/)
    {
	return $1.$2;
    }
    elsif ($a =~ /^(COMPLEX_SET_\d+_\d+_)COM_DEL_(\d+)$/)
    {
	return $1.$2;
    }
    else
    {
	die("Cannot get main com name of $a");
    }
    
}

#return 1 if true
sub complex_equal
{
    my ($a, $b) = @_;
    
    # Given two events, we want to knwo if they are the same. We could compare strings
    # except we want to consider SNPs and DEL within the same complex event as the same event.
    ## we just need to look at the number of the event
    if ($a eq $b)
    {
	return 1;
    }
    elsif (($a =~ /COM/) && ($b =~ /COM/)) ## can be COM_xxx or COM_SNP_xxx or COM_DE_xxx where xxx is the number of the event
    {
	if ($a =~ /(\d+)$/)
	{
	    my $a_num = $1;
	    my $b_num;
	    if ($b =~ /(\d+)$/)
	    {
		$b_num=$1;
	    }
	    else
	    {
		die("Cannot parse $b");
	    }
	    if ($a_num==$b_num)
	    {
		return 1;
	    }
	    else
	    {
		return 0;
	    }
	    
	}
	else
	{
	    die("Cannot parse $a");
	}
    }
    else
    {
	return 0;
    }

}



sub sanity_check_args
{
    my ($type, $I, $F5_start, $F5_end, $var_end, $F3_start, $F3_end) = @_;

    if (scalar (@_) != 7)
    {
	die("Calling sanity_check_args with not 7 arguments");
    }
    if ( (!defined $type)|| ($type eq "") )
    {
	die("Passed in empty type");
    }
    #sanity
    if ($F5_end != $I-1)
    {
	die("Seem to have 5prime not ending just before $type,. Ends at $F5_end, not one before $I.");
    }
    elsif ($F5_start>=$F5_end)
    {
	print "5prime flank is empty in $type. Starts at ".$F5_start." and ends at ".$F5_end;
	die();
    }
    elsif ($F5_end!=$i-1)
    {
	print "5prime flank not ending just before $type.. 5prime ends at ".$F5_end. " and SNP is at $i\n";
	die();
    }
    elsif($var_end+1>=$F3_start-1) 
    {
	print"$type at $i is after 3prime flank at ".$F3_start."\n";
	die();
    }
    elsif ($F3_start>=$F3_end)
    {
	print "$type has Null 3p flank. Starts at ".$F3_start." and ends at ".$F3_end."\n";
	die();
    }
    
    return;
}


sub get_reverse_complement_of_clean_sequence
{
    my $input_seq=$_[0];

    if (length($input_seq)<=0)
    {
        die("empty string in get_reverse_complement_of_clean_sequence");
    }

    my $reverse_seq= reverse $input_seq;
    my @reverse_seq_as_array = split(//,$reverse_seq);


    my @reverse_complement_as_array=();

    foreach my $character (@reverse_seq_as_array)
    {
	if (!exists $revcomp_hash{$character})
	{
            print "Canot rev comp this character $character\n";
            print "size of revcomp hash is ";
            print scalar(keys %revcomp_hash);
	}
        push @reverse_complement_as_array, $revcomp_hash{$character};
    }


    my $reverse_complement=join('',@reverse_complement_as_array);
    return $reverse_complement;
}



