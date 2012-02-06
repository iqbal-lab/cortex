#This script takes a chromosome 
#and output from ms, and uses the haplotypes ms has output to generate a set of chromosomes.


use strict;

use Getopt::Long;


my ($chromosome_fasta, $num_sites, $num_chroms,  $kmer_size, $samples);

&GetOptions(
            'chr|c:s'                => \$chromosome_fasta,
	    'num_sites|s:i'         => \$num_sites,      # number of variant sites
	    'num_chroms|n:i'         => \$num_chroms,      # number of haplotypes/chromosomes
	    'kmer_size|k:i'          => \$kmer_size,       # this is only used when printing the files of variants (flanks + branches) 
           );

my $space = 2*$kmer_size;



my ($href_ns, $chromosome_length)=get_coords_of_ns($chromosome_fasta);

my %mask_ref=(); ## coords of positions of SNPs---> ref/alt alleles
my %mask_alt=();

my %hash_site_number_to_freq=(); ## which-th site, counting from start of chrom, ---> number of haplotypes carrying ALT allele
pick_allele_freqs_according_to_sfs($num_chroms, $num_sites,\%hash_site_number_to_freq);

my %hash_site_to_which_haps_have_which_allele=(); ## site --> hap number --> 0 (ref) or 1(alt)
pick_which_haplotypes_have_alt_alleles($num_chroms, $num_sites,\%hash_site_number_to_freq, \%hash_site_to_which_haps_have_which_allele);

my $href_sites = allocate_positions_and_alleles(\%mask_ref, \%mask_alt, $chromosome_length, $num_sites, $space, $href_ns);

my $j;
for ($j=1; $j<=$num_chroms; $j++)
{
    print_fasta($num_sites, $href_sites, \%mask_ref, \%mask_alt, $j, $chromosome_fasta, $num_chroms, \%hash_site_to_which_haps_have_which_allele);    
}


print_coords_and_alleles($href_sites, \%mask_ref, \%mask_alt, $chromosome_length, \%hash_site_number_to_freq, \%hash_site_to_which_haps_have_which_allele, $num_chroms);



sub print_coords_and_alleles
{
    my ($href_sits, $href_ref, $href_alt, $len, $href_site_to_number_alt, $href_site_to_which_haps_have_which, $num_haps) = @_;

    my $out = "sites_and_alleles";
    open(OUT, ">".$out)||die("Cannot open $out");
    print OUT "Pos\tRef\tAlt\tNum_alt\tAlleles_in_haplotypes\n";
    my $i;
    my $count_sites=0;
    for ($i=1; $i<=$len; $i++)
    {
	if (exists $href_sits->{$i})
	{
	    $count_sites++;
	    print OUT "$i\t";
	    print OUT $href_ref->{$i};
	    print OUT "\t";
	    print OUT $href_alt->{$i};
	    print OUT "\t";
	    print OUT $href_site_to_number_alt->{$count_sites};
	    print OUT "\t";
	    my $j;
	    for ($j=1; $j<=$num_haps; $j++)
	    {
		print OUT $href_site_to_which_haps_have_which->{$count_sites}->{$j}; #0 for ref and 1 for alt
	    }
	    print OUT "\n";
	}
    }
    close(OUT);
}



sub get_coords_of_ns
{
    my ($fasta) = @_;
    open (FA, $fasta)||die("Cannot open $fasta");
    my $coord=0;

    my %hash=();

    while (<FA>)
    {
        my $line = $_;
        chomp $line;
        if ($line =~ /^>/)
        {
            next;
        }
        my @sp = split(//, $line);
        my $i;
        for ($i=0; $i<scalar(@sp); $i++)
        {
            $coord++;
            if ($sp[$i] !~ /[ACGT]/)
            {
                $hash{$coord}=1;
            }
        }
    }
    close(FA);
    
    return (\%hash, $coord);
}


sub pick_which_haplotypes_have_alt_alleles
{

    my ($num_haps, $num_sites, $href_site_number_to_freq, $href_site_to_which_haps_have_which_allele) = @_;


    my $i;
    for ($i=1; $i<=$num_sites; $i++)
    {
	my $number_alt = $href_site_number_to_freq->{$i};
	my %which_are_alt=();
	while (scalar(keys %which_are_alt) < $number_alt)
	{
	    my $n = int(rand($num_haps+1));## returs integer between 0 and $num_haps
	    if ($n>0)
	    {
		$which_are_alt{$n}=1;
	    }
	}
	my $j;
	for ($j=1; $j<=$num_haps; $j++)
	{
	    if (exists $which_are_alt{$j})
	    {
		$href_site_to_which_haps_have_which_allele->{$i}->{$j}=1;
	    }
	    else
	    {
		$href_site_to_which_haps_have_which_allele->{$i}->{$j}=0;
	    }
	}
    
    }


}

sub pick_allele_freqs_according_to_sfs
{
    my ($number_haps, $number_sites, $href_frequencies) = @_;

    my %freq_to_prob=(); ## frequency (number of haplotypes with the alt allele) --> probability of that freq
    my $i;
    my $sum=0;
    my @partition=();
    push @partition,0;

    for ($i=1; $i<$number_haps; $i++)
    {
	my $tmp = $i/$number_haps;
	$freq_to_prob{$i} = 1/($tmp*(1-$tmp));
	$sum += $freq_to_prob{$i};
    }
    for ($i=1; $i<$number_haps; $i++)
    {
	$freq_to_prob{$i} = $freq_to_prob{$i}/$sum;
	push @partition, $partition[$i-1] + $freq_to_prob{$i};
    }

    print "These are the prior probs fo freqs, and the partition\n";
    my $m;
    for ($m=1; $m<$number_haps; $m++)
    {
	print "$m\t";
	print $freq_to_prob{$m};
	print "\t";
	print $partition[$m];
	print "\n";
    }

    for ($i=1; $i<=$number_sites; $i++)
    {
	my $p=0;
	while ( ($p==0) || ($p==1) )
	{
	    $p = rand();
	}
	my $j;
	my $freq=-1;
	for ($j=1; $j<$number_haps; $j++)
	{
	    if ($p<=$partition[$j])
	    {
		$freq = $j;
		last;
	    }
	}
	if ($freq==-1)
	{
	    die("CODING ERROR");
	}
	$href_frequencies->{$i}=$freq;
    }
    

}
sub print_fasta
{
    my  ($number_sites, $href_coords, $href_ref, $href_alt, $number_of_hap, $chr_fa, $number_chromosomes, $href_site_to_num_alt_alleles) = @_;
    open(OUT, "> haplotype_".$number_of_hap)||die("Cannot opne output file");
    open(FASTA, $chr_fa)||die("Cannot open $chr_fa");
    

    my $coord=0;
    my $count_sites=-1;
    while (<FASTA>)
    {
	my $line = $_;
	chomp $line;
	
        if ($line =~ /^>/)
        {
	    print OUT "$line\n";
            next;
	}
        my @sp = split(//, $line);
        my $i;
        for ($i=0; $i<scalar(@sp); $i++)
        {
            $coord++;
	    if (!exists $href_coords->{$coord})
	    {
		print OUT $sp[$i];
	    }
	    else
	    {
		$count_sites++;
		my $ref_allele = $sp[$i];
		my $alt_allele;
		if (!exists $href_ref->{$coord})
		{
		    $alt_allele = get_alt_allele($ref_allele);
		    $href_ref->{$coord}=$ref_allele;
		    $href_alt->{$coord}=$alt_allele;
		}
		else
		{
		    $alt_allele = $href_alt->{$coord};
		    
		    ##double check
		    if ($href_ref->{$coord} ne $ref_allele)
		    {
			die("Some coord error at $coord, ref allele is $ref_allele, but we think otherwise");
		    }
		}

		## print out what the haplptype array tells us to 0 is ancestral (ref)
		if ($href_site_to_num_alt_alleles->{$count_sites}->{$number_of_hap}==0)
		{
		    print OUT $ref_allele;
		}
		else
		{
		    print OUT $alt_allele;
		}
	    }
	}
	print OUT "\n";

    }
    close(OUT);
    close(FASTA);
    
}


sub get_alt_allele
{
    my ($ref) = @_;
    my $which = int(rand(3)); # 0,1, or 2
    
    my @arr=("A","C","G","T");
    my @others=();
    my $i;
    for ($i=0; $i<4; $i++)
    {
	if ($ref eq $arr[$i])
	{
	}
	else
	{
	    push @others, $arr[$i];
	}
    }

    return $others[$which];
}


sub allocate_positions_and_alleles
{
    my ($href_ref, $href_alt, $chr_len, $num_sites, $space, $href_coords_of_ns)=@_;

    
    ### select $num_sites coords from 1 to $chr_len
    my %sites=();
    my $count=0;
    while ($count<$num_sites)
    {
	my $coord = int(rand($chr_len));
	if (check_space($coord, \%sites, $space, $href_coords_of_ns)==1)
	{
	    $sites{$coord}=1; 
	    $count++;
	}
    }

    return \%sites;
    
}

sub check_space
{
    my ($pos, $href, $space_to_leave, $href_coords_of_Ns) = @_;
    
    my $i;
    my $flag=1;
    for ($i=-$space_to_leave; $i<=$space_to_leave; $i++)
    {
	if ( (exists $href->{$pos+$i}) || (exists $href_coords_of_Ns->{$pos+$i}) )
	{
	    $flag=0;
	    last;
	}
    }
    return $flag;
}
