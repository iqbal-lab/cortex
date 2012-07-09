#!/usr/bin/perl -w
use strict;
use Getopt::Long;



### Takes an index file with callfile \t kmer \t cleaning
### If there is print_colour_coverages output, will ignore it. Similarly genotypes.
### Outputs a single set of calls, renamed UNION_kmer14-cleaning23_VAR_1, UNION_kmer_16_cleaing24_VAR_2,.. etc. ie var nums just increment. stuff before 


my ($index, $varname_stub, $outfile, $kmer);

&GetOptions(
     	    'index|f:s'                        => \$index,#each row is filename \t kmer \t cleaning

            'varname_stub|s:s'                 => \$varname_stub,                  
                                                     # variant names will start with this, 
                                                     # then _var_ then number
            'outfile|o:s'                      => \$outfile,
            'kmer|k:s'                         => \$kmer
           );


open(OUT, ">".$outfile)||die();
my %vars_hash=(); ## has branch1/branch2 =>1  
                  ## (when checking against hash, must check b1/b2 and b2/b1 and rev comps)
my %var_to_kmer=();

### fill hash 
open(FILE, $index)||die();
while (<FILE>)
{
    my $line = $_;
    chomp $line;
    my @sp = split(/\t/, $line);
    if (scalar(@sp) !=3)
    {
	die("Index file $index is supposed to have 3 columns, filename (of calls) , kmer and cleaning\n");
    }
    my $fname = $sp[0];
    my $km = $sp[1];
    if ($km != $kmer)
    {
	die("Calling make_union_varset_at_single_kmer.pl with an index containing a  line referring to a kmer of $km, but the --kmer option was $kmer\n");
    }
    my $cleaning = $sp[2];
    get_all_calls_from_this_file($fname, $kmer, $cleaning, \%vars_hash);
}
close(FILE);

print_unionset(\%vars_hash, $index, $varname_stub);
close(OUT);
#close(KMER_INDEX);

sub print_unionset
{
    my ($href, $index, $stub) = @_;

    my $global_ct=1;
    open (FH, $index)||die("Cannot open $index");
    while (<FH>)
    {
	my $index_line = $_;
	chomp $index_line;
	my @sp = split(/\t/, $index_line);
	my $f = $sp[0];
	my $km = $sp[1];
	my $cle = $sp[2];

	open(SOMEFILE, $f)||die();
	while (<SOMEFILE>)
	{
	    my $thisln = $_;
	    if ($thisln =~ /5p_flank/)
	    {
		my $vname = $stub."_var_".$global_ct;
		$var_to_kmer{$vname} = $km;
		my $f5_id = ">".$vname."_5p_flank";
		my $f5_seq=<SOMEFILE>;
		<SOMEFILE>;
		my $b1_id=">$vname"."_branch_1";
		my $b1_seq=<SOMEFILE>;
		<SOMEFILE>;
		my $b2_id=">$vname"."_branch_2";
		my $b2_seq=<SOMEFILE>;
		<SOMEFILE>;
		my $f3_id=">".$vname."_3p_flank";
		my $f3_seq=<SOMEFILE>;

		my $t1=$b1_seq;
		chomp $t1;
		my $t2=$b2_seq;
		chomp $t2;
		
		if (exists $href->{$t1."_".$t2})
		{
		    if ($href->{$t1."_".$t2} ne "0")
		    {
			my ($smallest_k, $cleaning_for_smallest_k) = get_smallest_k($vname, $href->{$t1."_".$t2});
			my $info =  "INFO:KMER=$smallest_k";
			print  OUT $f5_id."\t$info\n".$f5_seq.$b1_id."\n".$b1_seq.$b2_id."\n".$b2_seq.$f3_id."\n".$f3_seq;
			$href->{$t1."_".$t2}=0;
			$global_ct++;
		    }
		}
	    }
	    }
	close(SOMEFILE);

	
    }
    close(FH);
    
}

sub get_smallest_k
{
    my ($name, $str) = @_; ## format 31,12=41,4=51,5  which are kmer,cleaning  with = as separator
    my @sp = split(/_/, $str);
    my $i;
    my $min = 9999999999;
    my $clean_that_goes_with_min=999;
    for ($i=0; $i<scalar(@sp); $i++)
    {
	my @sp2 = split(/,/,$sp[$i]);
	if ($sp2[0]<$min)
	{
	    $min = $sp2[0];
	    $clean_that_goes_with_min = $sp2[1];
	}
	elsif (($sp2[0]==$min) && ($sp2[0]<$clean_that_goes_with_min) )## of multiple options with the same kmer, choose the one with least cleaning
	{
	    $min = $sp2[0];
	    $clean_that_goes_with_min = $sp2[1];
	}
    }
    if (scalar(@sp)>1)
    {
	print "This variant: $name was called with the same kmer with multiple cleanings,\n";
	print " specifically with kmer,cleaning =";
	print join(" and ", @sp);
	print "\n";
	print "So we choose this representative: kmer = $min, cleaning = $clean_that_goes_with_min\n";
    }
    return ($min, $clean_that_goes_with_min);
}


sub get_all_calls_from_this_file
{
    my ($file, $k, $cl, $href) = @_;
    open(F, $file)||die("Cannot find file $file");
    while (<F>)
    {
	my $ln = $_;
	chomp $ln;
	if ($ln =~ /5p_flank/)
	{
	    <F>; #ignore 5p
	    <F>; #ignore b1 read id
	    my $b1 = <F>;
	    chomp $b1;
	    <F>; #ignore b2 read id
	    my $b2 = <F>;
	    chomp $b2;
	    if (!exists $href->{$b1."_".$b2})
	    {
		$href->{$b1."_".$b2}=$k.",".$cl;
	    }
	    else
	    {
		$href->{$b1."_".$b2}=$href->{$b1."_".$b2}."_".$k.",".$cl;
	    }
	}
    }
    close(F);
}
    

sub checkhash
{
    my ($str, $href) = @_;
    my @sp = split(/_/, $str);
    if (scalar @sp != 2)
    {
	die("Cant split $str");
    }
    
    if (exists $href->{$str})
    {
	return 1;
    }
    elsif  (exists $href->{$sp[1]."_".$sp[0]})
    {
	return 1;
    }
    else
    {
	return 0;
    }
}




sub rev_comp
{
  my ($seq) = @_;

  my $r_seq = reverse($seq);
  $r_seq =~ tr/acgtACGT/tgcaTGCA/;
 # print join(" ",$seq,$r_seq),"\n";
  return $r_seq;

}

