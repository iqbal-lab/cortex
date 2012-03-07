#!/usr/bin/perl -w
use strict;
use Getopt::Long;



### Takes a filelist of call files 
### If there is print_colour_coverages output, will ignore it. Similarly genotypes.
### Outputs a single set of calls, renamed UNION_VAR_1, UNION_VAR_2,.. etc


my ($filelist, $varname_stub);

&GetOptions(
     	    'filelist|f:s'                        => \$filelist,              #file contains one big fasta entry
            'varname_stub|s:s'                    => \$varname_stub,                  #variant names will start with this, then _var_ then number
           );


my %vars_hash=(); ## has branch1/branch2 =>1  (when checking against hash, must check b1/b2 and b2/b1 and rev comps)


### fill hash 
open(FILE, $filelist)||die();
while (<FILE>)
{
    my $line = $_;
    chomp $line;
    get_all_calls_from_this_file($line, \%vars_hash);
}
close(FILE);

print_unionset(\%vars_hash, $filelist, $varname_stub);

sub print_unionset
{
    my ($href, $fylelist, $stub) = @_;

    my $global_ct=1;

    open (FH, $fylelist)||die();
    while (<FH>)
    {
	my $lyne = $_;
	chomp $lyne;
	
	open(SOMEFILE, $lyne)||die();
	while (<SOMEFILE>)
	{
	    my $thisln = $_;
	    if ($thisln =~ /5p_flank/)
	    {
		my $vname = $stub."_var_".$global_ct;
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
			print  $f5_id."\n".$f5_seq.$b1_id."\n".$b1_seq.$b2_id."\n".$b2_seq.$f3_id."\n".$f3_seq;
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

sub get_all_calls_from_this_file
{
    my ($file, $href) = @_;
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
	    $href->{$b1."_".$b2}=1;
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

