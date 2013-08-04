#!/usr/bin/perl -w
use strict;

my $vcf_file = shift;##just check overlapping-marked lines. If they are consistent, keep them

my $vcf_handle;

if(defined($vcf_file))
{
    open($vcf_handle, $vcf_file)
	or print_usage("Cannot open VCF file '$vcf_file'\n");
}
elsif(-p STDIN) {
    # STDIN is connected to a pipe
    open($vcf_handle, "<&=STDIN") or print_usage("Cannot read pipe");
}
else
{
    print_usage("Must specify or pipe in a VCF file");
}

my $line="";
my %sites=();
my $curr_num=0;
my $curr_far_right_of_overlapping_set=-1;
while (<$vcf_handle>)
{
    $line = $_;
    if ($line =~ /^\#/)
    {
	print $line;
	next;
    }
    chomp $line;
    my @sp = split(/\s+/, $line);
    if (scalar @sp <9)
    {
	print scalar @sp;
	die("This vcf $vcf_file contains a line with <9 columns: $line\n");
    }
    my $chr = $sp[0];
    my $pos = $sp[1];
    my $name = $sp[2];
    my $ref = $sp[3];
    my $alt = $sp[4];
    my $filter = $sp[6];
    my $info   = $sp[7];
    #$curr_far_right_of_overlapping_set = $pos+length($ref);
    if ($info=~/SVTYPE=SNP/)
    {
	$info="snp";
    }
    else
    {
	$info="indel";
	$pos++;
	$ref=substr($ref,1);
	$alt=substr($alt,1);
    }
    if ($filter !~ /OVERLAPPING_SITE/) #||($curr_far_right_of_overlapping_set==-1)) ##end of this collection of overlapping sites
    {
	print_prev_overlapping_sites_if_any(\%sites);
	%sites=();
	$curr_num=0;
	$curr_far_right_of_overlapping_set=-1;
	print "$line\n";
    }
    else#it is an overlapping site. Now must allowfor a set of consec sites in vcf all marked as overlapping,
       #but which are >1 overlapping loci (ie 2 overlapping sets with no PPASS between them
    {
	if ($pos>$curr_far_right_of_overlapping_set) 
	{
	    print_prev_overlapping_sites_if_any(\%sites);
	    %sites=();
	    $curr_num=0;
	    $curr_far_right_of_overlapping_set=-1;
	}
	push @{$sites{"chr"}},$chr;
	push @{$sites{"pos"}},$pos;
	push @{$sites{"name"}},$name;
	push @{$sites{"ref"}},$ref;
	push @{$sites{"alt"}},$alt;
	push @{$sites{"line"}},$line;
	push @{$sites{"filter"}},$filter;
	push @{$sites{"type"}}, $info;
	if ($pos+length($ref)>$curr_far_right_of_overlapping_set)
	{
	    $curr_far_right_of_overlapping_set=$pos+length($ref);
	}
	
	#else
	#{
	 #   print_prev_overlapping_sites_if_any(\%sites);
	  #  %sites=();
	   # $curr_num=0;
	    #$curr_far_right_of_overlapping_set=-1;
	    #print "BIF $line\n";
	#}
    }
}
close($vcf_handle);

print_prev_overlapping_sites_if_any(\%sites);


## note we have the true position and alleles by this stage, so SNP and indel treated the same
sub print_prev_overlapping_sites_if_any
{
    my ($href) = @_;
    if (!exists $href->{"chr"})
    {
	return;
    }
    elsif (scalar @{$href->{"chr"}} ==1)
    {
	print "There is a solitary overlapping site:\nZAM ";
	print ${$href->{"line"}}[0];
	die();
    }
    else
    {
	my $i;
	my @passes = ();
	my $curr_ref = ${$href->{"ref"}}[0];
	my $curr_alt = ${$href->{"alt"}}[0];
	my $curr_chr = ${$href->{"chr"}}[0];
	my $curr_pos = ${$href->{"pos"}}[0];
	my $curr_line= ${$href->{"line"}}[0];
	my $curr_filter =  ${$href->{"filter"}}[0];
	my $curr_index=0;
	my $pass_index=-1;

	for ($i=1; $i< scalar @{$href->{"chr"}}; $i++)
	{
	    my $ref2 = ${$href->{"ref"}}[$i];
	    my $alt2 = ${$href->{"alt"}}[$i];
	    my $chr2 = ${$href->{"chr"}}[$i];
	    my $pos2 = ${$href->{"pos"}}[$i];
	    my $line2 = ${$href->{"line"}}[$i];
	    my $filter2 = ${$href->{"filter"}}[$i];
	    my $name2 = ${$href->{"name"}}[$i];

	    if ($filter2 ne "OVERLAPPING_SITE") #then it must fail some other filter
	    {
		next;
	    }
	    if  (  ($curr_chr eq $chr2) && ($curr_pos==$pos2) )
	    {
		#start at same place
		if (length($curr_ref) < length($ref2))
		{
		    if ( ($ref2 =~ /^$curr_ref/) && ($alt2 =~ /^$curr_alt/))
		    {

			## set the longer one to PASS - new one is the only one to PASS
			$curr_ref=$ref2;
			$curr_alt=$alt2;
			$curr_chr=$chr2;
			$curr_pos=$pos2;
			$curr_line=$line2;
			$pass_index=$i;
			$curr_index=$i;
		    }
		    else
		    {
		    }
		    
		}
		else
		{
		    if ( ($curr_ref =~ /^$ref2/) && ($curr_alt =~ /^$alt2/))
		    {
			## set the longer one to PASS - ie current
			$pass_index=$curr_index;
		    }
		    else 
		    {
		    }

		}
	    }
	    elsif ( ($curr_chr eq $chr2)
		    &&
		    ($curr_pos+length($curr_ref)==$pos2+length($ref2) )
		)
		    
		{
		    #they meet at the end
		    if (length($curr_ref) < length($ref2) )
		    {
			if  ($ref2 =~ /$curr_ref$/)  
			{
			    ## set the longer one to PASS
			    $curr_ref=$ref2;
			    $curr_alt=$alt2;
			    $curr_chr=$chr2;
			    $curr_pos=$pos2;
			    $curr_line=$line2;
			    $pass_index=$i;
			    $curr_index=$i;
			}
			else
			{
			}
		    }
		    else
		    {
			if ( ($curr_ref =~ /$ref2$/) )
			{
			## set the longer one to PASS
			    $pass_index=$curr_index;
			}
			else
			{
			}
		    }
		}
		else
		{
		}
		
	}
	if ($pass_index!=-1)
	{
	    ${$href->{"line"}}[$pass_index] =~ s/OVERLAPPING_SITE/PASS/;
	}

	my $j;
	for ($j=0; $j< scalar @{$href->{"chr"}}; $j++)
	{
	    print ${$href->{"line"}}[$j];
	    print "\n";
	}
    }
}

