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
    my $ref = $sp[3];
    my $alt = $sp[4];
    my $filter = $sp[6];
    my $info   = $sp[7];
    if ($filter ne "OVERLAPPING_SITE") #then it is not overlaping, or it is overlapping, but fails some other filter also 
    {
	print_prev_overlapping_sites_if_any(\%sites);
	%sites=();
	$curr_num=0;
	print "$line\n";
    }
    else#it is an overlapping site, and fails no other filter
    {
	push @{$sites{"chr"}},$chr;
	push @{$sites{"pos"}},$pos;
	push @{$sites{"ref"}},$ref;
	push @{$sites{"alt"}},$alt;
	push @{$sites{"line"}},$line;
	push @{$sites{"filter"}},$filter;
	
    }
}
close($vcf_handle);

print_prev_overlapping_sites_if_any(\%sites);

sub print_prev_overlapping_sites_if_any
{
    my ($href) = @_;
    if (!exists $href->{"chr"})
    {
	return;
    }
    elsif (scalar @{$href->{"chr"}} ==1)
    {
	print "There is a solitary overlapping site:\n";
	print ${$href->{"line"}}[0];
	die();
    }
    
#    if (scalar @{$href->{"chr"}} >2)
#    {
#	my $i;
#	for ($i=0; $i< scalar @{$href->{"chr"}}; $i++)
#	{
#	    print ${$href->{"line"}}[$i];
#	    print "\n";
#	}
#   }
    else
    {
	my $i;
	my $to_print = "";
	for ($i=1; $i< scalar @{$href->{"chr"}}; $i++)
	{
	    my $ref1 = ${$href->{"ref"}}[$i-1];
	    my $ref2 = ${$href->{"ref"}}[$i];
	    my $alt1 = ${$href->{"alt"}}[$i-1];
	    my $alt2 = ${$href->{"alt"}}[$i];
	    my $chr1 = ${$href->{"chr"}}[$i-1];
	    my $chr2 = ${$href->{"chr"}}[$i];
	    my $pos1 = ${$href->{"pos"}}[$i-1];
	    my $pos2 = ${$href->{"pos"}}[$i];
	    my $line1 = ${$href->{"line"}}[$i-1];
	    my $line2 = ${$href->{"line"}}[$i];

	if  (  ($chr1 eq $chr2)
	       &&
	       ($pos1==$pos2) )
	{
	    #start at same place
	    if (length($ref1) < length($ref2) )
	    {
		if ( ($ref2 =~ /^$ref1/) && ($alt2 =~ /^$alt1/))
		{
		    ## set the longer one to PASS
		    my $newline = $line2;
		    $newline =~ s/OVERLAPPING_SITE/PASS/;
		    print "$line1\n$newline\n";
		}
		
	    }
	    elsif (length($ref1) >= length($ref2) )
	    {
		if ( ($ref1 =~ /^$ref2/) && ($alt1 =~ /^$alt2/))
		{
		    ## set the longer one to PASS
		    my $newline = $line1;
		    $newline =~ s/OVERLAPPING_SITE/PASS/;
		    print "$newline\n$line2\n";
		}
		
	    }
	}
	elsif ( ($chr1 eq $chr2)
		&&
		($pos1+length($ref1) == $pos2+length($ref2)) )

	{
	    #they meet at the end
	    if (length($ref1) < length($ref2) )
	    {
		if ( ($ref2 =~ /$ref1$/) && ($alt2 =~ /$alt1$/) )
		{
		    ## set the longer one to PASS
		    my $newline = $line2;
		    $newline =~ s/OVERLAPPING_SITE/PASS/;
		    print "$line1\n$newline\n";
		}
		
	    }
	    elsif (length($ref1) >= length($ref2) )
	    {
		if ( ($ref1 =~ /$ref2$/) && ($alt1 =~ /$alt2$/))
		{
		    ## set the longer one to PASS
		    my $newline = $line1;
		    $newline =~ s/OVERLAPPING_SITE/PASS/;
		    print "$newline\n$line2\n";
		}
		
	    }
	    
	}
	    else
	    {
	    #print them as they are
	    print "$line1\n$line2\n";
	}

    }
}
