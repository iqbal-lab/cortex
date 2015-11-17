#!/usr/bin/perl -w
use strict;

use Getopt::Long;

my $truth_vcf="";
my $test_vcf="";

&GetOptions(
    ##mandatory args
    'truth_vcf:s' => \$truth_vcf,
    'test_vcf:s'  => \$test_vcf,
    );


### Collect true genotypes
my %true_gt=(); #chr_pos -> gt
get_truth(\%true_gt, $truth_vcf);

## prepare to Collect stats
my %conf_bins=(); ## bin -> counts(TrueRR_calledRR, True RR called het, True RR called AA, True RR called ./., etc)
my @calls=("TrueRR_calledRR",
	   "TrueRR_calledAR",
	   "TrueRR_calledAA",
	   "TrueRR_calledmissing",
	   "TrueAR_calledRR",
	   "TrueAR_calledAR",
	   "TrueAR_calledAA",
	   "TrueAR_calledmissing",
	   "TrueAA_calledRR",
	   "TrueAA_calledAR",
	   "TrueAA_calledAA",
	   "TrueAA_calledmissing",
	   "Truemissing_calledRR",
	   "Truemissing_calledAR",
	   "Truemissing_calledAA",
	   "Truemissing_calledmissing");
my $i;
for ($i=0; $i<=100; $i+=1)
{
    foreach my $call (@calls)
    {
	$conf_bins{$i}{$call}=0;
    }
}

## Collect stats
get_data(\%conf_bins, $test_vcf, \%true_gt);

## print stats

my $j;
my $inc = 10;
for ($j=0; $j<=90; $j+=$inc)
{
    my $start = $j;
    my $end = $j+$inc-1;
    print_aggregates(\%conf_bins, \@calls, $start, $end);
}

sub print_aggregates
{
    my ($href_stats, $aref_calls, $s, $e)=@_;
    my %results=();
    foreach my $c (@$aref_calls)
    {
	$results{$c}=0;
    }
    my $k;
    for ($k=$s; $k<=$e; $k++)
    {
	foreach my $c (@$aref_calls)
	{
	    $results{$c} += $href_stats->{$k}->{$c};
	}
    }
    print "\n\nResults for Confidences between $s and $e\n";
    print "\tTrueRR\tTrueAR\tTrueAA\tTrueMissing\n";
    print "CallRR\t";
    print $results{"TrueRR_calledRR"}."\t";
    print $results{"TrueAR_calledRR"}."\t";
    print $results{"TrueAA_calledRR"}."\t";
    print $results{"Truemissing_calledRR"}."\n";
    print "CallAR\t";
    print $results{"TrueRR_calledAR"}."\t";
    print $results{"TrueAR_calledAR"}."\t";
    print $results{"TrueAA_calledAR"}."\t";
    print $results{"Truemissing_calledAR"}."\n";
    print "CallAA\t";
    print $results{"TrueRR_calledAA"}."\t";
    print $results{"TrueAR_calledAA"}."\t";
    print $results{"TrueAA_calledAA"}."\t";
    print $results{"Truemissing_calledAA"}."\n";
    print "Callmissing\t";
    print $results{"TrueRR_calledmissing"}."\t";
    print $results{"TrueAR_calledmissing"}."\t";
    print $results{"TrueAA_calledmissing"}."\t";
    print $results{"Truemissing_calledmissing"}."\n";


}




sub get_data
{
    my ($href_save, $vcf, $href_truth) = @_;

    open(VCF, $vcf)||die();
    my $which_field_is_conf=-1;
    my $which_field_is_gt=-1;
    while (<VCF>)
    {
	my $line = $_;
	chomp $line;
	if ($line =~ /^\#/)
	{
	}
	else
	{
	    my @sp = split(/\t/, $line);
	    my $id = $sp[0]."_".$sp[1];
	    if (! exists $href_truth->{$id})
	    {
		next;
	    }
	    my $format = $sp[8];
	    $which_field_is_conf = get_gt_conf_field_index($format);
	    $which_field_is_gt = get_gt_field_index($format);
	    
	    my @fields = split(/:/, $sp[9]);
	    my $conf = $fields[$which_field_is_conf];
	    my $gt = $fields[$which_field_is_gt];
	    my $rounded_conf = 0;
	    if ( ($conf>=1) & ($conf<=100) )
	    {
		$rounded_conf = int($conf);
	    }
	    elsif ($conf>=1)
	    {
		$rounded_conf = 100;
	    }
	    
	    if ( ($gt eq "0/0") && ($href_truth->{$id} eq "0/0") )
	    {
		$href_save->{$rounded_conf}->{"TrueRR_calledRR"} +=1; 
	    }
	    elsif ( ($gt eq "0/0") && ($href_truth->{$id} eq "0/1") )
	    {
		$href_save->{$rounded_conf}->{"TrueAR_calledRR"} +=1; 
	    }
	    elsif ( ($gt eq "0/0") && ($href_truth->{$id} eq "1/1") )
	    {
		$href_save->{$rounded_conf}->{"TrueAA_calledRR"} +=1; 
	    }
	    elsif ( ($gt eq "0/0") && ($href_truth->{$id} eq "./.") )
	    {
		$href_save->{$rounded_conf}->{"Truemissing_calledRR"} +=1; 
	    }

	    elsif ( ($gt eq "0/1") && ($href_truth->{$id} eq "0/0") )
	    {
		$href_save->{$rounded_conf}->{"TrueRR_calledAR"} +=1; 
	    }
	    elsif ( ($gt eq "0/1") && ($href_truth->{$id} eq "0/1") )
	    {
		$href_save->{$rounded_conf}->{"TrueAR_calledAR"} +=1; 
	    }
	    elsif ( ($gt eq "0/1") && ($href_truth->{$id} eq "1/1") )
	    {
		$href_save->{$rounded_conf}->{"TrueAA_calledAR"} +=1; 
	    }
	    elsif ( ($gt eq "0/1") && ($href_truth->{$id} eq "./.") )
	    {
		$href_save->{$rounded_conf}->{"Truemissong_calledAR"} +=1; 
	    }

	    elsif ( ($gt eq "1/1") && ($href_truth->{$id} eq "0/0") )
	    {
		$href_save->{$rounded_conf}->{"TrueRR_calledAA"} +=1; 
	    }
	    elsif ( ($gt eq "1/1") && ($href_truth->{$id} eq "0/1") )
	    {
		$href_save->{$rounded_conf}->{"TrueAR_calledAA"} +=1; 
	    }
	    elsif ( ($gt eq "1/1") && ($href_truth->{$id} eq "1/1") )
	    {
		$href_save->{$rounded_conf}->{"TrueAA_calledAA"} +=1; 
	    }
	    elsif ( ($gt eq "1/1") && ($href_truth->{$id} eq "./.") )
	    {
		$href_save->{$rounded_conf}->{"Truemissing_calledAA"} +=1; 
	    }

	    elsif ( ($gt eq "./.") && ($href_truth->{$id} eq "0/0") )
	    {
		$href_save->{$rounded_conf}->{"TrueRR_calledmissing"} +=1; 
	    }
	    elsif ( ($gt eq "./.") && ($href_truth->{$id} eq "0/1") )
	    {
		$href_save->{$rounded_conf}->{"TrueAR_calledmissing"} +=1; 
	    }
	    elsif ( ($gt eq "./.") && ($href_truth->{$id} eq "1/1") )
	    {
		$href_save->{$rounded_conf}->{"TrueAA_calledmissing"} +=1; 
	    }
	    elsif ( ($gt eq "./.") && ($href_truth->{$id} eq "./.") )
	    {
		$href_save->{$rounded_conf}->{"Truemissing_calledmissing"} +=1; 
	    }
	    else
	    {
		print "\n\nProblem in line $line\n Gt is $gt\n";
		print "Truth is ";
		print $href_truth->{$id};
		print "\n";
		die("dead");
	    }

	}
    }
    close(VCF);
}



sub get_truth
{
    my ($href, $vcf) = @_;
    open(FILE, $vcf)||die();
    my $which_field_is_gt=-1;
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
	    my $id = $sp[0]."_".$sp[1];
	    my $format = $sp[8];
	    $which_field_is_gt = get_gt_field_index($format);
	    my @fields = split(/:/, $sp[9]);
            my $gt = $fields[$which_field_is_gt];
	    if ($gt eq "1/0")
	    {
		$gt = "0/1";
	    }
	    $href->{$id}=$gt;
	}
    }
    close(FILE);
}


sub get_gt_field_index
{
    my ($str) = @_;
    my @spl = split(/:/, $str);
    my $i;
    for ($i=0; $i<scalar(@spl); $i++)
    {
	if ($spl[$i] eq "GT")
	{
	    return $i;
	}
    }
    die("Unable to find GT in $str\n");
}

   
sub get_gt_conf_field_index
{
    my ($str) = @_;
    my @spl = split(/:/, $str);
    my $i;
    for ($i=0; $i<scalar(@spl); $i++)
    {
	if ($spl[$i] eq "GT_CONF")
	{
	    return $i;
	}
    }
    die("Unable to find GT_CONF in $str\n");
}
