#!/usr/bin/perl -w 
use strict;

## Given name of directory containing a bunch of split flank fasta files + the .out files resulting from exonerate via LSF.
## Will work through all the .out files, and collect a big hash of the variants that are successful - ie for which >= one flank maps uniquely.


my $dir = shift;

if ($dir !~ /\/$/)
{
    $dir = $dir."/";
}


my $get_filenames_cmd = "ls ".$dir."\*.out";
my $files = qx{$get_filenames_cmd};
my @files = split(/\n/,$files);

my %threep_hash=();
my %fivep_hash=();
my %hash_of_vars_where_5p_flank_does_not_map=();

foreach my $f (@files)
{

    open(FILE,$f)||die("cannot open $f");

    print "LOOKING AT $f\n";
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

		    if (exists($hash_of_vars_where_5p_flank_does_not_map{$varid}))
			{
			    hash_of_vars_where_5p_flank_does_not_map{$varid}=0;
			}
		}
		
		else
		{
		    die("Flank is not 3p or 5p - $which_flank");
		}
	    }
	    else
	    {
		## if this is a 5p flank, this is an exmaple of it not mapping. Mark it, and will unmark it as soon as we find an example of it mapping
		if ($which_flank eq "5p_flank")
		{
		    if (!exists $fivep_hash{$varid})
		    {
			$hash_of_vars_where_5p_flank_does_not_map{$varid}=1;
		    }
		}
	    }
	    
	}
    }
	
	
    print "These vars have unique 5p and 3p anchors\n";

#    foreach my $k (keys %threep_hash)
#    {
#	if ($threep_hash{$k}==1)
#	{
#	    if (exists $fivep_hash{$k})
#	    {
#		if ($fivep_hash{$k}==1)
#		{
#		    print "$k\n";
#		    
#		    ##set value in hash so we don't count it again
#		    $threep_hash{$k}=0;
#		    $fivep_hash{$k}=0;
#		}
#	    }
#	}
#    }

#    print "These vars have uniqe 3p flanks but not 5p\n";
#    foreach my $k (keys %threep_hash)
#    {
#	if ($threep_hash{$k}==1)
#	{
#	    ## Then it can't also have uniqe 5p flank,as we have counted them already
#	    print "$k\n";
#	    $threep_hash{$k}=0;
#	}
#    }

#    print "These vars have unique 5p but not 3p flanks\n";
#    foreach my $k (keys %fivep_hash)
#       {
#	   if ($fivep_hash{$k}==1)
#	   {
#	       ## Then it can't also have uniqe 3pflank,as we have counted them already
#	       print "$k\n";
#	       $fivep_hash{$k}=0;
#	   }
#       }
	    
    print "These vars have 5p flank that does not map!\n";
    foreach my $k (keys %hash_of_vars_where_5p_flank_does_not_map)
    {
	if ($hash_of_vars_where_5p_flank_does_not_map{$k}==1)
	{
             
	    print "$k\n";
	    $hash_of_vars_where_5p_flank_does_not_map{$k}=0;
	}
    }



}
