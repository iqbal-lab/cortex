#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Long;

my $sample = "";
&GetOptions(
        'sample|s:s'                               => \$sample,
);

if ($sample eq "")
{
    die("You must enter the sample identifier, with --sample\n");
}


make_single_and_paired_fq_lists($sample, "/ddn/projects3/mcvean_res/zam/phase2_cortex/dummyrun/samples/".$sample, "20120522.sequence.index");

print "\n\n\n\n*********\nFinished!\n\n";



# Takes an individual identifier, and a sequence.index file, and will output 3 files, all with same root (ie path and start of filename)
# these filelists are what Cortex needs to run.
# recommend this filename root should include, or be, the individual's id - so we can run parallel processes.
# Assumes you have already downloaded the fastq for this sample 

sub make_single_and_paired_fq_lists
{
    my ($indiv, $output_filenames_root, $index) = @_;

    if ($output_filenames_root !~ /$indiv/)
    {
	die("Filename root given is $output_filenames_root and for safety we mandate this match $indiv - but it doesn't");
    }
    if (!(-d $output_filenames_root))
    {
	die("$output_filenames_root does not exist");
    }
    if ($output_filenames_root !~ /\/$/)
    {
	$output_filenames_root=$output_filenames_root.'/';
    }

    my $pout1 = $output_filenames_root.$indiv."_pe1";#left mates of paired end files
    my $pout2 = $output_filenames_root.$indiv."_pe2";#right mates of paired end files
    my $sout = $output_filenames_root.$indiv."_se";

    if ( (-e $pout1) || (-e $pout2) || (-e $sout) )
    {
	#print ("WARNING - Called make_single_and_paired_fq_lists with args $indiv, $output_filenames_root and $index, and one of these : $pout1, $pout2, $sout or the 454 ones already exists");
    }
    open(INDEX, $index)||die("Cannot ppen $index in make_single_and_paired_fq_lists");
    open(POUT1, ">".$pout1)||die("Cannot open $pout1 in make_single_and_paired_fq_lists");
    open(POUT2, ">".$pout2)||die("Cannot open $pout2 in make_single_and_paired_fq_lists");
    open(SOUT, ">".$sout)||die("Cannot open $sout in make_single_and_paired_fq_lists");
    
    my %files_covered_in_index=(); ## there is redundancy in the seq.index file - paired file sare seen twice. Use this to avoid double printing
    <INDEX>; #ignore header
    while(<INDEX>)
    {
	my $line = $_;
	chomp $line;
	
	if ($line =~ /SOLID/i)
	{
	    next;
	}
	elsif ($line =~ /exome/)
	{
	    next;
	}
	
	my @sp = split(/\t/, $line);


	if ($sp[20] != 1) ## not withdrawn
	{
	    if ( ($sp[9] eq $indiv) && (!exists $files_covered_in_index{basename($sp[0])}) )
	    {

		my $file1 = $sp[0];
		my $file2 = $sp[19];
		my $oldroot = "data/$indiv/sequence_read";
		my $newroot = $output_filenames_root;
		$file1 =~ s/$oldroot/$newroot/;
		$file2 =~ s/$oldroot/$newroot/;

		if ($file2 eq "")  ## single ended file
		{
		    $files_covered_in_index{basename($file1)}=1;
		    #single ended
		    if ($sp[12] =~ /454/)
		    {
			#ignore
		    }
		    else
		    {
			if (get_file_if_not_already_there($sp[0], $newroot, $file1)==1)
			{
			    print SOUT "$file1\n";
			}
			else
			{
			    print "WARNING - when making filelists, Ignore $file1 as it does not exist on our filesystem despite using Aspera - should have been downloaded - did something go wrong?\n";
			}
		    }
		}
		elsif ($file1 =~ /(\S+)_1(\S+)fastq/)## paired end
		{
		    my $left=$1;
		    my $right = $2;
		    if (($file2 =~ /$left/) && (($right eq ".filt.") || ($right eq ".")) )
		    {		    
			if ($sp[12] =~ /454/)
			{
			}
			else
			{
			    if ( (get_file_if_not_already_there($sp[0], $newroot, $file1)==1)
				 &&
				 (get_file_if_not_already_there($sp[19], $newroot, $file2)==1)
				)
                            {
				print POUT1 "$file1\n";
				print POUT2 "$file2\n";
			    }
			    else
			    {
				print "WARNING: One or both of $file1, $file2 is not present on the filesystem, and cannot be obtained with Aspera\n";
				print "Therefore ignore these files\n";
			    }
			    
			}
		    }
		    $files_covered_in_index{basename($file1)}=1;
		    $files_covered_in_index{basename($file2)}=1;
		    
		}
		else
		{
		}
		
	    }
	}
    }
    
    
    close(INDEX);
    close(POUT1);
    close(POUT2);
    close(SOUT);
    

    ## Double check

    my $cmd1 = "wc -l $pout1";
    my $ret1 = qx{$cmd1};
    my @sp1 = split(/\s+/, $ret1);
    my $cmd2 = "wc -l $pout2";
    my $ret2 = qx{$cmd2};
    my @sp2 = split(/\s+/, $ret2);

    if ($sp1[0] != $sp2[0]) 
    {
	die("Paired end files not same length. Check $pout1, $pout2");
    }
}



### file to get should be in sequence.index format - ie start with data/
sub get_file_if_not_already_there
{
    my ($file_to_get, $where_to_put_it, $file_you_should_have_after_getting) = @_;
    
    if (!(-e $file_you_should_have_after_getting || -l $file_you_should_have_after_getting ))
    {
	my $c1 = "/home/zam/bin/aspera/connect/bin/ascp -i /home/zam/bin/aspera/connect/etc/asperaweb_id_dsa.putty -Tr -Q -l 50M -L- fasp-g1k\@fasp.1000genomes.ebi.ac.uk:vol1/ftp/$file_to_get  $where_to_put_it";
	print "This file not present on filesystem: $file_to_get, so get it with Aspera:\n";
	my $cret1  = qx{$c1};
	print "$cret1\n";
	if (!(-e $file_you_should_have_after_getting))
	{
	    return 0;
	}
	else
	{
	    return 1;
	}
    }
    else
    {
	##already exists
	return 1;
    }
}
