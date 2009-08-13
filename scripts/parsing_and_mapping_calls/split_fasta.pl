#!/usr/bin/perl -w
use strict;


## Split a fasta file into N smaller fasta files.
## Mario hads a script that does something similar, but does not preserve read names.

my $file = shift;
my $number_of_subfiles = shift;
my $output_dir = shift;

open (IN, $file) || die("Cannot open $file");
if ($output_dir !~ /\/$/)
{
    $output_dir = $output_dir."/";
}



my $read_count_cmd = "grep -c \">\" $file";
my $read_count = qx{$read_count_cmd};

my $num_reads_per_subfile = $read_count/$number_of_subfiles + 1; ## except last subfile of course
my $count = 0;

my $outfile;

while (<IN>)
{
    if ($count % $num_reads_per_subfile==0)
    {
	my $num  = $count/ $num_reads_per_subfile;
	
	$outfile = $output_dir.$file.".".$num;
	open(OUT, " >".$outfile)||die("Cannot open $outfile");
    }

    my $line = $_;
    if ($line =~ /^>/)
    {
	print OUT $line;
	$line = <IN>;
	print OUT $line;
	$count++;

	if ($count % $num_reads_per_subfile==0)
	{
	    close(OUT);
	}

    }
    else
    {
	die("Bad format to file. Expected header line, got \n$line");
    }
}

close(IN);
