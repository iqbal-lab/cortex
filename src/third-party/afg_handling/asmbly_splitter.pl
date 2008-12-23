#!/usr/bin/perl -w
#
#   asmbly_splitter.pl
#
#	A script to split out a single contig's data from an AFG file
#	produced by Daniel Zerbino's Velvet assembler.  The output is 
#	in AFG format and is a subset of the data from the original file.
#	Because of the size of the AFG file, this script can take quite a
#	while to run.
#
#	Simon Gladman 2008: CSIRO Australia.
#
#	Usage:  ./asmbly_splitter.pl <contig number> <afg file>
#
#	Where: 	<contig number> is the number of the contig of interest.
#			<afg file> is the AFG file produced by Velvet.
#

use strict;

my $usage = "./asmbly_splitter.pl <contig number> <afg file>\n";

sub toggle {
	
	my $x = shift;
	if($x == 1){
		$x = 0;
	}
	else {
		$x = 1;
	}
	
	return $x;
}

my $contig = $ARGV[0];
my $file = $ARGV[1];

unless($contig){die "Usage: $usage \n";}

unless($file){die "Usage: $usage \n";}

unless(-e $file){die "$0: File $file doesn't exist.\n"};

my $outfile = "$file" . "_$contig.afg";

open IN, $file;
open OUT, ">$outfile";

my $incontig = 0;
my $correctContig = 0;
my $count = 0;
my $readCount = 0;
my $foundReads = 0;
my $inRead = 0;
my $inCorrectRead = 0;
my $currentRead = 0;
my $len = -1;
my @reads;

print STDERR "Searching for start of contig $contig\n";

while(<IN>){
	$count ++;
	if(/{CTG/ || ($incontig &&/{RED/)){
		$incontig = &toggle($incontig);
		if(!$incontig){
			#contig has finished.. new contig has begun.
			$incontig = &toggle($incontig);
			if($correctContig){
				$correctContig = &toggle($correctContig);
				print STDERR "\nFound the end of contig: $contig.\n";
				print STDERR "Number of reads in contig $contig = $readCount\n";
			}
		}
	}
	if(/^iid:$contig$/ && $incontig && !$correctContig && !$foundReads){
		#this is the contig we want...
		#print everything to output file until we get the end of the contig..
		print STDERR "\nFound contig: $contig\n";
		$correctContig = &toggle($correctContig);
		#must print the header line...
		print OUT "{CTG\n";
	}
	if($correctContig){
		#need to collect all the read numbers for later on..
		if(/src:/){
			my @temp = split /:/, $_;
			chomp($temp[1]);
			push @reads, $temp[1];
			$readCount ++;
		}
		print OUT $_;
	}
	if(($count % 1000000 == 0)&& !$foundReads){
		print STDERR ".";
	}
	if(!$foundReads && /{RED/){
		print STDERR "\nFound the start of the reads..\n";
		print STDERR "Number of reads in contig $contig = $readCount\n";
		print STDERR "Now sorting reads.\n";
		my @sorted = sort { $a <=> $b } @reads;
		@reads = @sorted;
		print STDERR "Finished sorting reads.\n";
		$foundReads = &toggle($foundReads);
		$currentRead = shift(@reads);
		print STDERR "Current read: $currentRead\n";
	}
	if($inRead && /}/){
		#end of the read..
		$inRead = 0;
		if($inCorrectRead){
			print OUT "}\n";
			$inCorrectRead = 0;
			$currentRead = 0;
		}
		
	}
	if($foundReads && /{RED/){
		$inRead = 1;
		if($currentRead == 0){
			$currentRead = shift(@reads);
			$len = @reads;
			if($len % 100 == 0){
				print STDERR "Current read: $currentRead\t";
				print STDERR "length \@reads: $len\n"
			}
		}
	}
	if($inRead && /^iid:$currentRead$/){
		print OUT "{RED\n";
		$inCorrectRead = 1;
	}
	if($inCorrectRead){
		print OUT $_;
	}
	if($foundReads && ($len == 0)){
		last;
	}

	
}
