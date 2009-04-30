#!/usr/local/bin/perl5.6.0 -w
use strict;

my ($i, $rp, $count, $name, $oname
	);
 
$name = "";
$oname = "";
my (%seq) = ();
my ($mindiff) = 0;
foreach $name (@ARGV) {
	open F,$name;
	my ($tag) = "";
	my ($tail) = "";
	my ($length) = 0;
	my ($pos) = -1;
	while (<F>) {
		chomp;
		my ($line) = $_;
		if ($line =~ /^\>(\S+)\s*(.*)/) {
			if($pos >= 0) {
				print "$tag $pos $length $name $tail\n";
				$length = 0;
			}
			$tag = $1;
			$tail = "";
			$tail = $2 if defined $2;
			$pos = tell(F) - length($line) - 1;
			next;
		}
		else{
#		if ($line =~ /^[acgtACGTnN]/) {
			$length += length($line);
			next;
		}
	}
	close F;
	print "$tag $pos $length $name $tail\n";
}
