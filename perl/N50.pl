#!/usr/bin/perl -w
use strict;

my ($sum, $n, $result, @fields, $val, $th);

$th = 0;
$n = 0;
$sum = 0;
my ($sum50) = 0;
my (@n50) = ();


while (@fields = split(/\s/,<>)) {
	if (($val = $fields[$#fields]) =~ /(-?\d+)/) {
		$sum += $val;
		$n50[$n] = $val;
		++$n;
    }
}
my($vp) = .5;
if($n > 0) {
	$result = $sum / $n;
	print "sum = $sum, n = $n, ave = $result, ";
	$n = 0;
	foreach (sort { $b <=> $a } @n50) {
		$sum50 += $_;
		print "largest = $_\n" if $n == 0;
		$n++;
		if($sum50/$vp > $sum) {
			printf "N%.0f = $_, n = $n\n", $vp*100;
			$vp += .1;
		}
	}
}
