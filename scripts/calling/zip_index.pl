#!/usr/bin/perl -w
use strict;

my $file = shift;

my $save = "cp $file $file".".saved";
qx{$save};

my $bgz = "bgzip -c $file > $file".".gz";
my $bgz_ret = qx{$bgz};
#print "$bgz\n$bgz_ret\n";
my $tab = "tabix -p vcf $file".".gz";
my $tab_ret = qx{$tab};
#print "$tab\n$tab_ret\n";
