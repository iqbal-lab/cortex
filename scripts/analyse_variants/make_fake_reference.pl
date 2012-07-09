#!/usr/bin/perl -w
use strict;
use Getopt::Long;


my $callfile="";
my $outfile="";
my $help='';

&GetOptions(
	'callfile|f:s'                         => \$callfile,
	'outfile|o:s'                           => \$outfile,
	'help'                                 => \$help,
);

if ($help)
{
    print "Script for making a \"fake\" reference, with one chromosome PER VARIANT CALL.\nThis is purely a technical trick to allow us to use the VCF file format even when we do not have a reference genome/coordinates\n";
    print "Usage: perl make_fake_reference.pl --callfile <filename of output of bubble caller or PD caller, does no matter if you used --print_colour_coverages>\n";
    print "                                   --outfile <filename of the file you want to generate>\n";
    print "\n\n Please note that this file expects to read a callfile with numerically named variants, containing \"var_\" and then a number which is unique to the callset\n such as var_1, var_2, or SOMENAME_var_1, SOMENAME_var_2, etc.\nThis is what the BC and PD do,\nand also what the run_calls pipeline generates too. So you only need to worry if you are trying something clever ;-) \n";
    exit(0);
}


check_args($callfile, $outfile);

print_reference($callfile, $outfile);

print "There you go. One fake reference file $outfile completed! Have fun :-) \n";


sub print_reference
{
    my ($callf, $outf) = @_;
    
    open(IN, $callf)||die("Unable to open the callfile $callf. But I have checked it exists, so maybe a permissions or network issue? Aborting. Please fix.\n");
    open(OUT, ">".$outf)||die("Unable to open the output file $outf. Maybe a permissions issue? Aborting. Please fix.\n");

    while (<IN>)
    {
	my $line = $_;
	chomp $line;
	if ($line =~ /var_(\d+)_5p_flank/)
	{
	    my $num = $1;
	    my $seq="";
	    $line = <IN>;
	    chomp $line;
	    $seq = $seq.$line;# now we have the 5p flank
	    <IN>; #ignore the branch one read id
	    $line = <IN>;
	    chomp $line;
	    $seq = $seq.$line;# now we have 5p flank + branch1
	    <IN>;
	    <IN>; #ignore branch2
	    $line = <IN>;
	    if ($line !~ /3p_flank/)
	    {
		die("SOme kind of formatting issue with the callfile, at line $line for variant number $num. Expected to see the 3p flank read id\n");
	    }
	    $line = <IN>;
	    chomp $line;
	    $seq = $seq.$line;# now we have 5p flank + branch1 + 3p flank
	    print OUT ">chr$num\n";
	    print OUT "$seq\n";
	    
	}

    }
    close(IN);
    close(OUT);

    if (!(-e($outf)))
	{
	    die("Somehow have failed to write the output file. The only thing I can imagine causing this error is the disk being full or a network issue\n");
	}
}

sub check_args
{
    my ($c, $f) = @_;
    
    if ($c eq "")
    {
	die("You  must specify --callfile and give a file as argument\n");
    }

    if ($f eq "")
    {
	die("You  must specify --outfile and give a file as argument\n");
    }

    if (!(-e $c))
    {
	die("The callfile you specify, $c, does not seem to exist\n");
    }
    if (-e($f))
    {
	die("The output filename you specify, $f, already exists, and I do not want to overwrite it. Delete it and then rerun this script\n");
    }

    return;
}
