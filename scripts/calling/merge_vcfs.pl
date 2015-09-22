#!/usr/bin/perl -w
use strict;
use Getopt::Long;

##basically, how many args can we pass vcflib\s vcfcombine on the command line?
my $T = 10;

my $list = "";
my $vcflib_path="";
my $outfile = "";

&GetOptions(
    ##mandatory args
    'list_vcfs:s' =>\$list,
    'vcflib_path:s' => \$vcflib_path,
    'outfile:s' => \$outfile);


if ($list eq "")
{
    die("You must specify --list_vcfs - a filelist of VCFs to merge\n");
}

if ($vcflib_path eq "")
{
    die("You must specify --vcflib_path\n");
}

if ($outfile eq "")
{
    die("You must specify --outfile - the name of the VCF file you want to produce\n");
}


my $num_cmd = "wc -l $list";
my $number_vcfs = qx{$num_cmd};
chomp $number_vcfs;
if ($number_vcfs =~ /^(\d+)/)
{
    $number_vcfs = $1;
}
else
{
    die("You must use --list_vcfs\n");
}


merge_vcfs($list, $number_vcfs, $outfile, $T, $vcflib_path."/bin/vcfcombine",1);


sub get_all
{
    my ($f, $aref) = @_;
    open(F, $f)||die("Cannot open $f");
    while (<F>)
    {
	my $ln = $_;
	chomp $ln;
	push @$aref, $ln;
    }
    close(F);
}
sub remove_tmps
{
    my ($aref) = @_;
   foreach my $tmp (@$aref)
   {
       my $c = "rm $tmp";
       qx{$c};
   }
}
sub get_min
{
    my ($a, $b) = @_;
    if ($a<=$b)
    {
	return $a;
    }
    return $b;
}
sub merge_vcfs
{
    my ($list_all, $num_vcfs, $ofile, $thresh, $vcf_com, $level) = @_;

    if ($num_vcfs<$thresh)
    {
	my $cmd = "cat $list_all | xargs $vcf_com > $ofile";
	my $ret = qx{$cmd};
	print "Final step in merge \n$cmd\n$ret\n";
    }
    else
    {
	my @list=();
	get_all($list_all, \@list);##get list of VCFs into an array
	my $i=0;
	my @tmpfiles=();
	while ($i<$num_vcfs)
	{
	    my $j;
	    my $str=" ";
	    for ($j=$i; $j<get_min($num_vcfs,$i+$thresh); $j++)
	    {
		$str .= $list[$j]." ";
	    }
	    my $tmp = $ofile.".tmp.level$level".".".$i;
	    push @tmpfiles, $tmp;
	    my $cmd = "$vcf_com $str > $tmp" ;
	    qx{$cmd};
	    $i +=$thresh;
	}
	my $nextlist =$list_all.".level2"; 
	open(OUT, ">".$nextlist)||die("Cannot open $nextlist");
	foreach my $t (@tmpfiles)
	{
	    print OUT $t."\n";
	}
	close(OUT);
	merge_vcfs($nextlist, scalar(@tmpfiles), $ofile, $thresh, $vcf_com, $level+1);
	remove_tmps(\@tmpfiles);
	
    }

}
