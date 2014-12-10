#!/usr/bin/perl -w
use strict;

my $list = shift; #list of VCFs

my %chrpos_to_stub=(); # chr_pos --> first 9 columns of VCF, which are same for all these VCFs
my %chrpos_to_sampleinfo=();
my %order = ();

### get stub - first  9 columns of the VCF
open(LIST, $list)||die();
my $first_file=<LIST>;
my ($header_stub, $num_lines) = get_stub($first_file, \%chrpos_to_stub, \%order);
close(LIST);


## get sample info
open(LIST, $list)||die();
my @samples=();
while (<LIST>)
{
    my $line = $_;
    chomp $line;
    my $sample = get_data($line, \%chrpos_to_sampleinfo);
    push @samples, $sample;
}
close(LIST);

## print combined VCF
print $header_stub."\t";
my $i;
for ($i=0; $i<scalar(@samples); $i++)
{
    print $samples[$i];
    if ($i<scalar(@samples)-1)
    {
	print "\t";
    }
    else
    {
	print "\n";
    }
}

my $j;
for ($j=1; $j<= $num_lines; $j++)
{
    my $chrpos = $order{$j};
    print $chrpos_to_stub{$chrpos};
    print "\t";
    #for each sample, print the genotype/etc info
    for ($i=0; $i<scalar(@samples); $i++)
    {
	print $chrpos_to_sampleinfo{$chrpos}{$samples[$i]};
	if ($i<scalar(@samples)-1)
	{
	    print "\t";
	}
	else
	{
	    print "\n";
	}
    }
}



sub get_data
{
    my ($file, $href) = @_;

    open(FILE, $file)||die();
    my $sample_id="";
    while (<FILE>)
    {
	my $line = $_;
	chomp $line;
	if ($line =~ /CHROM/)
	{
	    my @sp = split(/\s+/, $line);
	    $sample_id = $sp[9];
	}
	elsif ($line !~ /^\#/)
	{
	    my @sp = split(/\t/, $line);
	    my $sampleinfo = $sp[9];
	    my $chrpos = $sp[0]."_".$sp[1];
	    $href->{$chrpos}->{$sample_id} = $sampleinfo;
	}
	
    }
    close(FILE);
    return $sample_id;
}
sub get_stub
{
    my ($file, $href, $href_order) = @_;
    open(FILE, $file)||die();
    my $header = "";
    my $count=0;
    while (<FILE>)
    {
	my $lyne = $_;
	chomp $lyne;
	
	if ($lyne !~ /^\#/)
	{

	    my @sp = split(/\t/, $lyne);

	    $sp[7]=~ s/:/;/; ## fix a bug where the site confidence is preceded by colon

	    my $site_conf;
	    if ($sp[7]=~/SITE_CONF=(\S+)/)
	    {
		$site_conf=$1;
	    }
	    else
	    {
		die("Cannot parse $lyne properly - where is the site confidence?\n");
	    }
	    if ($site_conf<2)
	    {
		next;
	    }
	    ##if we get here, we are going to keep this call
	    $count++;
	    my @vals = @sp[0..8];
	    my $chrpos = $sp[0]."_".$sp[1];
	    $href->{$chrpos} = join("\t", @vals);
	    $href_order->{$count} = $chrpos;

	}
	else
	{
	    if ($lyne !~ /CHROM/)
	    {
		if ($lyne =~ /(models are variant, repeat and error)/)
		{
		    if ($lyne !~ /Calls are filtered so SITE_CONF is at least 2/)
		    {
			$lyne =~ s/\(models are variant, repeat and error\)/\(models are variant, repeat and error\). Calls are filtered so SITE_CONF is at least 2./;
		    }
		}
		$header = $header.$lyne."\n";
	    }
	    else
	    {
		my @sp = split(/\s+/, $lyne);
		my $j = join("\t", @sp[0..8]);
		$header = $header.$j;
	    }
	}
    }
    close(FILE);
    return ($header, $count);
}
