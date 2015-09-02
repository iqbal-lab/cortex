package ConfigMgmt;

use strict;
use warnings;

use base 'Exporter';
our @EXPORT = qw( get_from_config_if_undef print_to_config get_all_info_from_config_if_undef check_outdir);


## returns a pair (err, ret)
## err is EFILE if the file does not exist or cannot be opened, else EPASS
## ret is "" if the tag is not in the file, 
##      or is in the file but with no value (bad file)
##     else it returns the value.
sub get_from_config_if_undef
{
    my ($hashref, $str, $file) = @_;

    ## if the variable has a value, don't get a value from the config file
    if ($hashref->{$str} ne "")
    {
	return "EPASS";
    }
    if (!(-e $file))
    {
	return "EFILE";
    }
    open(FILE, $file)||die("Unable to open config file $file\n");
    while (<FILE>)
    {
	my $line = $_;
	chomp $line; 
	my @sp = split(/\t/, $line);
	my $tag = $sp[0];
	my $val = $sp[1];
	if ($tag eq $str)
	{
	    if (scalar @sp !=2)
	    {
		return "EFILE";
	    }
	    $hashref->{$str}=$val;
	    return "EPASS";
	}
    }
    close(FILE);
    return "EPASS";
}

sub print_to_config
{
    my ($str, $hashref, $fh) = @_;
    print $fh $str."\t".$hashref->{$str};
    print $fh "\n";
}

sub get_all_info_from_config_if_undef
{
    my ($hashref, $file) = @_;

    if (!( -e $file))
    {
	return;## just return if the file does not exist
    }
    open(FILE, $file) ||die("Cannot open file $file in get_all_from_config_if_undef\n");
    while (<FILE>)
    {
	my $line = $_;
	chomp $line;
	my @sp = split(/\t/, $line);
        my $tag = $sp[0];
	my $val = $sp[1];
	if (!exists $hashref->{$tag})
	{
	    $hashref->{$tag}=$val;
	}
	else
	{
	    if ($hashref->{$tag} eq "")
	    {
		$hashref->{$tag}=$val;
	    }
	}
    }
    close(FILE);
}


sub check_outdir
{
    my ($hashref) = @_;
    if ($hashref->{"outdir"} eq "")
    {
	my $errstr = "You must specify with --outdir the output directory, which \n";
	$errstr .= "must be the same directory used as output dir in prepare.pl\n";
	die($errstr);
    }
    $hashref->{"outdir"} = BasicUtils::add_slash($hashref->{"outdir"});
    if (!(-d $hashref->{"outdir"} ))
    {
	my $errstr = "Output directory does not exist:".$hashref->{"outdir"};
	$errstr .= "\n This should have been created previously by prepare.pl\n";
	die($errstr);
    }
}

1;
