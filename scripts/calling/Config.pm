package Config;

use strict;
use warnings;

use base 'Exporter';
our @EXPORT = qw( get_from_config_if_undef );


## returns a pair (err, ret)
## err is EFILE if the file does not exist or cannot be opened, else EPASS
## ret is "" if the tag is not in the file, 
##      or is in the file but with no value (bad file)
##     else it returns the value.
sub get_from_config_if_undef
{
    my ($var, $str, $file) = @_;

    ## if the variable has a value, don't get a value from the config file
    if ($var ne "")
    {
	return (EPASS, $var);
    }
    if (!(-e $file))
    {
	return (EFILE, "");
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
		return (EFILE, "");
	    }
	    return (EPASS, $val);
	}
    }
    close(FILE);
    return (EPASS, "");
}

1;
