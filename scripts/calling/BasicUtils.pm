package BasicUtils;

use strict;
use warnings;

use base 'Exporter';
our @EXPORT = qw(check_cortex_runnable add_slash is_fasta is_fastq create_dir_if_does_not_exist count_bases_in_fasta  );


sub check_cortex_runnable
{
    my ($ctx_bin) = @_;
    
    my $cmd = "$ctx_bin --help";
    my $ret = qx{$cmd};
    if ($ret =~ /Starting Cortex/)
    {
	return; 
    }
    else
    {
	die("Cortex has been compiled, but the binary ($ctx_bin) cannot be run. This probably means some library (probably GSL) could not be found\n");
    }
}


sub add_slash
{
    my ($d)  = @_;
    if ($d ne "")
    {
	if ($d !~ /\/$/)
	{
	    $d= $d.'/';
	}
    }
    return $d;
}


sub is_fasta
{
    my ($file) = @_;
    
    my $cmd = "head -n 3 $file\n";
    my $ret = qx{$cmd};
    
    my @lines = split("\n", $ret);
    if (substr($lines[0], 0,1) ne ">")
    {
	return "EMissingGtrSign";
    }
    elsif ($lines[1] !~ /^[ACGTNacgtn]+$/)
    {
	return "ENotBases";
    }
    elsif ($lines[2] =~ /\+/)
    {
	return "EFastq";
    }
    else
    {
	return "EFasta";
    }
}


sub is_fastq
{
    my ($file) = @_;
    
    my $cmd = "head -n 5 $file\n";
    my $ret = qx{$cmd};
    
    my @lines = split("\n", $ret);
    if (substr($lines[0], 0,1) ne "@")
    {
	return "EMissingAtSign";
    }
    elsif ($lines[1] !~ /^[ACGTNacgtn]$/)
    {
	return "ENotBases";
    }
    elsif ($lines[2] !~ /[+]/)
    {
	return "ENotFastq";
    }
    elsif (length($lines[3]) != length($lines[0]))
    {
	return "EQualLineNotMatchSeqLineInFirstRead";
    }
    elsif (substr($lines[4], 0,1) ne "@")
    {
	return "ENotFastq";
    }
    else
    {
	return "EFastq";
    }
}

sub create_dir_if_does_not_exist
{
    my ($dir, $funcname) = @_;

    if (!(-d $dir))
    {
	my $cmd = "mkdir -p $dir";
	my $ret = qx{$cmd};
	
	if (!(-d $dir))
    	{
	    my $str = "Unable to create directory ";
	    $str .= "$dir in script/function $funcname\n";
	    $str .= "Error output:\n$ret\n";
	    die($str);
    	}	
    }
}

sub count_bases_in_fasta
{
	my ($file) = @_;
	my $tot=0;
	open(FILE, $file)||die("Function count_bases_in_fasta cannot open $file\n");
	while(<FILE>)
	{
		my $line = $_;
			chomp $line;
			if ($line !~ /^>/)
			{
				$tot += length($line);
			}
	}
	close(FILE);
	return $tot;
}

1;
