#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use VCFFile;

# config
my $csvsep = ",";
#

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_print_svlen.pl [--abs|--tag <t>] [file.vcf]\n";
  print STDERR "  Print size distribution of SVLENs\n";
  print STDERR "  --tag specify which INFO tag to use\n";
  exit;
}

## Test for filtering
my $skip_failed_vars = 0;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  $skip_failed_vars = 1;
}
##

if(@ARGV > 4)
{
  print_usage();
}

my $abs = 0;
my $tag = 'SVLEN';

while(@ARGV > 0)
{
  if($ARGV[0] =~ /^-?-abs$/i)
  {
    shift;
    $abs = 1;
  }
  elsif($ARGV[0] =~ /^-?-t(ag)?$/i)
  {
    shift;
    $tag = shift;
  }
  else
  {
    last;
  }
}

if(@ARGV > 1)
{
  print_usage();
}

my $vcf_file = shift;

#
# Open VCF Handle
#
my $vcf_handle;

if(defined($vcf_file) && $vcf_file ne "-")
{
  open($vcf_handle, $vcf_file)
    or print_usage("Cannot open VCF file '$vcf_file'\n");
}
elsif(-p STDIN) {
  # STDIN is connected to a pipe
  open($vcf_handle, "<&=STDIN") or print_usage("Cannot read pipe");
}
else
{
  print_usage("Must specify or pipe in a VCF file");
}

#
# Read VCF
#
my $vcf = new VCFFile($vcf_handle);

# Skip non-PASS variants if -p passed
if($skip_failed_vars) { $vcf->set_filter_failed(undef); }

my $vcf_entry;

my @ins = ();
my @del = ();

while(defined($vcf_entry = $vcf->read_entry()))
{
  my $svlen = $vcf_entry->{'INFO'}->{$tag};

  if($abs)
  {
    $svlen = abs($svlen);
  }

  if($svlen >= 0)
  {
    $ins[$svlen]++;
  }
  else
  {
    # deletion
    $svlen = abs($svlen);
    $del[$svlen]++;
  }
}

close($vcf_handle);

# Print header

print "size".$csvsep."count\n";

if(!$abs)
{
  print join("", map { "-".$_.$csvsep.(defined($del[$_]) ? $del[$_] : 0)."\n" }
                   reverse(1..$#del));
}

print join("", map {$_ . $csvsep . (defined($ins[$_]) ? $ins[$_] : 0) . "\n" }
                 1..$#ins);

