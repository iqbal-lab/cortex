#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use VCFFile;
use IntervalList;
use UsefulModule; # for num2str

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./vcf_filter_by_regions.pl [OPTIONS] <region1 ..> [in.vcf]
  Filter by regions.  Range is inclusive. " .
               "Does not assume VCF is sorted.  

  Options:
   --file         file containing one region per line
   --invert       invert selection
   --flag <flag>  Print all entries, flag hits

  Regions:
   chr:*          an entire chromosome
   chr:start-end  a region of a chromosome

  Examples:
  \$ ./vcf_regions.pl --flag IN_CHR1 chr1:10,000-12,000 data.vcf
  \$ ./vcf_regions.pl --file --invert regions.txt data.vcf
  \$ cat data.vcf | ./vcf_regions.pl chr1:* chr2:500-1000\n";

  exit;
}

## Test for filtering
my $failed_vars_out = undef;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  open($failed_vars_out, ">-");
}
##

if(@ARGV == 0)
{
  print_usage();
}

my $vcf_file;

my $invert = 0;
my $flag;

my %print_all_chrs = ();
my %regions_by_chr = ();

# Get options
while(@ARGV > 1 && $ARGV[0] =~ /^-/)
{
  my $arg = shift;

  if($arg =~ /^-?-invert$/i)
  {
    $invert = 1;
  }
  elsif(@ARGV == 1)
  {
    # All the oth
    print_usage("Invalid or missing arguments '$arg'");
  }
  elsif($arg =~ /^-?-file$/i)
  {
    my $file = shift;
    open(FILE, $file) or print_usage("Cannot open region file '$file'");
  
    my $line;
    while(defined($line = <FILE>))
    {
      chomp($line);
      # Check line is not empty or a comment line
      if($line !~ /^\s*$/ && $line !~ /^#/ && !parse_region($line))
      {
        print_usage("Unexpected region in file '$file': '$line'");
      }
    }

    close(FILE);
  }
  elsif($arg =~ /^-?-flag$/i)
  {
    $flag = shift;
  }
}

# Get regions

while(@ARGV > 0)
{
  my $arg = shift;

  if(!parse_region($arg))
  {
    if(@ARGV == 0)
    {
      # Last argument
      $vcf_file = $arg;
    }
    else
    {
      print_usage("Invalid commandline options\n");
    }
  }
}

# Create intervals
my %interval_lists_by_chr = ();

for my $chr (keys %regions_by_chr)
{
  $interval_lists_by_chr{$chr} = new IntervalList(@{$regions_by_chr{$chr}});
}

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

# Print non-PASS variants straight to stdout if -p passed
if(defined($failed_vars_out)) { $vcf->set_filter_failed($failed_vars_out);}

$vcf->print_header();

my $num_of_filtered_entries = 0;
my $total_num_entries = 0;

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  $total_num_entries++;

  my $chr = $vcf_entry->{'CHROM'};
  my $pos = $vcf_entry->{'true_POS'};
  my $len = length($vcf_entry->{'true_REF'});

  my $print = 0;

  if(defined($print_all_chrs{$chr}))
  {
    $print = 1;
  }
  elsif(defined($interval_lists_by_chr{$chr}))
  {
    my @hits = $interval_lists_by_chr{$chr}->fetch($pos, $pos+$len);
    $print = (@hits > 0);
  }

  if(defined($flag))
  {
    # Print all, flag those that 'hit'
    if($print != $invert)
    {
      $vcf_entry->{'INFO_flags'}->{$flag} = 1;
    }

    $vcf->print_entry($vcf_entry);
  }
  elsif($print != $invert)
  {
    $num_of_filtered_entries++;
    $vcf->print_entry($vcf_entry);
  }
}

# Print filtered rate
print STDERR "vcf_filter_by_regions.pl: " .
             pretty_fraction($num_of_filtered_entries, $total_num_entries) .
             " variants printed\n";

close($vcf_handle);



#
# Methods
#

sub parse_region
{
  my ($arg) = @_;

  if($arg =~ /^(.+):(?:([0-9,]+)-([0-9,]+)|\*)$/)
  {
    my $chr = $1;
    my $start = $2;
    my $end = $3;

    if(defined($start))
    {
      # Remove commas
      $start =~ s/,//g;
      $end =~ s/,//g;
    
      if($end < $start)
      {
        print_usage("end position is less than start position in: '$arg'\n");
      }

      if($start <= 0 || $end <= 0)
      {
        print_usage("'$arg' is not a valid region (1-based coords)");
      }

      if(!defined($regions_by_chr{$chr}))
      {
        $regions_by_chr{$chr} = [];
      }

      push(@{$regions_by_chr{$chr}}, [$start, $end+1]);
    }
    else
    {
      $print_all_chrs{$chr} = 1;
    }

    return 1;
  }
  else
  {
    return 0;
  }
}
