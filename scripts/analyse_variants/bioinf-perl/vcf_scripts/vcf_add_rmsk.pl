#!/usr/bin/perl

use strict;
use warnings;

use List::MoreUtils qw(uniq);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use VCFFile;
use GeneticsModule;
use IntervalList;
use UsefulModule; # num2str

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "Usage: ./vcf_add_rmsk.pl <rmsk.txt> [in.vcf]\n";
  print STDERR "  Add repeat annoations to variants, using INFO tags:\n";
  print STDERR "  - rmsk_CLASS_left, rmsk_CLASS_right\n";
  print STDERR "    (distances (bp) to the left/right to nearest element\n";
  exit;
}

## Test for filtering
my $failed_vars_out = undef;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  open($failed_vars_out, ">-");
}
##

if(@ARGV < 1 || @ARGV > 2)
{
  print_usage();
}

my $rmsk_file = shift;
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
# Load rmsk.txt file
#

open(RMSK, $rmsk_file) or die("Cannot read rmsk file '$rmsk_file'");

my %num_per_class = ();
my %repeat_elements_by_chr = ();

my $rmsk_line;

while($rmsk_line = <RMSK>)
{
  my @rmsk_cols = split(/\t/, $rmsk_line);

  my $rmsk_chr = get_clean_chr_name($rmsk_cols[5]);
  my $rmsk_start = $rmsk_cols[6];
  my $rmsk_end = $rmsk_cols[7];
  #my $rmsk_class = $rmsk_cols[11];
  my $rmsk_class = $rmsk_cols[12];

  if($rmsk_class =~ /\?$/)
  {
    next;
  }

  $num_per_class{$rmsk_class}++;

  my $repeat = {'chr' => $rmsk_chr,
                'start' => $rmsk_start,
                'end' => $rmsk_end,
                'class' => $rmsk_class};

  if(!defined($repeat_elements_by_chr{$rmsk_chr}))
  {
    $repeat_elements_by_chr{$rmsk_chr} = [];
  }

  push(@{$repeat_elements_by_chr{$rmsk_chr}},
       [$rmsk_start, $rmsk_end+1, $repeat]);
}

close(RMSK);

my @rmsk_classes = sort {$a cmp $b} keys %num_per_class;

my %intervals_per_chrom = ();

for my $chr (keys %repeat_elements_by_chr)
{
  $intervals_per_chrom{$chr}
    = new IntervalList(@{$repeat_elements_by_chr{$chr}});
}

#
# Read VCF
#
my $vcf = new VCFFile($vcf_handle);

# Print non-PASS variants straight to stdout if -p passed
if(defined($failed_vars_out)) { $vcf->set_filter_failed($failed_vars_out);}

# Add header tags
$vcf->add_header_tag("INFO", "rmsk", 1, "String", "RMSK elements a variant is inside");

# Print VCF header
$vcf->print_header();

my $total_num_entries = 0;
my %missing_chrs = ();

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  $total_num_entries++;

  my $chr = $vcf_entry->{'CHROM'};

  my $var_start = $vcf_entry->{'true_POS'};
  my $var_end = $var_start + length($vcf_entry->{'true_REF'});

  my $interval_list = $intervals_per_chrom{$chr};

  if(!defined($interval_list))
  {
    if(!defined($missing_chrs{$chr}))
    {
      print STDERR "Chromosome '$chr' is missing from rmsk file '$rmsk_file'\n";
      $missing_chrs{$chr} = 1;
    }
    next;
  }

  my @hits = $interval_list->fetch($var_start, $var_end);
  
  if(@hits > 0)
  {
    $vcf_entry->{'INFO'}->{'rmsk'} = join(",",map {$_->{'class'}} @hits);
  }

  $vcf->print_entry($vcf_entry);
}

close($vcf_handle);
