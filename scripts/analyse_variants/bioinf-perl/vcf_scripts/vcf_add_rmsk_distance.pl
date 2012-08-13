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

  print STDERR "" .
"Usage: ./vcf_add_rmsk_distance.pl <rmsk.txt> [in.vcf]
  Add repeat annoations to variants, using INFO tags:
  - rmsk_CLASS_left, rmsk_CLASS_right
    (distances (bp) to the left/right to nearest element\n";

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
  my $rmsk_class = $rmsk_cols[11];

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
    $repeat_elements_by_chr{$rmsk_chr} = {};
  }
  
  if(!defined($repeat_elements_by_chr{$rmsk_chr}->{$rmsk_class}))
  {
    $repeat_elements_by_chr{$rmsk_chr}->{$rmsk_class} = [];
  }

  push(@{$repeat_elements_by_chr{$rmsk_chr}->{$rmsk_class}},
       [$rmsk_start, $rmsk_end+1, $repeat]);
}

close(RMSK);

my @rmsk_classes = sort {$a cmp $b} keys %num_per_class;

my %interval_lists = ();

for my $rmsk_chr (keys %repeat_elements_by_chr)
{
  for my $rmsk_class (@rmsk_classes)
  {
    $interval_lists{$rmsk_chr}->{$rmsk_class}
      = new IntervalList(@{$repeat_elements_by_chr{$rmsk_chr}->{$rmsk_class}});
  }
}

#
# Read VCF
#
my $vcf = new VCFFile($vcf_handle);

# Print non-PASS variants straight to stdout if -p passed
if(defined($failed_vars_out)) { $vcf->set_filter_failed($failed_vars_out);}

# Add header tags
for my $class (@rmsk_classes)
{
  my $descrip = "Distance (bp) left to nearset element '$class' (0 => inside)";
  $vcf->add_header_tag("INFO", "rmsk_".$class."_left", 1, "Integer", $descrip);
  $descrip = "Distance (bp) right to nearset element '$class' (0 => inside)";
  $vcf->add_header_tag("INFO", "rmsk_".$class."_right", 1, "Integer", $descrip);
}

# Print VCF header
$vcf->print_header();

my $total_num_entries = 0;
my %num_in_repeat_class = ();
my %missing_chrs = ();

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  $total_num_entries++;

  my $chr = $vcf_entry->{'CHROM'};

  my $var_start = $vcf_entry->{'true_POS'};
  my $var_end = $var_start + length($vcf_entry->{'true_REF'});

  my $interval_lists_hashref = $interval_lists{$chr};

  if(!defined($interval_lists_hashref))
  {
    if(!defined($missing_chrs{$chr}))
    {
      print STDERR "Chromosome '$chr' is missing from rmsk file '$rmsk_file'\n";
      $missing_chrs{$chr} = 1;
    }
    next;
  }

  for my $rmsk_class (@rmsk_classes)
  {
    my ($hits_arr, $left_arr, $right_arr)
      = $interval_lists_hashref->{$rmsk_class}->fetch_nearest($var_start,
                                                              $var_end);

    
    my $left_dist = 0;
    my $right_dist = 0;
    
    if(@$hits_arr == 0)
    {
      # Get distance to hits either side, or '.' if there were none
      $left_dist = @$left_arr > 0 ? $var_start - $left_arr->[0]->{'end'} : ".";
      $right_dist = @$right_arr > 0 ? $right_arr->[0]->{'start'} - $var_end : ".";
    }
    else
    {
      $num_in_repeat_class{$rmsk_class}++;
    }

    $vcf_entry->{'INFO'}->{'rmsk_'.$rmsk_class.'_left'} = $left_dist;
    $vcf_entry->{'INFO'}->{'rmsk_'.$rmsk_class.'_right'} = $right_dist;
  }

  $vcf->print_entry($vcf_entry);
}


print STDERR "vcf_add_repeat_masker.pl: of " .
             num2str($total_num_entries) . " VCF entries:\n";

my @sorted_classes = sort {$num_in_repeat_class{$b} <=> $num_in_repeat_class{$a}}
                     keys %num_in_repeat_class;

@sorted_classes = uniq(@sorted_classes, @rmsk_classes);

for my $class (@sorted_classes)
{
  my $count = $num_in_repeat_class{$class};

  if(!defined($count))
  {
    $count = 0;
  }
  
  my $percent = 100 * $count / $total_num_entries;

  print STDERR "  " . num2str($count) .
               " (" . sprintf("%.2f", $percent) . "%) " .
               " are in repeat class $class\n";
}

close($vcf_handle);
