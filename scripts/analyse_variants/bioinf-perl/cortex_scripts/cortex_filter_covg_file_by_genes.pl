#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(min);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use CortexCovgFile;
use UsefulModule;

## Config
my $csvsep = ",";
my $padding = 100; # bp away from a gene
#

sub print_usage
{
  for my $err (@_)
  {
    print STDERR "Error: $err\n";
  }
  
  print STDERR "Usage: ./cortex_filter_covg_file_by_genes.pl <ensGene.txt> " .
               "<5p_mapping.sam> <outgenes.out> <ingenes.out> <unmapped.out> [.colour_covgs]\n";
  print STDERR "  Print colour covg entries that map $padding\bp away from genes\n";
  exit;
}

if(@ARGV < 5 || @ARGV > 6)
{
  print_usage();
}

my $gene_file = shift;
my $mapping_file = shift;
my $outgenes_out = shift;
my $ingenes_out = shift;
my $unmapped_out = shift;
my $covg_file = shift;

#
# Open .colour_covgs handle
#
my $covg_handle;

if(defined($covg_file) && $covg_file ne "-")
{
  open($covg_handle, $covg_file)
    or print_usage("Cannot open .colour_covgs file '$covg_file'\n");
}
elsif(-p STDIN) {
  # STDIN is connected to a pipe
  open($covg_handle, "<&=STDIN") or print_usage("Cannot read pipe");
}
else
{
  print_usage("Must specify or pipe in a Cortex .colour_covgs file");
}

my $covgfile = new CortexCovgFile($covg_handle);

#
# Open mapping file
#
my $mapping_handle;
open($mapping_handle, $mapping_file)
  or die("Cannot open mapping file '$mapping_file'");

#
# Load txstart, txend coords from ensGene.txt file
#
#  `bin` smallint(5) unsigned NOT NULL,
#  `name` varchar(255) NOT NULL,
#  `chrom` varchar(255) NOT NULL,
#  `strand` char(1) NOT NULL,
#  `txStart` int(10) unsigned NOT NULL,
#  `txEnd` int(10) unsigned NOT NULL,
#  ...

my %gene_starts_by_chr = ();
my %gene_ends_by_chr = ();

open(GENES, $gene_file) or die("Couldn't open gene file '$gene_file'");

my $gene_line;
while(defined($gene_line = <GENES>))
{
  my ($b,$name,$chrom,$strand,$txstart,$txend) = split(/\t/, $gene_line);

  if(!defined($gene_starts_by_chr{$chrom}))
  {
    $gene_starts_by_chr{$chrom} = [];
    $gene_ends_by_chr{$chrom} = [];
  }

  push(@{$gene_starts_by_chr{$chrom}}, $txstart);
  push(@{$gene_ends_by_chr{$chrom}}, $txend);
}

close(GENES);

# sort by gene start
for my $chrom (keys %gene_starts_by_chr)
{
  my @starts = @{$gene_starts_by_chr{$chrom}};
  my @ends = @{$gene_ends_by_chr{$chrom}};

  my @order = sort {$starts[$a] <=> $starts[$b]} 0..$#starts;

  @starts = @starts[@order];
  @ends = @ends[@order];

  $gene_starts_by_chr{$chrom} = \@starts;
  $gene_ends_by_chr{$chrom} = \@ends;
}

#
# Open output files
#
open(OUTGENES, ">$outgenes_out")
  or print_usage("Cannot open filtered output file '$outgenes_out'");

open(INGENES, ">$ingenes_out")
  or print_usage("Cannot open genes output file '$ingenes_out'");

open(UNMAPPED, ">$unmapped_out")
  or print_usage("Cannot open genes output file '$unmapped_out'");

#
# Start reading bubbles
#
my $bubble_entry = $covgfile->grab_bubble_entry();
my $mapping_entry;

my $num_bubbles_outside_genes = 0;
my $num_bubbles_inside_genes = 0;
my $num_bubbles_unmapped = 0;
my $num_of_bubbles = 0;

while(defined($mapping_entry = <$mapping_handle>) && $mapping_entry =~ /^@/) {}

while(defined($bubble_entry))
{
  if(!defined($mapping_entry))
  {
    print STDERR "bubble: '$bubble_entry'\n";
    print_usage("Ran out of mapping entries (on bubble $num_of_bubbles)!");
  }

  my ($name,$flags,$ref_name,$pos,$mapQ,,,,,$seq) = split(/\t/, $mapping_entry);

  if($bubble_entry =~ /^>(var_\w+)/i)
  {
    if($1 ne $name)
    {
      die("Bubble names don't match ('$1' in colour_covgs, '$name' in mapping)");
    }
  }
  else
  {
    print STDERR "bubble: '$bubble_entry'\n";
    print_usage("Unexpected bubble");
  }

  my $print_bubble = 1;

  if(!($flags & 0x4))
  {
    # Flank is mapped
    
    if(!($flags & 0x10))
    {
      # If fragment hasn't been inverted, then this pos of variant is
      $pos += length($seq);
    }
    
    my $gene_starts = $gene_starts_by_chr{$ref_name};
    my $gene_ends = $gene_starts_by_chr{$ref_name};

    # check if $ref_name:$pos is within $padding bases of a gene
    my $nearest_gene = binary_search_nearest($gene_starts, $pos);
    my $dist;

    if($gene_starts->[$nearest_gene] <= $pos)
    {
      if($gene_ends->[$nearest_gene] > $pos)
      {
        # In gene
        $dist = -1;
      }
      else
      {
        # Check dist from end of gene
        $dist = $pos - $gene_ends->[$nearest_gene];
      
        if($nearest_gene+1 < @$gene_starts)
        {
          # Check dist before next gene start
          $dist = min($dist, $gene_starts->[$nearest_gene+1] - $pos);
        }
      }
    }
    else
    {
      # Check dist before start of gene
      $dist = $gene_starts->[$nearest_gene] - $pos;

      if($nearest_gene > 0)
      {
        # Check previous gene
        if($gene_ends->[$nearest_gene-1] > $pos)
        {
          # In previous gene
          $dist = -1;
        }
        else
        {
          # check distance after end of previous gene
          $dist = min($dist, $pos - $gene_ends->[$nearest_gene-1]);
        }
      }
    }

    if($dist > $padding)
    {
      # outside of $padding bases of a transcribed region
      print OUTGENES $bubble_entry;
      $num_bubbles_outside_genes++;
    }
    else
    {
      # within $padding bases of a transcribed region
      print INGENES $bubble_entry;
      $num_bubbles_inside_genes++;
    }
  }
  else
  {
    # unmapped
    print UNMAPPED $bubble_entry;
      $num_bubbles_unmapped++;
  }
  
  $bubble_entry = $covgfile->grab_bubble_entry();
  $mapping_entry = <$mapping_handle>;
  $num_of_bubbles++;
}

close(OUTGENES);
close(INGENES);
close(UNMAPPED);

close($covg_handle);
close($mapping_handle);

print "Outside genes: $num_bubbles_outside_genes / $num_of_bubbles bubbles " .
      "(" . ($num_bubbles_outside_genes / $num_of_bubbles) . ")\n";

print "Inside genes:  $num_bubbles_inside_genes / $num_of_bubbles bubbles " .
      "(" . ($num_bubbles_inside_genes / $num_of_bubbles) . ")\n";

print "unmapped:      $num_bubbles_unmapped / $num_of_bubbles bubbles " .
      "(" . ($num_bubbles_unmapped / $num_of_bubbles) . ")\n";

if(defined($bubble_entry))
{
  print STDERR "Warning: still more bubbles entries\n";
  print STDERR "entry: '$bubble_entry'\n";
}

if(defined($mapping_entry))
{
  print STDERR "Warning: still more mapping entries\n";
  print STDERR "entry: '$mapping_entry'\n";
}
