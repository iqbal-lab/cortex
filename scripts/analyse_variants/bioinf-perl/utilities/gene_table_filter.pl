#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(min max shuffle);

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"usage: ./gene_table_filter.pl <random|longest|first|last> [file]\n";

  exit;
}

if(@ARGV < 1 || @ARGV > 2)
{
  print_usage();
}

my $operation = shift;
my $file = shift;

if(!grep {/^(random|longest|first|last)$/i} $operation)
{
  print_usage("Invalid operation: $operation");
}

#
# Open Gene file
#
my $handle;

if(defined($file) && $file ne "-")
{
  open($handle, $file) or die("Cannot open VCF file '$file'\n");
}
elsif(-p STDIN) {
  # STDIN is connected to a pipe
  open($handle, "<&=STDIN") or die("Cannot read pipe");
}
else
{
  print_usage("Must specify or pipe in a gene file");
}

#
# Load all gene entries
#
my %genes_by_chr = ();

my $line;

my @columns = qw(bin key chr strand tx_start tx_end cds_start cds_end
                 exon_count exon_starts exon_ends score name2
                 cds_start_stat cds_end_stat exon_frames);

while(defined($line = <$handle>))
{
  chomp($line);
  my @gene_row = split(/\t/, $line);

  if(@gene_row != @columns)
  {
    print STDERR join(", ", @gene_row)."\n";
    die("Gene row has ".@gene_row." columns, expected ".@columns."\n");
  }

  my %gene = ();
  @gene{@columns} = @gene_row;

  my $chr = $gene{'chr'};

  if(!defined($genes_by_chr{$chr}))
  {
    $genes_by_chr{$chr} = [];
  }

  push(@{$genes_by_chr{$chr}}, \%gene);
}

close($handle);

#
# Sort
#

sub gene_cmp
{
  my $cmp = $a->{'tx_start'} - $b->{'tx_start'};

  if($cmp != 0)
  {
    return $cmp;
  }

  return ($a->{'tx_end'} - $b->{'tx_end'});
}

while(my ($chr,$list) = each(%genes_by_chr))
{
  my @sorted_genes = sort gene_cmp @$list;
  $genes_by_chr{$chr} = \@sorted_genes;
}

#
# Run and print!
#

for my $chr (sort keys %genes_by_chr)
{
  my $genes_list = $genes_by_chr{$chr};

  my $gene = $genes_list->[0];
  my $last_tx_end = $gene->{'tx_end'};

  my @overlapping_genes = ($gene);

  for(my $i = 1; $i < @$genes_list; $i++)
  {
    $gene = $genes_list->[$i];

    if($gene->{'tx_start'} <= $last_tx_end)
    {
      push(@overlapping_genes, $gene);
    }
    else
    {
      # Print a gene
      print_gene(\@overlapping_genes);

      @overlapping_genes = ($gene);
    }

    $last_tx_end = max($last_tx_end, $gene->{'tx_end'});
  }

  # Print a gene
  print_gene(\@overlapping_genes);
}


sub print_gene
{
  my ($genes_list) = @_;

  my $gene;

  if($operation =~ /first/i)
  {
    $gene = $genes_list->[0];
  }
  elsif($operation =~ /last/i)
  {
    $gene = $genes_list->[scalar(@$genes_list)-1];
  }
  elsif($operation =~ /random/i)
  {
    $gene = $genes_list->[int(rand(scalar(@$genes_list)))];
  }
  elsif($operation =~ /longest/i)
  {
    my @shuffled = shuffle(@$genes_list);
    my $idx_max = 0;
    $shuffled[$idx_max] > $shuffled[$_] or $idx_max = $_ for 1 .. $#shuffled;
    $gene = $genes_list->[$idx_max];
  }

  my %g = %$gene;

  print join("\t", @g{@columns}) . "\n";
}
