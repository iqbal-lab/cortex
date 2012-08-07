#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use VCFFile;
use IntervalList;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./vcf_add_ensgene.pl <ensgene.txt> [file.vcf]\n" .
"  Label VCF entries if in a gene.  Adds GENE, GENES_STRAND & GENES_POS tags\n" .
"  to the INFO field.  If no [file.vcf] given (or given '-'), reads STDIN\n";
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

my $ensgene_file = shift;
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
# Open ensgene file
#

my %genes_by_chr = ();

my $ensgene_handle;
open($ensgene_handle, $ensgene_file)
  or print_usage("Cannot open ensgene file '$ensgene_file'");

my $ensgene_line;

while(defined($ensgene_line = <$ensgene_handle>))
{
  my (undef,$key,$chr,$strand,$tx_start,$tx_end,$cds_start,$cds_end,undef,
      $exon_starts_str,$exon_ends_str) = split(/\t/, $ensgene_line);

  my @exon_starts = split(",", $exon_starts_str);
  my @exon_ends = split(",", $exon_ends_str);

  my @exons = ();

  for(my $i = 0; $i < @exon_starts; $i++)
  {
    push(@exons, $exon_starts[$i], $exon_ends[$i]);
  }

  my $gene = {'key' => $key, 'chr' => $chr, 'strand' => $strand,
              'tx_start' => $tx_start, 'tx_end' => $tx_end,
              'cds_start' => $cds_start, 'cds_end' => $cds_end,
              'exons' => \@exons};

  if(!defined($genes_by_chr{$chr}))
  {
    $genes_by_chr{$chr} = [];
  }

  push(@{$genes_by_chr{$chr}}, [$tx_start,$tx_end+1,$gene]);
}

close($ensgene_handle);

my %gene_lists = ();

while(my ($chr, $genes) = each(%genes_by_chr))
{
  $gene_lists{$chr} = new IntervalList(@$genes);
}

#
# Read VCF
#
my $vcf = new VCFFile($vcf_handle);

# Print non-PASS variants straight to stdout if -p passed
if(defined($failed_vars_out)) { $vcf->set_filter_failed($failed_vars_out);}

$vcf->add_header_tag("INFO", "GENES", 1, "String", "Ensembl gene name/key");
$vcf->add_header_tag("INFO", "GENES_STRAND", 1, "String", "Strand gene is on +/-");
$vcf->add_header_tag("INFO", "GENES_POS", 1, "String", "Parts of gene in variant");
$vcf->print_header();

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  my $chr = $vcf_entry->{'CHROM'};
  my $start = $vcf_entry->{'true_POS'};
  my $end = $start + length($vcf_entry->{'true_REF'});

  my @genes = $gene_lists{$chr}->fetch($start, $end);

  if(@genes)
  {
    $vcf_entry->{'INFO'}->{'GENES'} = join(",", map {$_->{'key'}} @genes);
    $vcf_entry->{'INFO'}->{'GENES_STRAND'} = join(",", map {$_->{'strand'}} @genes);

    $vcf_entry->{'INFO'}->{'GENES_POS'}
      = join(",", map {get_gene_pos($start, $end, $_)} @genes);
  }

  $vcf->print_entry($vcf_entry);
}

close($vcf_handle);

sub get_gene_pos
{
  my ($start, $end, $gene) = @_;

  my @pos = ();

  if($end >= $gene->{'tx_start'} && $start < $gene->{'cds_start'})
  {
    push(@pos, $gene->{'strand'} eq "+" ? "TX_START" : "TX_END");
  }

  my $num_exons_pos = @{$gene->{'exons'}};

  if($num_exons_pos == 0)
  {
    if($end >= $gene->{'cds_start'} && $start <= $gene->{'cds_end'})
    {
      push(@pos, "CDS");
    }
  }
  else
  {
    if($start < $gene->{'exons'}->[0] && $end >= $gene->{'cds_start'})
    {
      push(@pos, "CDS_START");
    }
  
    my $num_exons = $num_exons_pos / 2;

    for(my $i = 0; $i < $num_exons_pos; $i += 2)
    {
      if($end >= $gene->{'exons'}->[$i] && $start <= $gene->{'exons'}->[$i+1])
      {
        my $exon_num = $gene->{'strand'} eq "+" ? ($i/2) : $num_exons - 1 - ($i/2);
        push(@pos, "EXON".$exon_num);
      }
      
      if($i+2 < $num_exons_pos && # Not the last exon
         $end >= $gene->{'exons'}->[$i+1] && $start <= $gene->{'exons'}->[$i+2])
      {
        my $intron_num = $gene->{'strand'} eq "+" ? ($i/2) : $num_exons - 2 - ($i/2);
        push(@pos, "INTRON".$intron_num);
      }
    }

    if($start < $gene->{'cds_end'} && $end >= $gene->{'exons'}->[$num_exons_pos-1])
    {
      push(@pos, "CDS_END");
    }
  }

  if($end > $gene->{'cds_end'} && $start <= $gene->{'tx_end'})
  {
    push(@pos, $gene->{'strand'} eq "+" ? "TX_END" : "TX_START");
  }

  return join(",", @pos);
}
