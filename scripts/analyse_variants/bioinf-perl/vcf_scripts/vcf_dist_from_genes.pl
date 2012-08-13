#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use VCFFile;
use UsefulModule;

## Config
my $csvsep = "\t";
##

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }
  
  print STDERR "" .
"Usage: ./vcf_dist_from_genes.pl <ensGene.txt> <CDS/TX> " .
  "<chr_lengths.csv> <binsize> <out> [in.vcf] 
  Get distribution of variants around genes
  <CDS/TX> whether to use coding or transcription start
  <binsize> is the size of bins in kbp
  Saves: <out>.upstream.csv, <out>.downstream.csv,
         <out>.intron.csv, <out>.exon.csv,
         <out>.preexon.csv & <out>.postexon.csv>\n";
  exit;
}

## Test for filtering
my $skip_failed_vars = 0;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  $skip_failed_vars = 1;
}
##

if(@ARGV < 5 || @ARGV > 6)
{
  print_usage();
}

my $genes_file = shift;
my $cds_or_tx = shift;
my $chr_lengths_file = shift;
my $bin_size_kbp = shift;
my $out_base = shift;
my $vcf_file = shift;

my $use_cds;

if($cds_or_tx =~ /^CDS$/i)
{
  $use_cds = 1;
}
elsif($cds_or_tx =~ /^TX$/i)
{
  $use_cds = 0;
}
else
{
  print_usage("Invalid cds/tx option: '$cds_or_tx' - should be 'CDS' or 'TX'");
}


if($bin_size_kbp !~ /^\d+$/)
{
  print_usage("Invalid <binsize> value ('$bin_size_kbp')");
}
elsif($bin_size_kbp >= 10000)
{
  warn("<binsize> units are kbp - $bin_size_kbp kbp is probably too large");
}

#
# Load chr lengths
#
my %chr_lengths = ();

open(CHRS, $chr_lengths_file)
  or print_usage("Cannot open chr_lengths CSV file '$chr_lengths_file'");

my $line;

if(!defined($line = <CHRS>))
{
  print_usage("Empty rmsk.txt file! (file: $chr_lengths_file)");
}

if($line !~ /^\w+,\s*\d+$/)
{
  # Skip header?
  $line = <CHRS>;
}

while(defined($line = <CHRS>))
{
  if($line =~ /^(\w+),\s*(\d+)$/)
  {
    $chr_lengths{$1} = $2;
  }
  else {
    print_usage("Unexpected line in chr_lengths CSV file: '$line'");
  }
}

close(CHRS);

my $bin_size = $bin_size_kbp * 1000;
# TESTING
#$bin_size = 10;

print "Bin size ".num2str($bin_size)." bp\n";

#
# Open VCF Handle
#
my $vcf_handle;

if(defined($vcf_file) && $vcf_file ne "-") {
  open($vcf_handle, $vcf_file) or die("Cannot open VCF file '$vcf_file'\n");
}
elsif(-p STDIN) {
  # STDIN is connected to a pipe
  open($vcf_handle, "<&=STDIN") or die("Cannot read pipe");
}
else
{
  print_usage("Must specify or pipe in a VCF file");
}

#
# Open output files
#
my $upstream_file = $out_base.'.upstream.csv';
my $downstream_file = $out_base.'.downstream.csv';
my $ingenes_file = $out_base.'.ingenes.csv';
my $exon_file = $out_base.'.exon.csv';
my $intron_file = $out_base.'.intron.csv';
my $preexon_file = $out_base.'.preexon.csv';
my $postexon_file = $out_base.'.postexon.csv';

open(UPSTREAM, ">$upstream_file")
  or print_usage("Cannot open output file '$upstream_file'");

open(DOWNSTREAM, ">$downstream_file")
  or print_usage("Cannot open output file '$downstream_file'");

open(INGENES, ">$ingenes_file")
  or print_usage("Cannot open output file '$ingenes_file'");

open(EXONS, ">$exon_file")
  or print_usage("Cannot open output file '$exon_file'");

open(INTRONS, ">$intron_file")
  or print_usage("Cannot open output file '$intron_file'");

open(PREEXON, ">$preexon_file")
  or print_usage("Cannot open output file '$preexon_file'");

open(POSTEXON, ">$postexon_file")
  or print_usage("Cannot open output file '$postexon_file'");

#
# Load gene / exon file
#
my $genes_handle;

open($genes_handle, $genes_file)
  or print_usage("Cannot open genes file '$genes_file'");

print "Reading genes file...\n";

my $genes_line;

my %chr_genes = ();

while(defined($genes_line = <$genes_handle>))
{
  chomp($genes_line);

  my ($bin, $name, $chrom, $strand, $txStart, $txEnd, $cdsStart, $cdsEnd,
      $exonCount, $exonStarts, $exonEnds, $score, $name2,
      $cdsStartStat, $cdsEndStat, $exonFrames) = split(/\t/, $genes_line);

  if(!defined($exonFrames))
  {
    die("Not enough lines in gene file '$genes_file' - " .
        "sure it's ensGene.txt from UCSC?");
  }

  if(!defined($chr_genes{$chrom}))
  {
    $chr_genes{$chrom} = [];
  }

  my @exon_starts = split(/,/, $exonStarts);
  my @exon_ends = split(/,/, $exonEnds);

  my %gene_entry = ('start' => $use_cds ? $cdsStart : $txStart,
                    'end' => $use_cds ? $cdsEnd : $txEnd,
                    'strand' => $strand,
                    'exon_starts' => \@exon_starts,
                    'exon_ends' => \@exon_ends);

  push(@{$chr_genes{$chrom}}, \%gene_entry);
}

close($genes_handle);

#
# Read VCF
#
print "Reading VCF...\n";

my $vcf = new VCFFile($vcf_handle);

# Skip non-PASS variants if -p passed
if($skip_failed_vars) { $vcf->set_filter_failed(undef); }

my @positions = ();

my @upstream_bins = ();
my @upstream_denominator = ();
my @downstream_bins = ();
my @downstream_denominator = ();

# Positions in genes
my @pos_in_genes = ();

my @pre_exons = ();
my @exon_distances = ();
my @intron_distances = ();
my @post_exons = ();

my $vcf_entry = $vcf->read_entry();
my $curr_chrom = $vcf_entry->{'CHROM'};

do
{
  if($vcf_entry->{'CHROM'} ne $curr_chrom)
  {
    # Analyse positions on given chromosome
    if(defined($chr_genes{$curr_chrom}))
    {
      add_to_gene_bins($curr_chrom);
    }

    $curr_chrom = $vcf_entry->{'CHROM'};
    @positions = ();
  }

  push(@positions, $vcf_entry->{'true_POS'});
}
while(defined($vcf_entry = $vcf->read_entry()));

if(defined($chr_genes{$curr_chrom}))
{
  add_to_gene_bins($curr_chrom);
}

close($vcf_handle);

#
# SAVE!
#

# Save upstream histogram
print "Saving upstream results...\n";

# print header
print UPSTREAM "bin_index,bin_start,variant_count,bin_count\n";

for(my $i = 0; $i < @upstream_denominator; $i++)
{
  my $count = defined($upstream_bins[$i]) ? $upstream_bins[$i] : 0;
  my $denom = $upstream_denominator[$i];

  print UPSTREAM "$i,".($i*$bin_size).",$count,$denom\n";
}

close(UPSTREAM);

# Save downstream histogram
print "Saving downstream results...\n";

# print header
print DOWNSTREAM "bin_index,bin_start,variant_count,bin_count\n";

for(my $i = 0; $i < @downstream_denominator; $i++)
{
  my $count = defined($downstream_bins[$i]) ? $downstream_bins[$i] : 0;
  my $denom = $downstream_denominator[$i];

  print DOWNSTREAM "$i,".($i*$bin_size).",$count,".sprintf("%.2g",$denom)."\n";
}

close(DOWNSTREAM);

# Save in genes
print "Saving in gene results...\n";

@pos_in_genes = sort {$a <=> $b} @pos_in_genes;

for my $dist (@pos_in_genes) {
  print INGENES "$dist\n";
}

close(INGENES);

# Save exons
print "Saving exon results...\n";

@exon_distances = sort {$a <=> $b} @exon_distances;

for my $dist (@exon_distances) {
  print EXONS "$dist\n";
}

close(EXONS);

# Save introns
print "Saving intron results...\n";

@intron_distances = sort {$a <=> $b} @intron_distances;

for my $dist (@intron_distances) {
  print INTRONS "$dist\n";
}

close(INTRONS);

# Save pre-exon
print "Saving pre-exon results...\n";

@pre_exons = sort {$a <=> $b} @pre_exons;

for my $dist (@pre_exons) {
  print PREEXON "$dist\n";
}

close(PREEXON);

# Save post-exon
print "Saving post-exon results...\n";

@post_exons = sort {$a <=> $b} @post_exons;

for my $dist (@post_exons) {
  print POSTEXON "$dist\n";
}

close(POSTEXON);

#
# That's it folks!
#

print "Done.\n";



sub add_to_gene_bins
{
  my ($chr) = @_;

  print "Processing $chr\n";

  my $chr_length = $chr_lengths{$chr};

  # TESTING
  #$chr_length = 200;

  #print "$chr: ".join(",",@positions)."\n";

  for my $gene (@{$chr_genes{$chr}})
  {
    my $strand = $gene->{'strand'};

    # Increment bins
    my ($upstream_dist, $downstream_dist);
    
    if($strand eq "+")
    {
      $upstream_dist = $gene->{'start'} - 1;
      $downstream_dist = $chr_length - $gene->{'end'};
    }
    else
    {
      $upstream_dist = $chr_length - $gene->{'end'};
      $downstream_dist = $gene->{'start'} - 1;
    }

    #print "upstream: $upstream_dist; downstream: $downstream_dist;\n";

    # Full bins upstream
    my $full_bins_upstream = int($upstream_dist/$bin_size);

    for(my $i = 0; $i < $full_bins_upstream; $i++) {
      $upstream_denominator[$i]++;
    }
    
    my $part_upstream_bin = ($upstream_dist/$bin_size) - $full_bins_upstream;
    
    if($part_upstream_bin > 0)
    {
      # Add part upstream bin
      $upstream_denominator[$full_bins_upstream] += $part_upstream_bin;
    }

    # Full bins downstream
    my $full_bins_downstream = int($downstream_dist/$bin_size);

    for(my $i = 0; $i < $full_bins_downstream; $i++) {
      $downstream_denominator[$i]++;
    }
  
    my $part_downstream_bin = ($downstream_dist/$bin_size) -
                              $full_bins_downstream;
  
    if($part_downstream_bin > 0)
    {
      # Add part downstream bin
      $downstream_denominator[$full_bins_downstream] += $part_downstream_bin;
    }
    
    #
    # Iterate through positions
    #
    for my $pos (@positions)
    {
      my $bin = -1;

      if($pos < $gene->{'start'})
      {
        my $dist = $gene->{'start'} - $pos - 1;
        my $bin = int($dist / $bin_size);

        #print "$pos in bin $bin ".($strand eq "+"?"up":"down")."stream\n";

        if($strand eq "+") {
          $upstream_bins[$bin]++;
        }
        else {
          $downstream_bins[$bin]++;
        }
      }
      elsif($pos > $gene->{'end'})
      {
        my $dist = $pos - $gene->{'end'} - 1;
        my $bin = int($dist / $bin_size);
        
        #print "$pos in bin $bin ".($strand eq "+"?"down":"up")."stream\n";
        
        if($strand eq "+") {
          $downstream_bins[$bin]++;
        }
        else {
          $upstream_bins[$bin]++;
        }
      }
      else
      {
        # Tweak using tx_start/end and cds_start/end at the top of this file
        # (can be outside cds but in exons... don't understand ensGene.txt)

        # In gene
        my $gene_width = $gene->{'end'} - $gene->{'start'} + 1;
        my $dist_gene = ($strand eq "+") ? $pos - $gene->{'start'}
                                         : $gene->{'end'} - $pos;

        my $gene_dist_prop = $gene_width == 1 ? 0.5 : $dist_gene / ($gene_width-1);
        push(@pos_in_genes, $gene_dist_prop);
        
        # Check exons/introns
        my $exon_starts = $gene->{'exon_starts'};
        my $exon_ends = $gene->{'exon_ends'};
        my $num_of_exons = @$exon_starts;
        
        my $in_exon = 0;

        while($in_exon < $num_of_exons && $pos > $exon_ends->[$in_exon])
        {
          $in_exon++;
        }

        # $in_exon is the index of the exon or the intron just before it
        # that the variant is in

        if($in_exon == 0 && $pos < $exon_starts->[$in_exon])
        {
          # Before first exon
          # dist starting from 0
          my $dist = ($strand eq "+") ? $pos - $gene->{'start'}
                                      : $exon_starts->[0] - $pos - 1;

          my $width = $exon_starts->[0] - $gene->{'start'};
          my $dist_prop = $width == 1 ? 0.5 : $dist / ($width-1);

          if($strand eq "+") {
            push(@pre_exons, $dist_prop);
            #print "In preexon: $pos ($dist_prop; $dist / $width)\n";
          }
          else {
            push(@post_exons, $dist_prop);
            #print "In postexon: $pos ($dist_prop; $dist / $width)\n";
          }
        }
        elsif($in_exon == $num_of_exons)
        {
          # Past last exon
          # dist starting from 0
          my $dist = ($strand eq "+") ? $pos - $exon_ends->[$in_exon-1] - 1
                                      : $gene->{'end'} - $pos;

          my $width = $gene->{'end'} - $exon_ends->[$in_exon-1];
          my $dist_prop = $width == 1 ? 0.5 : $dist / ($width-1);

          if($strand eq "+") {
            push(@post_exons, $dist_prop);
            #print "In postexon: $pos ($dist_prop; $dist / $width)\n";
          }
          else {
            push(@pre_exons, $dist_prop);
            #print "In preexon: $pos ($dist_prop; $dist / $width)\n";
          }
        }
        elsif($pos >= $exon_starts->[$in_exon])
        {
          # In Exon
          # dist starting from 0
          my $dist = ($strand eq "+") ? $pos - $exon_starts->[$in_exon]
                                      : $exon_ends->[$in_exon] - $pos;

          my $exon_length = $exon_ends->[$in_exon] - $exon_starts->[$in_exon] + 1;
          my $dist_prop = $exon_length == 1 ? 0.5 : $dist / ($exon_length-1);
          push(@exon_distances, $dist_prop);

          #print "In exon: $pos ($dist_prop; $dist / $exon_length)\n";
        }
        else
        {
          # In Intron (just before exon $in_exon)
          # dist starting from 0
          my $dist = ($strand eq "+") ? $pos - ($exon_ends->[$in_exon-1] + 1)
                                      : $exon_starts->[$in_exon] - $pos - 1;

          my $intron_length = $exon_starts->[$in_exon] - ($exon_ends->[$in_exon-1] + 1);
          my $dist_prop = $intron_length == 1 ? 0.5 : $dist / ($intron_length-1);
          push(@intron_distances, $dist_prop);

          #print "In intron: $pos ($dist_prop; $dist / $intron_length)\n";
        }
      }
    }
  }
}
