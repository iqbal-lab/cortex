#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(min max);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use GeneticsModule;
use UsefulModule; # num2str
use VCFFile;
use FASTNFile;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./vcf_align.pl [--tag <name>] <LEFT|RIGHT> <in.vcf> <ref1.fa ..>\n" .
"  Shift clean indel variants to the left/right.  Variants that do not match\n".
"  the reference and those that are not clean indels are printed unchanged.\n" .
"  FASTA entry names must match VCF CHROM column.  If <in.vcf> is '-', reads\n".
"  from STDIN.\n".
"  --tag <tag_name>   INFO tag to label variation in position.\n";
  exit;
}

if(@ARGV < 3)
{
  print_usage();
}

my $tag;

if($ARGV[0] =~ /^-?-tag$/i)
{
  shift;
  $tag = shift;
}

my $justify = lc(shift);
my $vcf_file = shift;
my @ref_files = @ARGV;

if($justify ne "left" && $justify ne "right")
{
  print_usage("neither 'left' nor 'right' given");
}

#
# Open VCF Handle
#
my $vcf_handle;

if($vcf_file ne "-")
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
# Load reference
#
my ($ref_genomes_hashref, $qual) = read_all_from_files(@ref_files);

my %ref_genomes = %$ref_genomes_hashref;

# Correct '1..' to 'chr1...' etc
# Change chromosome to uppercase
while(my ($key,$value) = each(%ref_genomes))
{
  my $new_key = get_clean_chr_name($key);

  if($new_key ne $key)
  {
    $ref_genomes{$new_key} = uc($ref_genomes{$key});
    delete($ref_genomes{$key});
  }
  else {
    $ref_genomes{$key} = uc($ref_genomes{$key});
  }
}

#
# Read VCF
#
my $vcf = new VCFFile($vcf_handle);

# Add justify info to header and print
$vcf->add_header_metainfo("variants_justified",$justify);

if(defined($tag))
{
  # Add tag
  my $description = "Possible variation in clean indel position";
  $vcf->add_header_tag("INFO", $tag, 1, "Integer", $description);
}

$vcf->print_header();

my $vcf_entry;

my $num_ref_mismatch = 0;
my $num_of_variants = 0;

my %missing_chrs = ();

while(defined($vcf_entry = $vcf->read_entry()))
{
  $num_of_variants++;

  my $chr = $vcf_entry->{'CHROM'};

  # Get coordinates, convert to zero based
  my $var_start = $vcf_entry->{'true_POS'}-1;
  my $ref_allele = $vcf_entry->{'true_REF'};
  my $alt_allele = $vcf_entry->{'true_ALT'};

  # Get inserted or deleted sequence
  my $indel = get_clean_indel($vcf_entry);

  if(!defined($ref_genomes{$chr}))
  {
    $missing_chrs{$chr} = 1;
  }
  elsif(substr($ref_genomes{$chr}, $var_start, length($ref_allele)) ne
        uc($ref_allele))
  {
    $num_ref_mismatch++;
  }
  elsif(defined($indel))
  {
    my ($new_pos, $new_indel);

    if($justify eq "left")
    {
      # align to the left
      ($new_pos, $new_indel) = get_left_aligned_position($chr, $var_start,
                                                         $indel);
    }
    else
    {
      # align to the right
      ($new_pos, $new_indel) = get_right_aligned_position($chr, $var_start,
                                                          length($ref_allele),
                                                          $indel);
    }

    if(defined($tag))
    {
      my ($pos_left, $pos_right);
      
      if($justify eq "left")
      {
        $pos_left = $new_pos;
        ($pos_right) = get_right_aligned_position($chr, $var_start,
                                                  length($ref_allele), $indel);
      }
      else
      {
        ($pos_left) = get_left_aligned_position($chr, $var_start, $indel);
        $pos_right = $new_pos;
      }

      if($pos_left != $pos_right)
      {
        $vcf_entry->{'INFO'}->{$tag} = $pos_right - $pos_left;
      }
    }

    # Update VCF entry values
    # $new_pos is 0-based, VCF POS is 1-based
    $vcf_entry->{'true_POS'} = $new_pos+1;
    $vcf_entry->{'POS'} = $new_pos;
    
    $vcf_entry->{length($ref_allele) > 0 ? 'true_REF' : 'true_ALT'} = $new_indel;

    my $prior_base = substr($ref_genomes{$chr}, $new_pos-1, 1);
    $vcf_entry->{'REF'} = $prior_base.$vcf_entry->{'true_REF'};
    $vcf_entry->{'ALT'} = $prior_base.$vcf_entry->{'true_ALT'};
  }

  $vcf->print_entry($vcf_entry);
}

my @missing_chr_names = sort keys %missing_chrs;

if(@missing_chr_names > 0)
{
  print STDERR "vcf_add_flanks.pl: Missing chromosomes: " .
               join(", ", @missing_chr_names) . "\n";
}

if($num_ref_mismatch > 0)
{
  print STDERR "vcf_add_flanks.pl: " .
               pretty_fraction($num_ref_mismatch, $num_of_variants) . " " .
               "variants removed for not matching the reference\n";
}

close($vcf_handle);

# $pos is 0-based
sub get_left_aligned_position
{
  my ($chr, $pos, $indel) = @_;

  # align to the left
  # while base before variant (on the reference) equals last base of indel
  while(substr($ref_genomes{$chr}, $pos-1, 1) eq uc(substr($indel, -1)))
  {
    $indel = substr($indel, -1) . substr($indel, 0, -1);
    $pos--;
  }

  return ($pos, $indel);
}

# $pos is 0-based
sub get_right_aligned_position
{
  my ($chr, $pos, $ref_allele_length, $indel) = @_;

  # align to the right
  # while base after variant (on the reference) equals first base of indel
  while(substr($ref_genomes{$chr}, $pos+$ref_allele_length, 1) eq
        uc(substr($indel, 0, 1)))
  {
    $indel = substr($indel, 1) . substr($indel, 0, 1);
    $pos++;
  }

  return ($pos, $indel);
}
