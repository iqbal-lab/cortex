#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(min max);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use UsefulModule; # num2str
use VCFFile;
use RefGenome;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./vcf_align.pl [--tag <name>|--remove_ref_mismatch] <LEFT|RIGHT> " .
  "<in.vcf> <ref1.fa ..>
  Shift clean indel variants to the left/right.  Variants that do not match
  the reference and those that are not clean indels are printed unchanged.
  FASTA entry names must match VCF CHROM column.  If <in.vcf> is '-', reads
  from STDIN.

  --tag <tag_name>       INFO tag to label variation in position
  --remove_ref_mismatch  Remove variants that do not match the reference\n";

  exit;
}

## Test for filtering
my $failed_vars_out = undef;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  open($failed_vars_out, ">-");
}
##

if(@ARGV < 3)
{
  print_usage();
}

my $tag;
my $remove_ref_mismatch = 0;

while(@ARGV > 3)
{
  if($ARGV[0] =~ /^-?-tag$/i)
  {
    shift;
    $tag = shift;
  }
  elsif($ARGV[0] =~ /^-?-remove_ref_mismatch$/i)
  {
    shift;
    $remove_ref_mismatch = 1;
  }
  else
  {
    last;
  }
}

my $justify = lc(shift);
my $vcf_file = shift;
my @ref_files = @ARGV;

if(@ref_files == 0)
{
  print_usage("Not enough arguments.");
}

if($justify ne "left" && $justify ne "right")
{
  print_usage("neither 'left' nor 'right' given");
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
# Load reference
#
my $genome = new RefGenome(uppercase => 1);
$genome->load_from_files(@ref_files);

#
# Read VCF
#
my $vcf = new VCFFile($vcf_handle);

# Print non-PASS variants straight to stdout if -p passed
if(defined($failed_vars_out)) { $vcf->set_filter_failed($failed_vars_out);}

# Add justify info to header and print
$vcf->add_header_metainfo("variants_justified", $justify);

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
my $num_of_variants_moved = 0;
my $num_of_clean_indels_matching_ref = 0;

my $max_position_diff = 0;
my $total_position_diff = 0;

my %missing_chrs = ();

while(defined($vcf_entry = $vcf->read_entry()))
{
  $num_of_variants++;

  my $print = 1;

  my $chr = $vcf_entry->{'CHROM'};

  # Get coordinates, convert to zero based
  my $var_start = $vcf_entry->{'true_POS'}-1;
  my $ref_allele = uc($vcf_entry->{'true_REF'});
  my $alt_allele = uc($vcf_entry->{'true_ALT'});

  # Get inserted or deleted sequence
  my $indel = vcf_get_clean_indel($vcf_entry);

  my $ref_length = $genome->get_chr_length($chr);

  if(!defined($ref_length))
  {
    $missing_chrs{$chr} = 1;
  }
  elsif($ref_length < $var_start+length($ref_allele))
  {
    my $ref_chrom_name = $genome->guess_chrom_fasta_name($chr);
    my $real_pos = $vcf_entry->{'POS'};

    print STDERR "vcf_align.pl Warning: var " . $vcf_entry->{'ID'} . " at " .
                  "$chr:$real_pos:" . length($ref_allele) . " is out of " .
                  "bounds of $ref_chrom_name:1:" . $ref_length . "\n";
  }
  elsif($genome->get_chr_substr($chr, $var_start, length($ref_allele)) ne
        $ref_allele)
  {
    $num_ref_mismatch++;

    if($remove_ref_mismatch)
    {
      $print = 0;
    }
  }
  elsif(defined($indel))
  {
    $indel = uc($indel);

    $num_of_clean_indels_matching_ref++;

    # align to the left
    my ($pos_left, $indel_left) = get_left_aligned_position($chr,
                                                            $var_start,
                                                            $indel);

    my $ref_allele_len = length($ref_allele);

    # align to the right
    my ($pos_right, $indel_right) = get_right_aligned_position($chr,
                                                               $ref_length,
                                                               $var_start,
                                                               $ref_allele_len,
                                                               $indel);

    my $position_diff = $pos_right - $pos_left;

    if(defined($tag) && $position_diff > 0)
    {
      $vcf_entry->{'INFO'}->{$tag} = $position_diff;
    }

    $max_position_diff = max($position_diff, $max_position_diff);
    $total_position_diff += $position_diff;

    my ($new_pos, $new_indel);

    if($justify eq "left")
    {
      $new_pos = $pos_left;
      $new_indel = $indel_left;
    }
    else
    {
      $new_pos = $pos_right;
      $new_indel = $indel_right;
    }

    # Update VCF entry values
    if($var_start != $new_pos)
    {
      $num_of_variants_moved++;

      # $new_pos is 0-based, VCF POS is 1-based
      $vcf_entry->{'true_POS'} = $new_pos+1;
      $vcf_entry->{'POS'} = $new_pos;
    
      my $indel_allele = length($ref_allele) > 0 ? 'true_REF' : 'true_ALT';
      $vcf_entry->{$indel_allele} = $new_indel;

      my $prior_base = $genome->get_chr_substr($chr, $new_pos-1, 1);
      $vcf_entry->{'REF'} = $prior_base.$vcf_entry->{'true_REF'};
      $vcf_entry->{'ALT'} = $prior_base.$vcf_entry->{'true_ALT'};
    }
  }

  if($print)
  {
    $vcf->print_entry($vcf_entry);
  }
}

my @missing_chr_names = sort keys %missing_chrs;

if(@missing_chr_names > 0)
{
  print STDERR "vcf_align.pl Warning: Missing chromosomes: " .
               join(", ", @missing_chr_names) . "\n";
}

if($num_ref_mismatch > 0)
{
  print STDERR "vcf_align.pl: " .
               pretty_fraction($num_ref_mismatch, $num_of_variants) . " " .
               "variants did not match the reference " .
               ($remove_ref_mismatch ? "(removed)" : "(left in)") . "\n";
}

print STDERR "vcf_align.pl: " .
             pretty_fraction($num_of_clean_indels_matching_ref,
                             $num_of_variants) . " " .
             "of variants are clean indels\n";

print STDERR "vcf_align.pl: " .
             pretty_fraction($num_of_variants_moved, $num_of_variants) . " " .
             "variants moved to the $justify\n";

print STDERR "vcf_align.pl: average ambiguity per variant: " .
             sprintf("%.3f", $total_position_diff / $num_of_variants) .
             " bp\n";

if($num_of_clean_indels_matching_ref > 0)
{
  my $clean_indel_rate
    = $total_position_diff / $num_of_clean_indels_matching_ref;

  print STDERR "vcf_align.pl: average ambiguity per clean indel matching " .
               "the ref: " . sprintf("%.3f", $clean_indel_rate) . " bp\n";
}

print STDERR "vcf_align.pl: maximum variant ambiguity: " .
             $max_position_diff . " bp\n";

close($vcf_handle);

# $pos is 0-based
sub get_left_aligned_position
{
  my ($ref_chr, $pos, $indel) = @_;

  # align to the left
  # while base before variant (on the reference) equals last base of indel
  while($pos > 0 &&
        $genome->get_chr_substr($ref_chr, $pos-1, 1) eq substr($indel, -1))
  {
    $indel = substr($indel, -1) . substr($indel, 0, -1);
    $pos--;
  }

  return ($pos, $indel);
}

# $pos is 0-based
sub get_right_aligned_position
{
  my ($ref_chr, $chr_length, $pos, $ref_allele_length, $indel) = @_;

  # align to the right
  # while base after variant (on the reference) equals first base of indel
  while($pos < $chr_length &&
        $genome->get_chr_substr($ref_chr, $pos+$ref_allele_length, 1) eq
        substr($indel, 0, 1))
  {
    $indel = substr($indel, 1) . substr($indel, 0, 1);
    $pos++;
  }

  return ($pos, $indel);
}
