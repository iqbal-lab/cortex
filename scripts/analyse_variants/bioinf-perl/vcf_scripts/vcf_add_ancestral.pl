#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use Fcntl qw(SEEK_CUR);

use VCFFile;
use RefGenome;
use UsefulModule;

use List::Util qw(min max);

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR
"Usage: ./vcf_add_ancestral.pl <stampy.sam> <min_MAPQ> <outgroup_name> " .
  "<in.vcf> [ancestral_ref.fa1 ..]
  Assigns ancestral allele to VCF.  Adds AA and AALEN INFO tags. If <in.vcf>
  is '-' reads from stdin.

  Ancestral reference genome is optional - if it is not given, only variants
  with SVLEN != 0 will be polarised.  

  For variants with both alleles the same length, uses ancestral reference.
  For indels, uses mapping of 5'+allele+3' flank to outgroup genome.
  * AA = [. -> unknown, 0 -> ref, 1 -> alt allele]
  * AALEN = (derived length - ancestral length)\n";
  exit;
}

## Test for filtering
my $failed_vars_out = undef;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  open($failed_vars_out, ">-");
}
##

if(@ARGV < 4)
{
  print_usage();
}

# Ancestral Reference from:
# ftp://ftp.ebi.ac.uk/pub/databases/ensembl/jherrero/ancestral/
#
# It uses the same coordinates as PanTro2 / hg19
#
# Contains:
# ACGT => confident
# acgt => less confident
# . => unknown
# - => deleted

my $mapping_file = shift;
my $min_mapq = shift;
my $outgroup_name = shift;
my $vcf_file = shift;
my @ancestral_ref_files = @ARGV;

# Commandline arg checks
if($min_mapq !~ /^\d+$/)
{
  print_usage("Min. MAPQ must be a positive integer (>=0)");
}

#
# Open VCF Handle
#
my $vcf_handle;

if(defined($vcf_file) && $vcf_file ne "-")
{
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

my $vcf = new VCFFile($vcf_handle);

# Print non-PASS variants straight to stdout if -p passed
if(defined($failed_vars_out)) { $vcf->set_filter_failed($failed_vars_out);}

#
# Open mapping file
#
my $mapping_handle;

open($mapping_handle, $mapping_file)
  or die("Cannot open mapping file '$mapping_file'");

# skip header lines starting @...
my $mapping_line;
while(defined($mapping_line = <$mapping_handle>) && $mapping_line =~ /^@/) {}

seek($mapping_handle,-length($mapping_line),SEEK_CUR);

my @mapping_columns = qw(name flags chr pos MAPQ CIGAR
                         rnext pnext tlen seq quality);

#
# Load Ancestral Reference
#
print STDERR "vcf_add_ancestral_to_indels.pl: Loading ancestral ref...\n";

my $ancestral;

if(@ancestral_ref_files > 1 ||
   (@ancestral_ref_files > 0 && $ancestral_ref_files[0] ne "-"))
{
  $ancestral = new RefGenome(uppercase => 1);
  $ancestral->load_from_files(@ancestral_ref_files);
}

#
# Print VCF header
#
my $tag_description = "Ancestral Allele using outgroup $outgroup_name " .
                      "(.=unknown, 0=reference, 1=alternate allele)";

$vcf->add_header_tag("INFO", "AA", 1, "String", $tag_description);
$vcf->print_header();

# Read VCF entry
my $vcf_entry;

my $num_of_variants = 0;
my $num_of_np = 0;
my $num_of_indels = 0;

my $num_of_np_polarised = 0;
my $num_of_np_ref = 0;

my $num_of_indels_polarised = 0;
my $num_of_indels_ref = 0;

while(defined($vcf_entry = $vcf->read_entry()))
{
  $num_of_variants++;

  my $var_id = $vcf_entry->{'ID'};
  my $chr = $vcf_entry->{'CHROM'};
  my $pos = $vcf_entry->{'true_POS'};
  my $svlen = $vcf_entry->{'INFO'}->{'SVLEN'};

  my $ref_allele = uc($vcf_entry->{'true_REF'});
  my $alt_allele = uc($vcf_entry->{'true_ALT'});

  my $aa = '.';

  if($svlen == 0)
  {
    # Nucleotide Polymorphism (NP) - use ancestral allele
    $num_of_np++;

    if(defined($ancestral))
    {
      # $pos-1 because perl uses 0-based
      my $anc_ref = $ancestral->get_chr_substr($chr, $pos-1, length($ref_allele));
      my $anc_alt = $ancestral->get_chr_substr($chr, $pos-1, length($alt_allele));

      if(!defined($anc_ref))
      {
        print STDERR "vcf_add_ancestral.pl Error: " .
                     "Ancestor lacks chromosome '$chr'\n";
        die();
      }

      my $ref_matches = ancestral_match($anc_ref, $ref_allele);
      my $alt_matches = ancestral_match($anc_alt, $alt_allele);

      if($ref_matches && !$alt_matches)
      {
        $aa = '0';
        $num_of_np_polarised++;
      }
      elsif(!$ref_matches && $alt_matches)
      {
        $aa = '1';
        $num_of_np_polarised++;
        $num_of_np_ref++;
      }
    }
  }
  else
  {
    # Indel - use stampy mapping
    $num_of_indels++;

    # Read corresponding mapping
    my ($ref_mapping_line, $alt_mapping_line);

    my %ref_mapping;
    my %alt_mapping;

    while(defined($ref_mapping_line = <$mapping_handle>) &&
          defined($alt_mapping_line = <$mapping_handle>))
    {
      chomp($ref_mapping_line);
      chomp($alt_mapping_line);

      if(!defined($ref_mapping_line) || !defined($alt_mapping_line))
      {
        die("Ran out of mappings (var: $var_id)");
      }

      %ref_mapping = ();
      %alt_mapping = ();

      @ref_mapping{@mapping_columns} = split(/\t/, $ref_mapping_line);
      @alt_mapping{@mapping_columns} = split(/\t/, $alt_mapping_line);

      if(lc($ref_mapping{'name'}) eq lc($var_id."_ref"))
      {
        last;
      }
    }

    if(!defined($ref_mapping_line))
    {
      die("vcf_add_ancestral.pl: Couldn't find mapping for var '$var_id'\n");
    }
    elsif(lc($alt_mapping{'name'}) ne lc($var_id."_alt"))
    {
      die("Alt alleles' names don't match " .
          "(var '$var_id' => '$ref_mapping{'name'}')");
    }

    if($ref_mapping{'CIGAR'} eq "*" || $alt_mapping{'CIGAR'} eq "*" ||
       $ref_mapping{'chr'} ne $alt_mapping{'chr'} || 
       abs($ref_mapping{'pos'} - $alt_mapping{'pos'}) > 32)
    {
      # Ref/alt don't map, map to different chromosomes or map too far apart
      # 32 used here because it is the flank size
      $aa = '.';
    }
    else
    {
      if($ref_mapping{'CIGAR'} !~ /^(?:\d+[MID])+$/i)
      {
        # Full range MIDNSHP=X
        die("Unexpected cigar entry " .
            "(var: $var_id; cigar: '$ref_mapping{'CIGAR'}')");
      }
      elsif($alt_mapping{'CIGAR'} !~ /^(?:\d+[MID])+$/i)
      {
        die("Unexpected cigar entry " .
            "(var: $var_id; cigar: '$alt_mapping{'CIGAR'}')");
      }

      my $long_cigar = $ref_mapping{'CIGAR'};
      my $short_cigar = $alt_mapping{'CIGAR'};

      my $long_mapq = $ref_mapping{'MAPQ'};
      my $short_mapq = $alt_mapping{'MAPQ'};

      if($svlen > 0)
      {
        # alt is longer than ref -> swap
        ($short_cigar, $long_cigar) = ($long_cigar, $short_cigar);
        ($short_mapq, $long_mapq) = ($long_mapq, $short_mapq);
      }

      $aa = get_indel_ancestor($long_cigar, $short_cigar,
                               $long_mapq, $short_mapq,
                               $svlen);
      
      if($aa ne '.')
      {
        $num_of_indels_polarised++;

        if($aa eq '0')
        {
          $num_of_indels_ref++;
        }
      }
    }
  }

  $vcf_entry->{'INFO'}->{'AA'} = $aa;

  if($aa eq '.')
  {
    $vcf_entry->{'INFO'}->{'AALEN'} = '.';
  }
  else
  {
    $vcf_entry->{'INFO'}->{'AALEN'} = ($aa eq '0' ? $svlen : -$svlen);
  }

  $vcf->print_entry($vcf_entry);

  # DEBUG
  #if($num_of_variants > 4) {
  #  die("DEBUG END");
  #}
}

print STDERR "vcf_add_ancestral.pl: " .
             pretty_fraction($num_of_np_polarised, $num_of_np) . " " .
             "NP polarised\n";

print STDERR "vcf_add_ancestral.pl: of which " .
             pretty_fraction($num_of_np_ref, $num_of_np_polarised) . " " .
             "were ref\n";

print STDERR "vcf_add_ancestral.pl: " .
             pretty_fraction($num_of_indels_polarised, $num_of_indels) . " " .
             "indels polarised\n";

print STDERR "vcf_add_ancestral.pl: of which " .
             pretty_fraction($num_of_indels_ref, $num_of_indels_polarised) .
             " were ref\n";

my $total_polarised = $num_of_np_polarised + $num_of_indels_polarised;

print STDERR "vcf_add_ancestral.pl: Total variants polarised = " .
             pretty_fraction($total_polarised, $num_of_variants) . "\n";

# Done
close($vcf_handle);
close($mapping_handle);


#
# Methods
#

sub get_indel_ancestor
{
  my ($long_cigar, $short_cigar, $longer_mapq, $shorter_mapq, $svlen) = @_;
  my $size = abs($svlen);

  my $long_cigar_match = ($long_cigar =~ /^\d+M$/i);
  my $short_cigar_match = ($short_cigar =~ /^\d+M$/i);

  if($short_cigar_match && !$long_cigar_match && $shorter_mapq >= $min_mapq)
  {
    # shorter allele matches
    # If (alt-ref) > 0 then reference (0) is the shorter allele
    #                  otherwise alt (1)
    return $svlen > 0 ? '0' : '1';
  }
  elsif($long_cigar_match && !$short_cigar_match && $longer_mapq >= $min_mapq)
  {
    # longer allele matches
    # If (alt-ref) > 0 then alt (1) is the longer allele
    #                  otherwise ref (0)
    return $svlen > 0 ? '1' : '0';
  }
  else
  {
    return '.';
  }
}

sub ancestral_match
{
  my ($anc, $allele) = @_;

  $anc = uc($anc);
  $allele = uc($allele);

  my $len = length($anc);
  my $unknowns = 0;

  for(my $i = 0; $i < $len; $i++)
  {
    my $anc_char = substr($anc,$i,1);

    if($anc_char eq ".")
    {
      $unknowns++;
    }
    elsif($anc_char ne substr($allele,$i,1))
    {
      return 0;
    }
  }

  return $unknowns < ($len / 2);
}

