#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(min max);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use GeneticsModule;
use VCFFile;
use FASTNFile;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./vcf_print_allele_branches.pl <K> <in.vcf> [ref.files ..]
  Prints alleles plus flanking sequence from VCF file in FASTA format. If
  <in.vcf> is '-', reads from stdin.  If [ref.files] omitted, uses INFO flank
  tags.

  <K> is flank size

  OTPIONS:
    --trim <t>          trim sequences longer than <t>\n";

  exit;
}

if(@ARGV < 2 || @ARGV > 5)
{
  print_usage();
}

my $trim;

while(@ARGV > 2 && $ARGV[0] =~ /^--/)
{
  my $arg = shift(@ARGV);

  if($arg =~ /^-?-trim$/i)
  {
    $trim = shift(@ARGV);

    if($trim !~ /^\d+$/ || $trim == 0)
    {
      print_usage("Invalid trim value ('$trim') - " .
                  "must be positive integer (>0)");
    }
  }
  else
  {
    last;
  }
}

my $max_flank_size = shift;
my $vcf_file = shift;
my @ref_files = @ARGV;

my $vcf_handle;

if($vcf_file ne "-") {
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

my %ref_genomes;

if(@ref_files > 0)
{
  # Load reference
  my ($ref_genomes_hashref) = read_all_from_files(@ref_files);

  %ref_genomes = %$ref_genomes_hashref;
}

#
# Read VCF
#
my $vcf = new VCFFile($vcf_handle);

if(@ref_files > 0)
{
  # Set the reference genome chromosome names
  $vcf->set_ref_chrom_names(keys %ref_genomes);
}

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  my $ref_allele = $vcf_entry->{'true_REF'};
  my $alt_allele = $vcf_entry->{'true_ALT'};

  my $left_flank = "";
  my $right_flank = "";

  if($max_flank_size > 0)
  {
    if(@ref_files == 0)
    {
      if(!defined($vcf_entry->{'INFO'}->{'left_flank'}) ||
         !defined($vcf_entry->{'INFO'}->{'right_flank'}))
      {
        print_usage("Flanks missing on entry $vcf_entry->{'ID'}");
      }

      $left_flank = $vcf_entry->{'INFO'}->{'left_flank'};
      $right_flank = $vcf_entry->{'INFO'}->{'right_flank'};
  
      if(length($left_flank) > $max_flank_size)
      {
        $left_flank = substr($left_flank, 0, $max_flank_size);
      }
    
      if(length($left_flank) > $max_flank_size)
      {
        $right_flank = substr($right_flank, -$max_flank_size);
      }
    }
    else
    {
      my $chr = $vcf_entry->{'CHROM'};
      my $ref_chr = $vcf->guess_ref_chrom_name($chr);

      if(defined($ref_chr))
      {
        ($left_flank, $right_flank) = get_flanks_from_ref_genome($vcf_entry,
                                                                 $ref_chr);
      }
      else
      {
        die("Chromosome '$vcf_entry->{'CHROM'}' not in reference " .
            "(".join(",",sort keys %ref_genomes).")");
      }
    }
  }

  print_to_fasta($vcf_entry->{'ID'}."_ref",
                 $left_flank . $ref_allele . $right_flank);

  print_to_fasta($vcf_entry->{'ID'}."_alt",
                 $left_flank . $alt_allele . $right_flank);
}

close($vcf_handle);

sub print_to_fasta
{
  my ($fasta_name,$fasta_seq) = @_;

  print ">$fasta_name\n";
  
  if(defined($trim))
  {
    print substr($fasta_seq,0,$trim) . "\n";
  }
  else
  {
    print "$fasta_seq\n";
  }
}

sub get_flanks_from_ref_genome
{
  my ($vcf_entry, $ref_chr) = @_;

  # Convert from 1-based to 0-based
  my $var_start = $vcf_entry->{'true_POS'}-1;
  
  my $var_length = length($vcf_entry->{'true_REF'});

  my $left_flank_start = max(0, $var_start-$max_flank_size);
  my $left_flank_length = $var_start - $left_flank_start;

  my $right_flank_start = $var_start + $var_length;
  my $right_flank_length = min($max_flank_size,
                               length($ref_genomes{$ref_chr}) - $right_flank_start);
  
  my $left_flank = substr($ref_genomes{$ref_chr},
                          $left_flank_start, $left_flank_length);
  my $right_flank = substr($ref_genomes{$ref_chr},
                           $right_flank_start, $right_flank_length);

  return ($left_flank, $right_flank);
}
