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
"Usage: ./vcf_correct_strand.pl [OPTIONS] <in.vcf> <ref1.fa ..>
  If the REF allele doesn't match the reference but the ALT allele does, switch
  them and correct the genotypes.
  
  --remove_mismatches       Remove mismatches that cannot be fixed
  --flag_mismatches <tag>   Add an INFO flag to variants that don't match REF
  --filter_mismatches <tag> Add a FILTER flag to variants that don't match REF
  --flag_swapped <tag>      Add an INFO flag to variants that were swapped
  --filter_swapped <tag>    Add a FILTER flag to variants that were swapped

  --keep_phased   Keep phasing when switching ref/alt (e.g. 1|1 -> 0|0) [default]
  --remove_phased Variants lose phasing if switched (e.g. 1|1 -> 0/0)\n";

  exit;
}

## Test for filtering
my $failed_vars_out = undef;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  open($failed_vars_out, ">-");
}
##

if(@ARGV < 2)
{
  print_usage();
}

my $remove_ref_mismatch = 0;
my ($flag_mismatches, $flag_swapped);
my ($filter_mismatches, $filter_swapped);

my $keep_phased = 1;
my ($set_keep_phased, $set_remove_phased);

while(@ARGV > 0)
{
  if($ARGV[0] =~ /^-?-remove_mismatches$/i)
  {
    shift;
    $remove_ref_mismatch = 1;
  }
  elsif($ARGV[0] =~ /^-?-keep_phased?$/i)
  {
    shift;
    $keep_phased = 1;
    $set_keep_phased = 1;
  }
  elsif($ARGV[0] =~ /^-?-remove_phased?$/i)
  {
    shift;
    $keep_phased = 0;
    $set_remove_phased = 1;
  }
  elsif($ARGV[0] =~ /^-?-flag_mismatches$/i)
  {
    shift;

    if(!defined($flag_mismatches = shift))
    {
      print_usage("--flag_mismatches <tag> needs a tag name!");
    }
  }
  elsif($ARGV[0] =~ /^-?-filter_mismatches$/i)
  {
    shift;

    if(!defined($filter_mismatches = shift))
    {
      print_usage("--filter_mismatches <tag> needs a tag name!");
    }
  }
  elsif($ARGV[0] =~ /^-?-flag_swapped$/i)
  {
    shift;

    if(!defined($flag_swapped = shift))
    {
      print_usage("--flag_swapped <tag> needs a tag name!");
    }
  }
  elsif($ARGV[0] =~ /^-?-filter_swapped$/i)
  {
    shift;

    if(!defined($filter_swapped = shift))
    {
      print_usage("--filter_swapped <tag> needs a tag name!");
    }
  }
  elsif(substr($ARGV[0],0,1) eq "-")
  {
    print_usage("Unknown option '$ARGV[0]'");
  }
  else
  {
    last;
  }
}

my $vcf_file = shift;
my @ref_files = @ARGV;

if(@ref_files == 0)
{
  print_usage("Not enough arguments.");
}

if(defined($set_keep_phased) && defined($set_remove_phased))
{
  print_usage("Cannot use --keep_phased and --remove_phased together");
}

if($remove_ref_mismatch)
{
  if(defined($flag_mismatches))
  {
    print_usage("Cannot set --remove_ref_mismatch and --flag_mismatches... " .
                "it just doesn't make sense!");
  }
  elsif(defined($filter_mismatches))
  {
    print_usage("Cannot set --remove_ref_mismatch and --filter_mismatches... " .
                "it just doesn't make sense!");
  }
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
elsif(-p STDIN)
{
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

$vcf->print_header();

my @samples = $vcf->get_list_of_sample_names();

my $vcf_entry;

my $num_of_variants = 0;
my $num_of_ref_mismatch = 0;
my $num_of_switched = 0;

my %missing_chrs = ();
my %format_fields = ();

while(defined($vcf_entry = $vcf->read_entry()))
{
  $num_of_variants++;

  my $chr = $vcf_entry->{'CHROM'};

  # Get coordinates, convert to zero based
  my $var_start = $vcf_entry->{'true_POS'}-1;

  my $ref_allele = uc($vcf_entry->{'true_REF'});
  my $alt_allele = uc($vcf_entry->{'true_ALT'});

  my $max_allele_len = max(length($ref_allele), length($alt_allele));

  my $unfixed_mismatch = 0;

  my $genome_seq;

  if(!$genome->chr_exists($chr))
  {
    $missing_chrs{$chr} = 1;
  }
  else
  {
    $genome_seq = $genome->get_chr_substr($chr, $var_start, $max_allele_len);

    if($genome->get_chr_length($chr) < $var_start + $max_allele_len)
    {
      print STDERR "vcf_correct_strand.pl: var " . $vcf_entry->{'ID'} . " at " .
                   "$chr:$var_start:$max_allele_len is out of " .
                   "bounds of $chr:1:" . $genome->get_chr_length($chr) . "\n";
      $genome_seq = undef;
    }
  }

  if(defined($genome_seq))
  {
    my $genome_seq_ref = uc(substr($genome_seq, 0, length($ref_allele)));
    my $genome_seq_alt = uc(substr($genome_seq, 0, length($alt_allele)));

    if($ref_allele ne $genome_seq_ref)
    {
      # Mismatch
      $num_of_ref_mismatch++;

      if($alt_allele eq $genome_seq_alt)
      {
        # Switch!
        $num_of_switched++;

        if(defined($flag_swapped))
        {
          $vcf_entry->{'INFO_flags'}->{$flag_swapped} = 1;
        }

        if(defined($filter_swapped))
        {
          vcf_add_filter_txt($vcf_entry, $filter_swapped);
        }

        # Switch alleles
        my $tmp = $vcf_entry->{'REF'};
        $vcf_entry->{'REF'} = $vcf_entry->{'ALT'};
        $vcf_entry->{'ALT'} = $tmp;

        # Switch genotypes
        if(!defined($vcf_entry->{'FORMAT'}))
        {
          die("No format column");
        }

        my @formats = @{$vcf_entry->{'FORMAT'}};
        @format_fields{@formats} = 1;

        for my $sample (@samples)
        {
          if(defined($vcf_entry->{$sample}->{'GT'}))
          {
            if($vcf_entry->{$sample}->{'GT'} =~ /^([.01])([\/\|])([.01])$/)
            {
              my ($gt1,$gt_sep,$gt2) = ($1,$2,$3);

              # Switch genotype [01] => [10], . => .
              my %switch_gt = ('1' => '0', '0' => '1', '.' => '.');

              if(!$keep_phased)
              {
                # Change separator to / (i.e. not genotyped anymore)
                $gt_sep = '/';
              }

              $vcf_entry->{$sample}->{'GT'}
                = $switch_gt{"$gt1"}.$gt_sep.$switch_gt{"$gt2"};
            }
            else
            {
              die("Invalid GT format '".$vcf_entry->{$sample}->{'GT'}."' " .
                  "[var:$vcf_entry->{'ID'}; sample:$sample]");
            }
          }

          if(defined($vcf_entry->{$sample}->{'COV'}))
          {
            if($vcf_entry->{$sample}->{'COV'} =~ /^(\d+),(\d+)$/)
            {
              $vcf_entry->{$sample}->{'COV'} = $2.",".$1;
            }
            else
            {
              die("Invalid COV format [var:$vcf_entry->{'ID'}; sample:$sample]");
            }
          }
        }
      }
      else
      {
        # Couldn't fix mismatch
        $unfixed_mismatch = 1;

        if(defined($flag_mismatches))
        {
          $vcf_entry->{'INFO_flags'}->{$flag_mismatches} = 1;
        }

        if(defined($filter_mismatches))
        {
          vcf_add_filter_txt($vcf_entry, $filter_mismatches);
        }
      }
    }
  }

  if(!$unfixed_mismatch || !$remove_ref_mismatch)
  {
    $vcf->print_entry($vcf_entry);
  }
}

my @format_field_names
  = grep {uc($_) ne "GT" && uc($_) ne "COV"}
    sort keys %format_fields;

if(@format_field_names > 0)
{
  print STDERR "vcf_correct_strand.pl: Ignored format fields: " .
               join(",", @format_field_names) . "\n";
  print STDERR "                       (I only manipulate GT & COVG)\n";
}

my @missing_chr_names = sort keys %missing_chrs;

if(@missing_chr_names > 0)
{
  print STDERR "vcf_correct_strand.pl Warning: Missing chromosomes: " .
               join(", ", @missing_chr_names) . "\n";
}

print STDERR "vcf_correct_strand.pl: " .
             pretty_fraction($num_of_ref_mismatch, $num_of_variants) . " " .
             "variants did not match the reference\n";

print STDERR "vcf_correct_strand.pl: " .
             pretty_fraction($num_of_switched,
                             $num_of_ref_mismatch) . " " .
             "mismatches were switched successfully\n";

my $num_of_remaining_mismatches = $num_of_ref_mismatch - $num_of_switched;

print STDERR "vcf_correct_strand.pl: " .
             pretty_fraction($num_of_remaining_mismatches,
                             $num_of_ref_mismatch) . " " .
             "mismatches couldn't be fixed and were " .
             ($remove_ref_mismatch ? "removed" : "left in") . "\n";

if(defined($flag_mismatches))
{
  print STDERR "vcf_correct_strand.pl: variants that don't match the ref genome " .
               "were labelled with '$flag_mismatches' INFO flag\n";
}

if(defined($filter_mismatches))
{
  print STDERR "vcf_correct_strand.pl: variants that don't match the ref genome " .
               "were labelled with FILTER '$filter_mismatches'\n";
}

if(defined($flag_swapped))
{
  print STDERR "vcf_correct_strand.pl: variants that had alleles swapped were " .
               "labelled with '$flag_swapped' INFO flag\n";
}

if(defined($filter_swapped))
{
  print STDERR "vcf_correct_strand.pl: variants that had alleles swapped were " .
               "labelled with FILTER '$filter_swapped'\n";
}

close($vcf_handle);
