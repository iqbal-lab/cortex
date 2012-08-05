#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(min max reduce);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use VCFFile;
use UsefulModule; # for pretty_fraction

use constant {DUPE_SELECT_NONE => 0,
              DUPE_SELECT_LOWEST_TAG => 1,
              DUPE_SELECT_HIGHEST_TAG => 2,
              DUPE_SELECT_FIRST => 3,
              DUPE_SELECT_LAST => 4,
              DUPE_SELECT_RANDOM => 5};

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./vcf_remove_dupes.pl [OPTIONS] [file.vcf]
  Remove entries that match position, REF and ALT alleles. Assumes sorted VCF.

  Order duplicates by an INFO tag and take the variant with highest/lowest value:
  --take_lowest <tag>  OR
  --take_highest <tag> OR
  --take_first         OR  [default]
  --take_last          OR
  --take_random        OR
  --take_none

  --removed_tag <tag>  Used with --take_[lowest|highest], adds the removed
                       values as a list to the selected variant

  --filter_txt <txt>   Add to / set the filter column instead of removing\n";

  exit;
}

## Test for filtering
my $failed_vars_out = undef;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  open($failed_vars_out, ">-");
}
##

if(@ARGV > 5)
{
  print_usage();
}

# How to select a single variant from a set of dupes
my $select_dupe;
# Tag to use in selection
my $select_tag;
my $list_removed_tag;

my $filter_txt;

while(@ARGV >= 1)
{
  if($ARGV[0] =~ /^-?-take_lowest$/i)
  {
    shift;
    $select_dupe = DUPE_SELECT_LOWEST_TAG;
    $select_tag = shift
      or print_usage("Missing INFO tag argument to --take_lowest <tag>");
  }
  elsif($ARGV[0] =~ /^-?-take_highest$/i)
  {
    shift;
    $select_dupe = DUPE_SELECT_HIGHEST_TAG;
    $select_tag = shift
      or print_usage("Missing INFO tag argument to --take_highest <tag>");
  }
  elsif($ARGV[0] =~ /^-?-take_first$/i)
  {
    shift;
    $select_dupe = DUPE_SELECT_FIRST;
  }
  elsif($ARGV[0] =~ /^-?-take_last$/i)
  {
    shift;
    $select_dupe = DUPE_SELECT_LAST;
  }
  elsif($ARGV[0] =~ /^-?-take_none$/i)
  {
    shift;
    $select_dupe = DUPE_SELECT_NONE;
  }
  elsif($ARGV[0] =~ /^-?-take_random$/i)
  {
    print "TAKE RANDOM\n";
    shift;
    $select_dupe = DUPE_SELECT_RANDOM;
  }
  elsif($ARGV[0] =~ /^-?-removed_tag$/i)
  {
    shift;
    $list_removed_tag = shift
      or print_usage("Missing INFO tag argument to --removed_tag <tag>");
  }
  elsif($ARGV[0] =~ /^-?-filter_txt$/i)
  {
    shift;
    $filter_txt = shift or print_usage("Missing FILTER text for --filter_txt");
  }
  elsif(@ARGV > 1)
  {
    print_usage("Unknown option '$ARGV[0]'");
  }
  else
  {
    last;
  }
}

my $vcf_file = shift;

# Set default behaviour
if(!defined($select_dupe))
{
  if(defined($list_removed_tag))
  {
    print_usage("Must use --take_lowest or --take_highest with --removed_tag");
  }

  $select_dupe = defined($filter_txt) ? DUPE_SELECT_NONE : DUPE_SELECT_FIRST;
}

#
# Open VCF Handle
#
my $vcf_handle;

if(defined($vcf_file) && $vcf_file ne "-")
{
  open($vcf_handle, $vcf_file)
    or print_usage("Cannot open VCF file '$vcf_file'");
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
# Read VCF
#
my $vcf = new VCFFile($vcf_handle);

my $num_of_duplicate_sets = 0;
my $num_of_duplicate_vars = 0;

# Print non-PASS variants straight to stdout if -p passed
if(defined($failed_vars_out)) { $vcf->set_filter_failed($failed_vars_out);}

if(defined($filter_txt))
{
  my $description;
  
  if($select_dupe == DUPE_SELECT_FIRST)
  {
    $description = "Duplicated variant that was not the first seen";
  }
  elsif($select_dupe == DUPE_SELECT_LAST)
  {
    $description = "Duplicated variant that was not the last seen";
  }
  elsif($select_dupe == DUPE_SELECT_LOWEST_TAG)
  {
    $description = "Duplicated variant that didn't have the lowest value of " .
                   $select_tag;
  }
  elsif($select_dupe == DUPE_SELECT_HIGHEST_TAG)
  {
    $description = "Duplicated variant that didn't have the highest value of " .
                   $select_tag;
  }
  elsif($select_dupe == DUPE_SELECT_RANDOM)
  {
    $description = "Duplicated variant that wasn't randomly selected";
  }
  else
  {
    $description = "Duplicated variant";
  }

  $vcf->add_header_tag("FILTER", $filter_txt, 0, undef, $description);
}

if(defined($list_removed_tag))
{
  $vcf->add_header_tag("INFO", $list_removed_tag, ".", "Integer",
                       "Comma-separated list of removed values of '$select_tag'");
}

$vcf->print_header();

my $num_of_entries = 0;
my $num_of_printed = 0;

my @variants = ();

my $curr_variant = $vcf->read_entry();

if(defined($curr_variant))
{
  $num_of_entries++;
  push(@variants, $curr_variant);
}

while(defined($curr_variant = $vcf->read_entry()))
{
  $num_of_entries++;

  if($curr_variant->{'CHROM'} ne $variants[0]->{'CHROM'} ||
     $curr_variant->{'POS'} != $variants[0]->{'POS'})
  {
    print_variants();
    @variants = ();
  }

  push(@variants, $curr_variant);
}

print_variants();

print STDERR "vcf_remove_dupes.pl: " . num2str($num_of_duplicate_sets) .
             " sets of duplicates\n";

print STDERR "vcf_remove_dupes.pl: " .
      pretty_fraction($num_of_duplicate_vars, $num_of_entries) . " variants " .
      "involved in duplicates\n";

print STDERR "vcf_remove_dupes.pl: " .
      pretty_fraction($num_of_printed, $num_of_entries) . " variants printed\n";

close($vcf_handle);



sub print_variants
{
  # print "print_variants: " . join(",", map {$_->{'ID'}} @variants) . "\n";

  if(@variants == 1)
  {
    $vcf->print_entry($variants[0]);
    $num_of_printed++;
  }
  else
  {
    # vcf-sort doesn't sort by ALT allele!
    # Have to use VCFFile::vcf_sort_variants

    # Copy and sort variants (by chrom, pos, SVLEN, ref-allele, alt-allele)
    # so as not to upset original ordering
    my @tmp_variants = @variants;
    vcf_sort_variants(\@tmp_variants);

    # Consider lists of dupes together
    my @dupes = ($tmp_variants[0]);

    for(my $i = 1; $i < @tmp_variants; $i++)
    {
      if(uc($tmp_variants[$i]->{'ALT'}) ne uc($dupes[0]->{'ALT'}))
      {
        process_dupes(@dupes);
        @dupes = ();
      }

      push(@dupes, $tmp_variants[$i]);
    }

    process_dupes(@dupes);

    if(defined($filter_txt))
    {
      # If filter_txt is set then we print all
      for my $variant (@variants)
      {
        $vcf->print_entry($variant);
        $num_of_printed++;
      }
    }
    else
    {
      # Print only those that have print defined
      for my $variant (@variants)
      {
        if(defined($variant->{'print'}))
        {
          $vcf->print_entry($variant);
          $num_of_printed++;
        }
      }
    }
  }
}

# Two or more duplicate variants
sub process_dupes
{
  # print "process_dupes: " . join(",", map {$_->{'ID'}} @_) . "\n";

  if(@_ == 1)
  {
    # No dupes - nothing to be done
    $_[0]->{'print'} = 1;
    return;
  }

  $num_of_duplicate_sets++;
  $num_of_duplicate_vars += scalar(@_);

  if($select_dupe == DUPE_SELECT_NONE)
  {
    if(defined($filter_txt))
    {
      # Print all variants, but labelled in the FILTER column
      for my $variant (@_)
      {
        vcf_add_filter_txt($variant, $filter_txt);
      }
    }
  }
  else
  {
    # Select one variant
    my $selected_variant;

    # Check tag is defined for all entries
    if(($select_dupe == DUPE_SELECT_HIGHEST_TAG ||
        $select_dupe == DUPE_SELECT_LOWEST_TAG) &&
       (my @missing = grep {!defined($_->{'INFO'}->{$select_tag})} @_))
    {
      print STDERR "vcf_remove_dupes.pl Error: Variant '$missing[0]->{'ID'}' " .
                   "is missing INFO tag '$select_tag'\n";
      exit;
    }

    if($select_dupe == DUPE_SELECT_HIGHEST_TAG)
    {
      $selected_variant = reduce { $a->{'INFO'}->{$select_tag} >=
                                   $b->{'INFO'}->{$select_tag} ? $a : $b } @_;
    }
    elsif($select_dupe == DUPE_SELECT_LOWEST_TAG)
    {
      $selected_variant = reduce { $a->{'INFO'}->{$select_tag} <=
                                   $b->{'INFO'}->{$select_tag} ? $a : $b } @_;
    }
    elsif($select_dupe == DUPE_SELECT_FIRST)
    {
      $selected_variant = $_[0];
    }
    elsif($select_dupe == DUPE_SELECT_LAST)
    {
      $selected_variant = $_[$#_];
    }
    elsif($select_dupe == DUPE_SELECT_RANDOM)
    {
      $selected_variant = $_[int(rand() * scalar(@_))];
    }

    if(defined($list_removed_tag) && @_ > 1)
    {
      $selected_variant->{'INFO'}->{$list_removed_tag}
        = join(",", map {$_->{'INFO'}->{$select_tag}}
                    grep {$_ != $selected_variant} @_);
    }

    $selected_variant->{'print'} = 1;

    if(defined($filter_txt))
    {
      # Labelling all but the selected one in FILTER column
      for my $variant (@_)
      {
        if($selected_variant != $variant)
        {
          vcf_add_filter_txt($variant, $filter_txt);
        }
      }
    }
  }
}
