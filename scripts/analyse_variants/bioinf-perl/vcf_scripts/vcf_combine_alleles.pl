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
"Usage: ./vcf_combine_alleles.pl [in.vcf]
  Combines variants (and their alleles) that have matching ref sites.
  Prints to stdout.  If no in.vcf given, or '-', reads from STDIN\n";

  exit;
}

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

$vcf->print_header();

my @samples = $vcf->get_list_of_sample_names();

my $vcf_entry;

my @vars_at_same_pos = ();

# Stats
my $num_of_entries = 0;
my $num_of_printed = 0;
my $num_of_merges = 0;
my $num_of_vars_merged = 0;
my $biggest_merge = 0;

# Read first variant
if(defined($vcf_entry = $vcf->read_entry()))
{
  $num_of_entries++;
  @vars_at_same_pos = ($vcf_entry);
}

while(defined($vcf_entry = $vcf->read_entry()))
{
  $num_of_entries++;

  if($vcf_entry->{'CHROM'} ne $vars_at_same_pos[0]->{'CHROM'} ||
     $vcf_entry->{'POS'} != $vars_at_same_pos[0]->{'POS'})
  {
    # Parse current list
    print_list();

    @vars_at_same_pos = ();
  }

  push(@vars_at_same_pos, $vcf_entry);
}

# Parse remaining list
print_list();

# Print stats

my $mean_vars_per_merge
  = ($num_of_merges == 0 ? 0 : ($num_of_vars_merged / $num_of_merges));

print STDERR "vcf_combine_alleles.pl: mean number of variants per merge is " .
             num2str($mean_vars_per_merge) . "\n";

print STDERR "vcf_combine_alleles.pl: biggest merge has " .
             num2str($biggest_merge) . " variants\n";

print STDERR "vcf_combine_alleles.pl: " .
             pretty_fraction($num_of_printed, $num_of_entries) .
             " variants printed\n";

print STDERR "vcf_combine_alleles.pl: " .
             pretty_fraction($num_of_merges, $num_of_printed) .
             " printed variants are merges\n";

close($vcf_handle);





# Merge and print list of variants starting at the same position
sub print_list
{
  if(@vars_at_same_pos == 1)
  {
    $vcf->print_entry($vars_at_same_pos[0]);
    $num_of_printed++;
  }
  else
  {
    # Sort byt length of ref allele
    @vars_at_same_pos = sort {length($a->{'REF'}) <=> length($b->{'REF'})}
                             @vars_at_same_pos;

    for(my $i = 0; $i < @vars_at_same_pos; $i++)
    {
      if($i+1 == @vars_at_same_pos ||
         length($vars_at_same_pos[$i]->{'REF'}) !=
           length($vars_at_same_pos[$i+1]->{'REF'}))
      {
        $vcf->print_entry($vars_at_same_pos[$i]);
        $num_of_printed++;
      }
      else
      {
        # Two variants [$i,$i+1] share the same ref allele - look for more
        my $start = $i;
        my $end = $i+1;

        $num_of_merges++;

        while($end+1 < @vars_at_same_pos &&
              length($vars_at_same_pos[$end]->{'REF'}) ==
                length($vars_at_same_pos[$end+1]->{'REF'}))
        {
          $end++;
        }

        # @vars_at_same_pos[$start..$end] all share the same ref site
        # merge and print a single entry
        my $alts = join(",", map {$_->{'ALT'}} @vars_at_same_pos[$start..$end]);

        # Use $vars_at_same_pos[$start] for merge
        $vars_at_same_pos[$start]->{'ALT'} = $alts;

        # Reset samples
        for my $sample (@samples)
        {
          $vars_at_same_pos[$start]->{$sample} = '.';
        }

        # Set filter field
        vcf_add_filter_txt($vars_at_same_pos[$start], 'MULTIALLELIC');

        # Print single entry
        $vcf->print_entry($vars_at_same_pos[$start]);
        $num_of_printed++;

        my $vars_merged = $end-$start+1;

        $biggest_merge = max($biggest_merge, $vars_merged);
        $num_of_vars_merged += $vars_merged;

        # Update $i
        $i = $end;
      }
    }
  }
}

