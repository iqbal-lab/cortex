#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use VCFFile;

## Config
my $csvsep = ",";
my $info_prefix = "INFO_";
##

sub print_usage
{
  if(@_) {
    print STDERR "Error: " . join('; ', @_) . "\n";
  }
  
  print STDERR "Usage: ./vcf_to_csv.pl [vcf]\n";
  print STDERR "  Converts VCF to CSV\n";
  exit;
}

## Test for filtering
my $skip_failed_vars = 0;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  $skip_failed_vars = 1;
}
##

if(@ARGV > 1)
{
  print_usage();
}

my $vcf_file = shift;

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
# Open VCF
#
my $vcf = new VCFFile($vcf_handle);

# Skip non-PASS variants if -p passed
if($skip_failed_vars) { $vcf->set_filter_failed(undef);}

# Miss out INFO and FORMAT for now
my @vcf_standard_cols = qw(CHROM POS ID REF ALT QUAL FILTER);
my %vcf_file_cols = $vcf->get_columns_hash();

my @vcf_cols = grep {defined($vcf_file_cols{$_})} @vcf_standard_cols;

# Add some extra columns
push(@vcf_cols, ('true_REF','true_ALT','true_POS'));

# Get info fields from the header
my %header_tags = $vcf->get_header_tags();

my @info_tags_from_header =
  map {$_->{'ID'}} grep {$_->{'column'} eq "INFO"} values %header_tags;

my %info_fields_hash = ();
@info_fields_hash{@info_tags_from_header} = 1;

# Read the first VCF entry
my $vcf_entry = $vcf->read_entry();

my @additional_info_fields = (keys %{$vcf_entry->{'INFO'}},
                              keys %{$vcf_entry->{'INFO_flags'}});

for my $additional_info_field (@additional_info_fields)
{
  $info_fields_hash{$additional_info_field} = 1;
}

my @info_fields = sort {$a cmp $b} keys %info_fields_hash;

# Get sample names
my @sample_names = $vcf->get_list_of_sample_names();

# Print header
my @all_cols = (@vcf_cols, map {$info_prefix.$_} @info_fields);

if(defined($vcf_file_cols{'FORMAT'}))
{
  push(@all_cols, 'FORMAT');
}

push(@all_cols, @sample_names);

print join($csvsep, @all_cols) . "\n";

#
# Read VCF
#
while(defined($vcf_entry))
{
  # Do simple columns
  print join($csvsep, map {$vcf_entry->{$_}} @vcf_cols);

#  print $vcf_entry->{$vcf_cols[0]};

#  for(my $i = 1; $i < @vcf_cols; $i++) {
#    print $csvsep . $vcf_entry->{$vcf_cols[$i]};
#  }

  # Print info columns
  for my $info_field (@info_fields)
  {
    print $csvsep;

    my $value;
    
    if(defined($value = $vcf_entry->{'INFO'}->{$info_field}))
    {
      print $value;
    }
    elsif(defined($vcf_entry->{'INFO_flags'}->{$info_field}))
    {
      # Print name of flag in info flag field
      print $info_field;
    }
  }
  
  # Print FORMAT and sample columns
  if(defined($vcf_entry->{'FORMAT'}))
  {
    my $format_arr = $vcf_entry->{'FORMAT'};
    print $csvsep . join(":", @$format_arr);

    # Print sample columns
    for my $sample (@sample_names)
    {
      print $csvsep . join(":", map {$vcf_entry->{$sample}->{$_}} @$format_arr);
    }
  }

  print "\n";

  # Read next entry
  $vcf_entry = $vcf->read_entry();
}

close($vcf_handle);
