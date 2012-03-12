#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use VCFFile;

## Config
my $csvsep = "\t";
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

if(@ARGV > 1) {
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
# Read VCF
#
my $vcf = new VCFFile($vcf_handle);

my @vcf_cols = $vcf->get_columns_array();

@vcf_cols = grep {$_ !~ /^INFO$/i} @vcf_cols;
push(@vcf_cols, ('true_REF','true_ALT','true_POS'));

# Get info fields
my %info_fields_hash = ();

my %header_tags = $vcf->get_header_tags();

while(my ($id, $tag) = each(%header_tags))
{
  if($tag->{'column'} eq "INFO")
  {
    $info_fields_hash{$id} = 1;
  }
}

# Read VCF entry
my $vcf_entry = $vcf->read_entry();

my @additional_info_fields = (keys %{$vcf_entry->{'INFO'}},
                              keys %{$vcf_entry->{'INFO_flags'}});

for my $additional_info_field (@additional_info_fields)
{
  $info_fields_hash{$additional_info_field} = 1;
}

my @info_fields = sort {$a cmp $b} keys %info_fields_hash;

print join($csvsep, (@vcf_cols, @info_fields)) . "\n";

while(defined($vcf_entry))
{
  print join($csvsep, map {$vcf_entry->{$_}} @vcf_cols);

#  print $vcf_entry->{$vcf_cols[0]};

#  for(my $i = 1; $i < @vcf_cols; $i++) {
#    print $csvsep . $vcf_entry->{$vcf_cols[$i]};
#  }

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
  
  print "\n";
  
  $vcf_entry = $vcf->read_entry();
}

close($vcf_handle);
