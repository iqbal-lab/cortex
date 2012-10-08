#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use VCFFile;

sub print_usage
{
  for my $err (@_) {
    print "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./vcf_header.pl [--entries|--noheader] <vcf> [header1 ..]
  Print vcf file header (reads from STDIN if <file> is '-').
  Adds extra headers if passed after file.

  [header] should be of the following forms:
  1) to add a header:
    +tag:INFO,ID,Number,Type,Description    e.g. +tag:INFO,T,0,Flag,T number
    +tag:FORMAT,ID,Number,Type,Description  e.g. +tag:FORMAT,T,1,String,T number
    +tag:FILTER,ID,Description              e.g. +tag:FILTER,T,T number
    +tag:ALT,ID,Description                 e.g. +tag:ALT,T,T number
  2) to remove a header:
    -tag:ID                                 e.g. -tag:DP
  3) to add metainfo:
    +meta:ID=value                          e.g. +meta:date,today
  4) to remove metainfo:
    -meta:ID                                e.g. -meta:date

  --entries  Print VCF entries as well as header
  --noheader Don't print the header\n";

  exit;
}

## Test for filtering
my $failed_vars_out = undef;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  open($failed_vars_out, ">-");
}
##

if(@ARGV == 0)
{
  print_usage();
}

my $print_entries = 0;
my $no_header = 0;

for my $arg (@ARGV)
{
  if($arg =~ /^-?-e(ntries)?$/i)
  {
    $print_entries = 1;
    shift;
  }
  elsif($arg =~ /^-?-n(oheader)?$/i)
  {
    $no_header = 1;
    shift;
  }
  else
  {
    last;
  }
}

my $vcf_file = shift;

my @extra_headers = @ARGV;

#
# Open VCF Handle
#
my $vcf_handle;

if(defined($vcf_file) && $vcf_file ne "-") {
  open($vcf_handle, $vcf_file) or print_usage("Cannot open VCF file '$vcf_file'\n");
}
elsif(-p STDIN) {
  # STDIN is connected to a pipe
  open($vcf_handle, "<&=STDIN") or print_usage("Cannot read pipe");
}
else
{
  print_usage("Must specify or pipe in a VCF file");
}

my $vcf = new VCFFile($vcf_handle);

# Print non-PASS variants straight to stdout if -p passed
if(defined($failed_vars_out)) { $vcf->set_filter_failed($failed_vars_out);}

for my $extra_header (@extra_headers)
{
  my $add_header;
  my $header_metainfo;
  my $header_txt;

  if($extra_header =~ /^([\-\+]?)(meta|tag):(.*)$/i)
  {
    $add_header = ($1 eq "+" || $1 eq "");
    $header_metainfo = (lc($2) eq "meta");
    $header_txt = $3;
  }
  else
  {
    print_usage("Invalid header argument '$extra_header'");
  }

  if($header_metainfo)
  {
    if($add_header)
    {
      my @parts = split("=", $header_txt);
      
      if(@parts != 2)
      {
        print_usage("Invalid add metainfo tag '$extra_header'");
      }

      $vcf->add_header_metainfo($parts[0], $parts[1]);
    }
    else
    {
      if(!defined(get_header_metainfo($header_txt)))
      {
        warn("Metainfo tag '$header_txt' in VCF file - cannot remove\n");
      }
    
      $vcf->remove_metainfo($header_txt);
    }
  }
  else
  {
    # Header tag
    ##INFO=<ID=id,...>
    if($add_header)
    {
      # Add header tag
      my @parts = split(",", $header_txt);
      my ($tag_column, $tag_id, $tag_number, $tag_type, $tag_description);

      if(@parts == 3)
      {
        ($tag_column,$tag_id,$tag_description) = @parts;
      }
      elsif(@parts == 5)
      {
        ($tag_column, $tag_id, $tag_number, $tag_type, $tag_description) = @parts;
      }
      else
      {
        print_usage("Invalid extra header argument: '$extra_header'");
      }
    
      if($tag_description =~ /^\".*\"$/)
      {
        # Trim first and last characters off
        $tag_description = substr($tag_description, 1, -1);
      }
      elsif($tag_description =~ /\"/)
      {
        print_usage("Description should not contain quotation marks");
      }

      $vcf->add_header_tag($tag_column, $tag_id, $tag_number,
                           $tag_type, $tag_description);
    }
    else
    {
      # Remove header tag
      if(!defined(get_header_tag($header_txt)))
      {
        warn("Header tag '$header_txt' in VCF file - cannot remove\n");
      }

      $vcf->remove_header_tag($header_txt);
    }
  }
}

if(!$no_header)
{
  $vcf->print_header();
}

if($print_entries)
{
  my $vcf_entry;

  while(defined($vcf_entry = $vcf->read_entry()))
  {
    $vcf->print_entry($vcf_entry);
  }
}

close($vcf_handle);
