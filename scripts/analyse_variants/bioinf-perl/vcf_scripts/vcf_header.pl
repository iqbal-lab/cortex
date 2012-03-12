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

  print STDERR "Usage: ./vcf_header.pl [--entries] <vcf> [header1 ..]\n";
  print STDERR "  Print vcf file header (reads from STDIN if <file> is '-').\n";
  print STDERR "  Adds extra headers if passed after file.\n";
  print STDERR "\n";
  print STDERR "  [header] should be of the following forms:\n";
  print STDERR "  1) to add a header:\n";
  print STDERR "    +tag:INFO,ID,Number,Type,Description    e.g. +tag:INFO,T,0,Flag,T number\n";
  print STDERR "    +tag:FORMAT,ID,Number,Type,Description  e.g. +tag:FORMAT,T,1,String,T number\n";
  print STDERR "    +tag:FILTER,ID,Description              e.g. +tag:FILTER,T,T number\n";
  print STDERR "    +tag:ALT,ID,Description                 e.g. +tag:ALT,T,T number\n";
  print STDERR "  2) to remove a header:\n";
  print STDERR "    -tag:ID                                 e.g. -tag:DP\n";
  print STDERR "  3) to add metainfo:\n";
  print STDERR "    +meta:ID=value                          e.g. +meta:date,today\n";
  print STDERR "  4) to remove metainfo:\n";
  print STDERR "    -meta:ID                                e.g. -meta:date\n";
  print STDERR "\n";
  print STDERR "  --entries  Print VCF entries as well as header\n";
  exit;
}

if(@ARGV == 0)
{
  print_usage();
}

my $print_entries = 0;

if(@ARGV > 1 && $ARGV[0] =~ /^-?-e(ntries)?$/i)
{
  $print_entries = 1;
  shift;
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

$vcf->print_header();

if($print_entries)
{
  my $vcf_entry;

  while(defined($vcf_entry = $vcf->read_entry()))
  {
    $vcf->print_entry($vcf_entry);
  }
}

close($vcf_handle);
