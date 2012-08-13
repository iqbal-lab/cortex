#!/usr/bin/perl

use strict;
use warnings;

use Switch;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use VCFFile;
use UsefulModule;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./vcf_filter_by_info.pl [--invert] <file.vcf> [<FIELD> <OP> <SEARCH>] ..
  <OP> must be one of: =~ =~i !~ !~i == != > >= < <= eq eqi ne gt ge lt le
       =~i, !~i & eqi are case-insensitive comparisons
       =~ peforms a regex and !~ is the same as using --invert with =~\n";
  exit;
}

## Test for filtering
my $failed_vars_out = undef;
if(scalar(@ARGV) != scalar(@ARGV = grep {$_ !~ /^-?-p(ass(es)?)?$/i} @ARGV))
{
  open($failed_vars_out, ">-");
}
##

# At least three arguments required for searching
if(@ARGV < 3)
{
  print_usage();
}

my $invert = 0;

if($ARGV[0] =~ /^-?-invert$/i)
{
  $invert = 1;
  shift;
}

my $vcf_file = shift;

if((@ARGV % 3) != 0)
{
  print_usage();
}

my @searches = ();

for(my $i = 0; $i < @ARGV; $i += 3)
{
  my ($field, $op, $search) = ($ARGV[$i], $ARGV[$i+1], $ARGV[$i+2]);

  if($op eq "=")
  {
    $op = "==";
  }

  if(!grep {$op eq $_} qw(=~ =~i !~ !~i == != > >= < <= eq eqi ne gt ge lt le))
  {
    print_usage("Unknown operator '$op'");
  }

  push(@searches, [$field, $op, $search]);
}

#
# Open VCF Handle
#
my $vcf_handle;

if(defined($vcf_file) && $vcf_file ne "-")
{
  open($vcf_handle, $vcf_file) or die("Cannot open VCF file '$vcf_file'\n");
}
elsif(-p STDIN)
{
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

# Print non-PASS variants straight to stdout if -p passed
if(defined($failed_vars_out)) { $vcf->set_filter_failed($failed_vars_out);}

$vcf->print_header();

my $num_of_filtered_entries = 0;
my $total_num_entries = 0;

my $vcf_entry;

while(defined($vcf_entry = $vcf->read_entry()))
{
  $total_num_entries++;

  my $info_hashref = $vcf_entry->{'INFO'};
  my $flags_hashref = $vcf_entry->{'INFO_flags'};

  my $match = 0;

  for my $search_arr (@searches)
  {
    my ($field, $op, $search) = @$search_arr;

    if(defined($flags_hashref->{$field}) &&
       cmp_search($flags_hashref->{$field}, $op, $search) ||
       defined($info_hashref->{$field}) &&
       cmp_search($info_hashref->{$field}, $op, $search))
    {
      $match = 1;
      last;
    }
  }

  if($match != $invert)
  {
    $num_of_filtered_entries++;
    $vcf->print_entry($vcf_entry);
  }
}

# Print stats
print STDERR "vcf_filter_by_info.pl: " .
             pretty_fraction($num_of_filtered_entries, $total_num_entries) .
             " variants printed\n";

close($vcf_handle);


sub cmp_search
{
  my ($value, $op, $search_str) = @_;

  switch($op)
  {
    # Regex
    case "=~" { return ($value =~ /$search_str/); }
    case "=~i" { return ($value =~ /$search_str/i); }
    case "!~" { return ($value !~ /$search_str/); }
    case "!~i" { return ($value !~ /$search_str/i); }

    # Numerical comparisons
    case "==" { return ($value == $search_str); }
    case "!=" { return ($value != $search_str); }
    case ">"  { return ($value >  $search_str); }
    case ">=" { return ($value >= $search_str); }
    case "<"  { return ($value <  $search_str); }
    case "<=" { return ($value <= $search_str); }

    # String comparisons
    case "eq" { return ($value eq $search_str); }
    case "eqi" { return (lc($value) eq lc($search_str)); }
    case "ne" { return ($value ne $search_str); }
    case "gt" { return ($value gt $search_str); }
    case "ge" { return ($value ge $search_str); }
    case "lt" { return ($value lt $search_str); }
    case "le" { return ($value le $search_str); }
  }
}
