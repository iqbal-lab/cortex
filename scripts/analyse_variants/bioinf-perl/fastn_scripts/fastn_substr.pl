#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use FASTNFile;
use RefGenome;
use UsefulModule;

sub print_usage
{
  for my $err (@_) {
    print STDERR "Error: $err\n";
  }

  print STDERR "" .
"Usage: ./fastn_substr.pl [--allow_dupes|--file <f.pos>] [chr:pos:len,chr:pos-end,.] " .
  "[in1.fa in2.fq ..]\n" .
"  Prints substrings from FASTA/Q files (or STDIN if '-').  Takes comma-separated\n" .
"  list of regions.\n" .
"  \n" .
"  * Regions must be 'chr:start:length,..' or 'chr:start-end,..' or 'chrX:*'\n" .
"  * (coordinates are 1-based)  \n" .
"  * Start/end positions <0 mean X bp from the end\n" .
"  * --file <f.pos>  is a file with each line chr:pos:len or chr:pos-end\n" .
"  * --allow_dupes   if multiple entries have the same name\n";

  exit;
}

if(@ARGV < 1)
{
  print_usage();
}

my @search_files = ();
my $searches_list;

my $allow_dupes = 0;

while(@ARGV >= 1)
{
  if(@ARGV >= 2 && $ARGV[0] =~ /^-?-f(ile)?$/i)
  {
    shift;
    push(@search_files, shift);
  }
  elsif($ARGV[0] =~ /^-?-allow_dupes?$/i)
  {
    $allow_dupes = 1;
  }
  else
  {
    last;
  }
}

if(@ARGV >= 0 && $ARGV[0] =~ /.+:(\d+|\*)/)
{
  $searches_list = shift;
}
elsif(@search_files == 0)
{
  print_usage();
}

my @ref_files = @ARGV;

# Check we can read ref files

for my $ref_file (@ref_files)
{
  if(!(-r $ref_file))
  {
    print_usage("Cannot read fasta/q file '$ref_file'");
  }
}

my %search_chrs = ();

# Load searches from files
for my $search_file (@search_files)
{
  open(SEARCH, $search_file)
    or print_usage("Couldn't open file of coords '$search_file'");

  my $search_line;
  my $search_line_number = 1;

  while(defined($search_line = <SEARCH>))
  {
    chomp($search_line);
    my @searches = split(",", $search_line);
  
    for my $search (@searches)
    {
      parse_search_term($search,
                        " [file:".$search_file.":".$search_line_number."]");
    }

    $search_line_number++;
  }

  close(SEARCH);
}

# Load searches from cmd line
if(defined($searches_list))
{
  my @searches = split(",", $searches_list);

  for my $search (@searches)
  {
    parse_search_term($search);
  }
}

#
# Open FASTA/Q handles
#
if(@ref_files == 0)
{
  # Add stdin
  push(@ref_files, "-");
}

for my $file (@ref_files)
{
  my $fastn_handle;

  if($file eq "-")
  {
    $fastn_handle = open_stdin("Cannot read ref -- need to pipe in fasta/fastq");
  }
  else
  {
    open($fastn_handle, $file)
      or print_usage("Cannot open FASTA/Q file '$file'");
  }
  
  search_file($fastn_handle);

  close($fastn_handle);

  if(keys(%search_chrs) == 0)
  {
    # All entries have been found - no need to open any more files
    last;
  }
}

# Done

if(keys(%search_chrs) > 0)
{
  # Print warnings for search results with no match
  for my $chr (sort {$a cmp $b} keys %search_chrs)
  {
    my $array_ref = $search_chrs{$chr};

    for(my $i = 0; $i < @$array_ref; $i++)
    {
      my $term_str = search_term_to_str($chr, $array_ref->[$i]);
      print STDERR "# fastn_substr.pl Warning: Couldn't find $term_str\n";
    }
  }
}


#
# Functions
#

sub search_file
{
  my ($handle) = @_;

  my $fastn = new FASTNFile($handle);

  # Read in
  my $array_ref;

  while(my ($fasta_name, $seq) = $fastn->read_next())
  {
    my $search_name = $fasta_name;
    my $array_ref = $search_chrs{$search_name};

    if(!defined($array_ref))
    {
      $search_name = guess_plain_name($fasta_name, keys %search_chrs);

      if(defined($search_name))
      {
        $array_ref = $search_chrs{$search_name};
      }
    }

    if(defined($array_ref))
    {
      # Print substrings on this chromosome
      for(my $i = 0; $i < @$array_ref; $i++)
      {
        my ($start, $length, $end, $all) = @{$array_ref->[$i]};

        if($all)
        {
          print ">$search_name\n";
          print "$seq\n";
        }
        else
        {
          my $search_str = search_term_to_str($search_name, $array_ref->[$i]);

          if($start < 0)
          {
            # still in 1-based coords
            $start += length($seq)+1;
          }

          if(!defined($length))
          {
            if($end < 0)
            {
              $end += length($seq)+1;
            }

            $length = $end-$start+1;
            print ">$search_str\n";
          }
          else
          {
            $end = $start+$length-1;
            print ">$search_str\n";
          }

          if($end > length($seq))
          {
            
            print STDERR "# Warning: $search_str is out of bounds of " .
                         "$fasta_name:1:".length($seq)."\n";
          }

          # Correct for 1-based coords here (convert to 0-based)
          print substr($seq, $start-1, $length)."\n";
        }
      }

      if(!$allow_dupes)
      {
        # Only return the first result for each search - so exit
        delete($search_chrs{$search_name});
      }
    }

    if(keys(%search_chrs) == 0)
    {
      last;
    }
  }

  close($handle);
}

sub search_term_to_str
{
  my $chrom_str = shift;
  my @results = ();

  foreach my $term (@_)
  {
    my ($start, $length, $end, $all) = @$term;

    if(defined($all))
    {
      push(@results, $chrom_str.":*");
    }
    elsif(defined($length))
    {
      push(@results, $chrom_str.":".$start.":".$length);
    }
    else
    {
      push(@results, $chrom_str.":".$start."-".$end);
    }
  }

  return @results == 1 && !wantarray() ? $results[0] : @results;
}

sub parse_search_term
{
  my ($search, $file_err) = @_;

  if(!defined($file_err))
  {
    $file_err = "";
  }

  # All initially undefined - some will be defined
  my ($chr, $start, $length, $end, $all);

  if($search =~ /^(.*):(-?\d+):(\d+)$/)
  {
    $chr = $1;
    $start = $2;
    $length = $3;
    
    if($start == 0) {
      print_usage("Start position is 1-based - cannot be 0");
    }
    elsif($length == 0) {
      print_usage("Substring length cannot be 0");
    }
    elsif($start < 0 && $length > -$start)
    {
      print_usage("Start position is " . (-$start) . "bp from the end of the " .
                  "read, but length is ".$length."bp (length too long!)" .
                  $file_err);
    }
  }
  elsif($search =~ /^(.*):(-?\d+)-(-?\d+)$/)
  {
    $chr = $1;
    $start = $2;
    $end = $3;
    
    if($start == 0) {
      print_usage("Start position is 1-based - cannot be 0");
    }
    elsif($end == 0) {
      print_usage("End position cannot be 0");
    }
    elsif($start < 0 && $end < $start)
    {
      print_usage("Start position is " . (-$start) . "bp from the end of the " .
                  "read, but end is ".(-$end)."bp from the end (negative length!)" .
                  $file_err);
    }
  }
  elsif($search =~ /^(.*):\*$/)
  {
    $chr = $1;
    $all = 1;
  }
  else
  {
    print_usage("Invalid position argument '$search'" .
                $file_err);
  }

  if(!defined($search_chrs{$chr}))
  {
    # First search position on this chromosome
    $search_chrs{$chr} = [];
  }
  
  push(@{$search_chrs{$chr}}, [$start, $length, $end, $all]);
}
