#!/usr/bin/perl

use strict;
use warnings;

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use CortexCovgFile;

sub print_usage
{
  for my $err (@_) {
    print "Error: $err\n";
  }

  print STDERR "usage: ./make_branch_file.pl [OPTIONS] [.colour_covgs]\n";
  print STDERR "  Print a FASTA of 5p_flank+branch1+3p_flank from each bubble\n";
  print STDERR "  If [.colour_covgs] omitted (or '-') read from STDIN\n";
  print STDERR "  When called with no OPTIONS branch1 is printed\n";
  print STDERR "\n";
  print STDERR "  OTPIONS:\n";
  print STDERR "    --print5p           Print 5p flank\n";
  print STDERR "    --print3p           Print 3p flank\n";
  print STDERR "    --ref <ref_col>     always prints a branch with covg on ref_col >= 1\n";
  print STDERR "    --branch <1|2|1,2>  print branch 1 or 2 or BOTH\n";
  print STDERR "    --trim <t>          trim sequences longer than <t>\n";
  print STDERR "\n";
  print STDERR "    Note: when both --ref and --branch are specified, will print --branch\n";
  print STDERR "          when no branch has covgs on ref_col >= 1\n";
  
  exit;
}

my $covgs_file;

my $trim;
my $ref_col;
my $branch_to_print;

my $print5p = 0;
my $print3p = 0;

if(@ARGV > 9)
{
  # Too many arguments
  print_usage();
}
elsif(@ARGV > 0 && $ARGV[0] =~ /^--/)
{
  # OPTIONS used

  for(my $i = 0; $i < @ARGV; $i++)
  {
    if($ARGV[$i] =~ /^--?print5p$/i)
    {
      $print5p = 1;
    }
    elsif($ARGV[$i] =~ /^--?print3p$/i)
    {
      $print3p = 1;
    }
    elsif($ARGV[$i] =~ /^-?-trim$/i)
    {
      if($i+1 == @ARGV) {
        print_usage("--trim <t> requires an argument");
      }

      $trim = $ARGV[$i+1];

      if($trim !~ /^\d+$/ || $trim == 0)
      {
        print_usage("Invalid trim value ('$trim') - " .
                    "must be positive integer (>0)");
      }
      
      $i++; # Argument used
    }
    elsif($ARGV[$i] =~ /^-?-ref$/i)
    {
      if($i+1 == @ARGV) {
        print_usage("--ref <col> requires an argument");
      }

      $ref_col = $ARGV[$i+1];
      
      if($ref_col !~ /^\d+$/)
      {
        print_usage("Invalid --ref value ('$ref_col') - " .
                    "must be positive integer (>=0)");
      }
      
      $i++; # Argument used
    }
    elsif($ARGV[$i] =~ /^-?-branch$/i)
    {
      if($i+1 == @ARGV) {
        print_usage("--branches <b> requires an argument");
      }

      $branch_to_print = $ARGV[$i+1];
      
      if($branch_to_print ne "1,2" &&
         $branch_to_print != 1 && $branch_to_print != 2)
      {
        print_usage("Invalid --branch value ('$branch_to_print') - " .
                    "must be either 1 or 2 or '1,2'");
      }
      
      $i++; # Argument used
    }
    elsif($i+1 == @ARGV)
    {
      # Last argument if filename
      $covgs_file = $ARGV[$#ARGV];
    }
    else
    {
      print_usage("Unknown argument '$ARGV[$i]'");
    }
  }
}

#print "trim: $trim; ref: $ref_col; branch_to_print: $branch_to_print;\n";

#
# Open .colour_covgs Handle
#
my $covg_handle;

if(defined($covgs_file) && $covgs_file ne "-") {
  open($covg_handle, $covgs_file)
    or die("Cannot open .colour_covgs file '$covgs_file'\n");
}
elsif(-p STDIN) {
  # STDIN is connected to a pipe
  open($covg_handle, "<&=STDIN") or die("Cannot read pipe");
}
else
{
  print_usage("Must specify or pipe in a .colour_covgs file");
}

my $covgfile = new CortexCovgFile($covg_handle);

# Start reading bubbles
my ($flank_5p, $flank_3p, $branches) = $covgfile->read_bubble_entry();

while(defined($flank_5p))
{
  if(defined($branch_to_print) && $branch_to_print eq "1,2")
  {
    # Just print both branches
    print_seq($branches->[0]->{'seq'});
    print_seq($branches->[1]->{'seq'});
  }
  elsif(defined($ref_col))
  {
    if(min($branches->[0]->{'min_coverage'} >= 1))
    {
      print_seq($branches->[0]->{'seq'});
    }
    elsif(min($branches->[1]->{'min_coverage'} >= 1))
    {
      print_seq($branches->[1]->{'seq'});
    }
    elsif(defined($branch_to_print))
    {
      # default to $branch_to_print
      print_seq($branch_to_print == 1 ? $branches->[0]->{'seq'} : $branches->[1]->{'seq'});
    }
    else
    {
      # die
      print_usage("No branch has ref covg >= 1 ($!)");
    }
  }
  elsif(defined($branch_to_print))
  {
    # print $branch_to_print
    print_seq($branch_to_print == 1 ? $branches->[0]->{'seq'} : $branches->[1]->{'seq'});
  }
  else {
    print_seq($branches->[0]->{'seq'});
  }

  ($flank_5p, $flank_3p, $branches) = $covgfile->read_bubble_entry();
}

close($covg_handle);

sub print_seq
{
  my ($seq) = @_;

  if($print5p)
  {
    $seq = $flank_5p->{'seq'}.$seq;
  }

  if($print3p)
  {
    $seq .= $flank_3p->{'seq'};
  }

  if(defined($trim) && length($seq) > $trim) {
    print substr($seq, 0, $trim);
  }
  else {
    print $seq;
  }

  print "\n";
}
