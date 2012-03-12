package CortexCovgFile;

use strict;
use warnings;

use Carp;

use List::Util qw(min max sum);


# perl library code for handling data from the cortex variant caller
# url: http://github.com/noporpoise/bioinf-perl
# author: Isaac Turner <turner.isaac@gmail.com>

# Cortex can be found here:
# http://cortexassembler.sourceforge.net/
# (by Zam Iqbal & Mario Caccamo)

# Please reference:
# "De novo assembly and genotyping of variants using colored de Bruijn graphs",
# Iqbal(*), Caccamo(*), Turner, Flicek, McVean (Nature Genetics) (2012)
# (doi:10.1038/ng.1028)


# Object methods
# grab_bubble_entry read_bubble_entry read_align_entry

use base 'Exporter';

our @EXPORT = qw{get_num_of_read_arrivals estimate_contig_rate};


sub new
{
  my $class = shift;
  my $handle = shift;

  my $next_line = <$handle>;

  my $self = {
      _handle => $handle,
      _next_line => $next_line
  };

  bless $self, $class;
  return $self;
}

sub peak_line
{
  my ($self) = @_;
  return $self->{_next_line};
}

sub read_line
{
  my ($self) = @_;
  my $temp_line = $self->{_next_line};
  my $handle = $self->{_handle};
  
  $self->{_next_line} = <$handle>;
  
  return $temp_line;
}

# Don't parse, just read
sub grab_bubble_entry
{
  my ($self) = @_;

  if(!defined($self->peak_line()))
  {
    return undef;
  }

  my $entry = $self->read_line(); # 5p_flank

  if(!defined($entry))
  {
    return undef;
  }

  my $peak;
  
  while(defined($peak = $self->peak_line()) && $peak !~ /^(?:Colour|>var_\d+_5p_flank)/i)
  {
    $entry .= $self->read_line();
  }

  return $entry;
}

sub read_bubble_entry
{
  my ($self) = @_;

  # Check there is another entry
  if(!defined($self->peak_line()))
  {
    return undef;
  }
  
  # Read likelihoods if they are there
  my $col_llk;

  if($self->peak_line() =~ /^Colour\/sample/i)
  {
    $col_llk = {};

    # e.g.:
    # Colour/sample	GT_call	llk_hom_br1	llk_hom_br2
    # 0	HOM2	-4.68	-2.45
    
    my $likelihood_header = $self->read_line();
    
    while(($self->peak_line() =~ /^\d+/))
    {
      my $likelihood_line = $self->read_line();
      my ($col, $gt_call, $llk_hom_br1, $llk_hom_br2)
        = split(/\s/, $likelihood_line);

      $col_llk->{$col} = [$gt_call, $llk_hom_br1, $llk_hom_br2];
    }
  }

  my $flank_5p = $self->parse_bubble_graphline();
  my $branch1 = $self->parse_bubble_graphline();
  my $branch2 = $self->parse_bubble_graphline();
  my $flank_3p = $self->parse_bubble_graphline();
  
  # Read empty line
  $self->read_line();
  
  # Read coverages for branch1
  $self->parse_bubble_covgs('1',$branch1);
  
  # Read coverages for branch2
  $self->parse_bubble_covgs('2',$branch2);
  
  # Read two empty lines
  $self->read_line();
  $self->read_line();
  
  # May add more branches later if needed
  my $branches = [$branch1, $branch2];

  return ($flank_5p, $flank_3p, $branches, $col_llk);
}

sub parse_bubble_covgs
{
  my ($self,$branch_num,$branch_hashref) = @_;

  # Array for colour coverages on this branch
  $branch_hashref->{'covgs'} = [];

  my $covg_line = $self->read_line();
  chomp($covg_line);

  if($covg_line ne "branch$branch_num coverages")
  {
    croak("Expected branch1 covgs line ('$covg_line')");
  }

  while($self->peak_line() =~ /^Covg in Colour (\d+):/)
  {
    my $colour = $1;
    $self->read_line();
    my $covgs_line;
    
    if(!defined($covgs_line = $self->read_line()))
    {
      croak("Missing branch$branch_num covg line (premature file end)");
    }
    elsif($covgs_line !~ /^\d+(?:\s+\d+)*\s*$/)
    {
      chomp($covgs_line);
      croak("Invalid branch$branch_num covg line ('$covgs_line')");
    }
    
    my @covgs = split(/\s+/, $covgs_line);
    $branch_hashref->{'covgs'}->[$colour] = \@covgs;
  }
}

sub parse_bubble_graphline
{
  my ($self) = @_;

  my $header_line;
  my $seq_line;
  
  if(!defined($header_line = $self->read_line()) ||
     !defined($seq_line = $self->read_line()))
  {
    croak("Missing cortex bubble lines");
  }

  chomp($header_line);
  chomp($seq_line);

  my %data = ('seq' => $seq_line);

  my $regex = '^>(\S+)\s+length:\s*(\d+)\s*average_coverage:\s*([\d\.]+)\s*' .
              'min_coverage:\s*(\d+)\s*max_coverage:\s*(\d+)\s*' .
              'fst_coverage:\s*(\d+)\s*fst_kmer:\s*([ACGT]*)\s*' .#fst
              'fst_r:\s*([ACGT]*)\s*fst_f:([ACGT]*)\s*' .
              'lst_coverage:\s*(\d+)\s*lst_kmer:\s*([ACGT]*)\s*' .#lst
              'lst_r:\s*([ACGT]*)\s*lst_f:([ACGT]*)';

  my @keys = qw{name length average_coverage min_coverage max_coverage
                fst_coverage fst_kmer fst_r fst_f
                lst_coverage lst_kmer lst_r fst_f};

  if($header_line =~ /$regex/)
  {
    @data{@keys} = ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13);
  }
  else
  {
    croak("Unexpected bubble line '$header_line'");
  }

  return \%data;
}

#
# Read align .colour_covgs output
#
sub read_align_entry
{
  my ($self) = @_;

  my $read_name = $self->read_line();
  
  if(!defined($read_name)) {
    return undef;
  }
  
  chomp($read_name);
  
  if($read_name =~ /^>(.*)/)
  {
    $read_name = $1;
  }
  else {
    croak("Unexpected .colour_covgs read name '$read_name': $!");
  }
  
  my $sequence = $self->read_line();
  
  if(!defined($sequence))
  {
    croak("Unexpected end of .colour_covgs file (read '$read_name'): $!");
  }
  
  chomp($sequence);
  
  my %colour_lines = ();
  
  my $expected_line = '>'.quotemeta($read_name).'_colour_(\d+)_kmer_coverages';
  my $peak_line;
  
  while(defined($peak_line = $self->peak_line()) &&
        $peak_line =~ /^$expected_line/)
  {
    my $cortex_colour = $1;

    # Read line and ignore
    $self->read_line();

    # read coverages for colour
    my $covg_line = $self->read_line();

    if(!defined($covg_line)) {
      croak("Unexpected end of .colour_covgs file (read '$read_name'): $!");
    }
    
    chomp($covg_line);
    my @kmer_covgs = split(/\s/, $covg_line);
    $colour_lines{$cortex_colour} = \@kmer_covgs;
  }
  
  return ($read_name, $sequence, \%colour_lines);
}


#
# static methods
#

sub estimate_contig_rate
{
  my ($covgs_ptr, $kmer_size, $read_length, $epsilon,
      $autosome_depth, $x_depth, $y_depth) = @_;
  
  if(!defined($y_depth)) {
    croak("estimate_contig_rate(..) Missing arguments $!");
  }
  
  my ($num_of_read_arrivals) = get_num_of_read_arrivals($covgs_ptr);

  my $contig_length_in_k = @$covgs_ptr;
  my $epsilon_correction = (1 - $kmer_size * $epsilon);
  
  my $auto_count = ($num_of_read_arrivals * $read_length) /
                   ($epsilon_correction * $contig_length_in_k * $autosome_depth);

  my $x_count = ($num_of_read_arrivals * $read_length) /
                ($epsilon_correction * $contig_length_in_k * $x_depth);

  my $y_count = ($num_of_read_arrivals * $read_length) /
                ($epsilon_correction * $contig_length_in_k * $y_depth);

  return ($auto_count, $x_count, $y_count);
}

sub get_num_of_read_arrivals
{
  my ($covgs_ptr) = @_;

  if(!defined($covgs_ptr)) {
    croak("get_num_of_read_arrivals_zam(..) Missing arguments $!");
  }

  my $total = $covgs_ptr->[1];

  my $num_of_kmers = @$covgs_ptr;

  #ignore first and last
  for(my $i = 2; $i < $num_of_kmers-1; $i++)
  {
    my $jump = $covgs_ptr->[$i] - $covgs_ptr->[$i-1];
    my $next_jump = $covgs_ptr->[$i+1] - $covgs_ptr->[$i];
    
    if($jump > 0 && -$next_jump != $jump)
    {
      # ie provided the jump was not a spike on a single kmer
      $total += $jump;
    }
  }
  
  return $total;
}

1;
