package RefGenome;

use strict;
use warnings;

use Carp;

use List::Util qw(min max sum);

# Use current directory to find modules
use FindBin;
use lib $FindBin::Bin;

use FASTNFile;

# All methods are object methods except these:
use base 'Exporter';
our @EXPORT = qw(guess_plain_name guess_fasta_name guess_name_pairs);

sub new
{
  my $class = shift;

  my %options = (
      uppercase => 0,
      lowercase => 0,
      first_base_index => 0,
      @_
    );

  if($options{'uppercase'} && $options{'lowercase'})
  {
    croak("Can't define upper and lowercase in new RefGenome options");
  }

  my $self = {
    _chroms => {},
    _quals => {},
    _lengths => {},
    _guess_names => {},
    _uppercase => $options{'uppercase'},
    _lowercase => $options{'lowercase'},
    _first_base_index => $options{'first_base_index'}
  };

  bless $self, $class;
  return $self;
}

sub remove_chromosomes
{
  my $self = shift;

  for my $chrom (@_)
  {
    $self->{_chroms}->{$chrom} = undef;
    $self->{_quals}->{$chrom} = undef;
    $self->{_lengths}->{$chrom} = undef;
  }
}

sub add_chr
{
  my ($self, $name, $seq, $quals) = @_;
  
  if($self->{'uppercase'})
  {
    $seq = uc($seq);
  }
  elsif($self->{'lowercase'})
  {
    $seq = lc($seq);
  }

  $self->{_chroms}->{$name} = $seq;
  $self->{_lengths}->{$name} = length($seq);
  $self->{_quals}->{$name} = $seq;
}

sub remove_chr
{
  my ($self, $name) = @_;

  my $chrom = guess_chrom_fasta_name($name);

  if(defined($chrom))
  {
    delete($self->{_chroms}->{$chrom});
    delete($self->{_lengths}->{$chrom});
    delete($self->{_quals}->{$chrom});

    return 1;
  }
  else
  {
    return 0;
  }
}

sub load_from_files
{
  my ($self, @files) = @_;

  my $uppercase = $self->{'uppercase'};
  my $lowercase = $self->{'lowercase'};

  my ($ref_genomes_hashref, $quals) = read_all_from_files(@files);

  while(my ($name, $seq) = each(%$ref_genomes_hashref))
  {
    if($uppercase)
    {
      $seq = uc($seq);
    }
    elsif($lowercase)
    {
      $seq = lc($seq);
    }

    $self->{_chroms}->{$name} = $seq;
    $self->{_lengths}->{$name} = length($seq);
    $self->{_quals}->{$name} = $quals->{$name};
  }

  # Reset guessed names
  $self->{_guess_names} = {};
}

# For a single FASTA chrom name and a list of plain chrom names,
# return the plain chrom name that matches
sub guess_plain_name
{
  my $fasta_name = shift;

  my ($plain_names, $fasta_names) = guess_name_pairs(\@_, [$fasta_name]);

  return $plain_names->[0];
}

# For a single VCF chrom name and a list of FASTA chrom names,
# return the FASTA name that matches
sub guess_fasta_name
{
  my $plain_name = shift;

  my ($plain_names, $fasta_names) = guess_name_pairs([$plain_name], \@_);

  return $fasta_names->[0];
}

sub guess_name_pairs
{
  my ($plain_names_arr, $fasta_names_arr) = @_;

  my @results_plain = ();
  my @results_fasta = ();

  for my $plain_name_full (@$plain_names_arr)
  {
    my $plain_name = $plain_name_full;
    $plain_name =~ s/^chr//g;

    my @regexes = ('^\s*(chr(om(osome)?)?)?\s*'.$plain_name.'\b',
                   '\b(chr(om(osome)?)?)?\s*'.$plain_name.'\b',
                   '\b(chr(om(osome)?)?):\w+:(chr(om(osome)?)?)?'.$plain_name.'\b',
                   '([_\-:]|chr|chrom|chromosome|\b)'.$plain_name.'\b');

    my $result_plain = $plain_name;
    my $result_fasta;

    for my $regex (@regexes)
    {
      for my $fasta_name (@$fasta_names_arr)
      {
        if($fasta_name =~ /$regex/i)
        {
          # Look for just the name at the start, e.g. '>10' '> chr10' '> chromosome 10'
          $result_plain = $plain_name_full;
          $result_fasta = $fasta_name;
          last;
        }
      }

      if(defined($result_fasta))
      {
        last;
      }
    }

    push(@results_plain, $result_plain);
    push(@results_fasta, $result_fasta);
  }

  return (\@results_plain, \@results_fasta);
}

sub guess_chrom_fasta_name
{
  my ($self, $plain_name) = @_;

  if(defined($self->{_chroms}->{$plain_name}))
  {
    return $plain_name;
  }

  my @ref_chroms = keys(%{$self->{_chroms}});

  if(@ref_chroms == 0)
  {
    return undef;
  }

  if(defined(my $guess = $self->{_guess_names}->{$plain_name}))
  {
    return $guess;
  }

  my $fasta_name = guess_fasta_name($plain_name, @ref_chroms);

  if(defined($fasta_name))
  {
    $self->{_guess_names}->{$plain_name} = $fasta_name;
    return $fasta_name;
  }
  else
  {
    return undef;
  }
}

sub get_guesses_hashref
{
  my ($self) = @_;

  return $self->{_guess_names};
}

sub get_chr_names
{
  my ($self) = @_;

  return keys %{$self->{_chroms}};
}

sub get_chr_lengths
{
  my ($self) = @_;

  return values %{$self->{_lengths}};
}

sub chr_exists
{
  my ($self, $chrom) = @_;

  my $name = $self->guess_chrom_fasta_name($chrom);

  return defined($name);
}

sub get_chr
{
  my ($self, $chrom) = @_;

  my $name = $self->guess_chrom_fasta_name($chrom);

  return defined($name) ? $self->{_chroms}->{$name} : undef;
}

# $offet-based -- last base not included in return
# chr1,2,3 returns base at index 2 only
sub get_chr_substr
{
  my ($self, $chrom, $start, $end) = @_;

  my $name = $self->guess_chrom_fasta_name($chrom);

  if(!defined($name))
  {
    return undef;
  }

  my $offset = $self->{_first_base_index};

  return substr($self->{_chroms}->{$name}, $start-$offset, $end-$offset)
    or carp("Chromosome position out of bounds of " . $offset . "-" . 
            ($self->{_lengths}->{$name}+$offset-1));
}

sub get_chr_substr1
{
  my ($self, $chrom, $start, $end) = @_;

  my $name = $self->guess_chrom_fasta_name($chrom);

  if(!defined($name))
  {
    return undef;
  }

  return substr($self->{_chroms}->{$name}, $start-1, $end-$start+1)
    or carp("Chromosome position out of bounds of 1-" .
            ($self->{_lengths}->{$name}));
}

sub get_chr_substr0
{
  my ($self, $chrom, $start, $len) = @_;

  my $name = $self->guess_chrom_fasta_name($chrom);

  if(!defined($name))
  {
    return undef;
  }

  return substr($self->{_chroms}->{$name}, $start, $len)
    or carp("Chromosome position out of bounds of 0-" .
            ($self->{_lengths}->{$name}-1));
}

sub get_chr_length
{
  my ($self, $chrom) = @_;

  my $name = $self->guess_chrom_fasta_name($chrom);

  return defined($name) ? $self->{_lengths}->{$name} : undef;
}

sub get_chr_quality_scores
{
  my ($self, $chrom) = @_;

  my $name = $self->guess_chrom_fasta_name($chrom);

  return defined($name) ? $self->{_quals}->{$name} : undef;
}

sub get_chr_data
{
  my ($self, $chrom) = @_;

  my $name = $self->guess_chrom_fasta_name($chrom);

  if(defined($name))
  {
    return ($self->{_chroms}->{$name},
            $self->{_lengths}->{$name},
            $self->{_quals}->{$name});
  }
  else
  {
    return undef;
  }
}

1;
