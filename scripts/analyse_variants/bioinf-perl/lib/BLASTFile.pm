package BLASTFile;

use strict;
use warnings;

use Carp;

# Object methods:
# new read_blast_entries

use base 'Exporter';

our @EXPORT = qw(get_blast_summary);

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

sub peek_line
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

sub get_blast_summary
{
  my ($b) = @_;

  return 'Query='.$b->{'align_query_start'}.'-'.$b->{'align_query_end'}.','.
         'Sbjct='.$b->{'align_subject_start'}.'-'.$b->{'align_subject_end'}.
         $b->{'expect'};#.','.$b->{'align_query'}.','.$b->{'align_query'};
}

sub read_blast_entries
{
  my ($self, $expected_name) = @_;

  my $line;

  while(defined($line = $self->read_line()) && $line !~ /^Query=/i) {}

  if(!defined($line))
  {
    return ();
  }
  elsif(defined($expected_name) && $line !~ /$expected_name/i)
  {
    chomp($line);
    croak("BLAST entry does not match '$expected_name': '$line'");
  }
  
  my $blast_txt;

  my $peek;
  while(defined($peek = $self->peek_line()) && $peek !~ /^Query=/i)
  {
    $line = $self->read_line();
    $blast_txt .= $line;
  }

  return ($blast_txt, split_blast_entry($blast_txt));
}

sub split_blast_entry
{
  my ($blast_entry) = @_;

  my @lines = split(/[\n\r]/, $blast_entry);
  my @hits = ();
  my %chr_lengths = (); # record this because we can..

  my $curr_chr;

  for(my $i = 0; $i < @lines; $i++)
  {
    chomp($lines[$i]);

    if($lines[$i] =~ /^>.*\|(.*)$/) {
      $curr_chr = $1;
      $i++;
      # Look at next line
      
      for(; $i < @lines; $i++)
      {
        if($lines[$i] =~ /^Length=(\d+)/i)
        {
          $chr_lengths{$curr_chr} = $1;
          last;
        }
        else
        {
          $curr_chr .= $lines[$i];
        }
      }

      if(!defined($chr_lengths{$curr_chr}))
      {
        croak("BLAST processing error : $!");
      }
    }
    elsif($lines[$i] =~ /Score = (.*) bits \((\d+)\),  Expect = (.*)/i)
    {
      # New hit
      my %new_hit = ();
  
      $new_hit{'score'} = $1;
      $new_hit{'bits'} = $2;
      $new_hit{'expect'} = $3;
      
      $new_hit{'chr'} = $curr_chr;
      
      $i++; # Look at next line
      if($i < @lines &&
         $lines[$i] =~ /Identities = (\d+)\/(\d+) \((.*)%\), Gaps = (\d+)\/(\d+) \((.*)%\)/i)
      {
        $new_hit{'identity_nom'} = $1;
        $new_hit{'identity_denom'} = $2;
        $new_hit{'identity_percent'} = $3;
        $new_hit{'gaps_nom'} = $4;
        $new_hit{'gaps_denom'} = $5;
        $new_hit{'gaps_percent'} = $6;
      }
      else
      {
        my $err_line = $i < @lines ? $lines[$i] : '';
        chomp($err_line);
        croak("BLAST processing error ($i/".scalar(@lines).": '$err_line') : $!");
      }
      
      $i++; # Look at next line
      if($i < @lines &&
         $lines[$i] =~ /Strand=(Plus|Minus)\/(Plus|Minus)/i)
      {
        $new_hit{'query_strand'} = $1;
        $new_hit{'subject_strand'} = $2;
      }
      else {
        croak("BLAST processing error : $!");
      }
      
      # Read alignments
      $new_hit{'align_query'} = "";
      $new_hit{'align_subject'} = "";

      $i++;
      for(; $i < @lines; $i++)
      {
        #print " DEBUG: $lines[$i]\n";

        if($lines[$i] =~ /^\s*$/) {
          next;
        }
        elsif($lines[$i] =~ /^Query\s+(\d+)\s+([acgt\-]+)\s+(\d+)/i)
        {
          if(!defined($new_hit{'align_query_start'})) {
            $new_hit{'align_query_start'} = $1;
          }

          $new_hit{'align_query'} .= $2;
          $new_hit{'align_query_end'} = $3;
        }
        elsif($lines[$i] =~ /Sbjct\s+(\d+)\s+([acgt\-]+)\s+(\d+)/i)
        {
          if(!defined($new_hit{'align_subject_start'})) {
            $new_hit{'align_subject_start'} = $1;
          }

          $new_hit{'align_subject'} .= $2;
          $new_hit{'align_subject_end'} = $3;
        }
        elsif($lines[$i] !~ /^\s+\|+/)
        {
          # line is not alignment between query and subject
          # - no other lines allowed
          last;
        }
      }

      # Save new hit
      push(@hits, \%new_hit);

      # Backtract by one line
      $i--;

    } # end of line begining ' Score...' marking beginning of entry
  }

  return @hits;
}

1;
