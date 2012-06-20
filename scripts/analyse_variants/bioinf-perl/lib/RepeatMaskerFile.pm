package RepeatMaskerFile;

use strict;
use warnings;

use Carp;

sub new
{
  my $class = shift;
  my $handle = shift;
  my $next_line = <$handle>;

  my $self = {
      _handle => $handle,
      _next_line => $next_line,
      _unread_entries => []
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

sub unread_entry
{
  my ($self, @entries) = @_;

  push(@{$self->{_unread_entries}}, @entries);
}

sub read_entry
{
  my ($self) = @_;

  if(@{$self->{_unread_entries}} > 0)
  {
    return pop(@{$self->{_unread_entries}});
  }

  my $line;
  while(defined($line = $self->read_line()) && $line !~ /^\s*\d+/)
  {
    $self->read_line();
  }

  if(!defined($line))
  {
    return undef;
  }

  # trim whitespace from ends
  $line =~ s/^\s+//;
  $line =~ s/\s+$//;

  my @cols = split(/\s+/, $line);

  if(@cols < 14)
  {
    print STDERR join(",", @cols)."\n";
    croak("Expected at leat 14 columns");
  }

  my @names = qw(sw_score perc_div perc_del perc_ins
                 query query_begin query_end query_left
                 strand repeat class
                 repeat_begin repeat_end repeat_left id
                 overlapping_domain);

  my %entry = ();
  @entry{@names} = @cols;

  return \%entry;
}

1;
