package IntervalList;

use 5.006;
use strict;
use warnings;

use Carp;

our $VERSION = '0.01';

=head1 NAME
 
InternalTree - Quickly find intersecting / overlapping intervals

=head1 SYNOPSIS

  use IntervalList;

  # Create a set with two intervals, each associated with an object
  my @interval_set = ([10,15,$item1], [20,25,$item2]);

  # Create intervals object from array of intervals
  my $intervals = new IntervalList(@intervals);

  # Different search methods
  my @hits = $intervals->fetch($position);
  my @hits = $intervals->fetch($start, $end);
  my ($hits_arr, $hits_left_arr, $hits_right_arr) = $intervals->fetch($start, $end);

  # Find regions contained within the search region
  my @hits = $intervals->fetch_contained($start, $end);

=head1 DESCRIPTION

Given a set of intervals, quickly find which intervals overlap a specified region.  
Ideal for bioinformatics (e.g. finding which annotations (repeats, genes etc.)
overlap a given region).  

See: http://www.dgp.toronto.edu/~jstewart/378notes/22intervals/ for a good
explanation of IntervalTree data structures.  

For dense interval sets, using IntervalTree is most effecient.
For sparse interval sets using IntervalList is probably faster / uses lower memory. 

=head2 Example

  use IntervalList;
  
  my @intervals = ([-12,3,'Zeroth'],
                   [10,20,'First'],
                   [15,45,'Second'],
                   [60,100,'Third']);
  
  my $interval_tree = new IntervalList(@intervals);

  my $search_start = 5;
  my $search_end = 25;

  my @hits = $interval_tree->fetch($search_start, $search_end);
  
  print "Hits ($search_start <= x <= $search_end):\n";
  for my $hit (@hits) {
    print " $hit\n";
  }

  # == Output ==
  # Hits ($search_start <= x <= $search_end):
  #  First
  #  Second
  
  my ($hits_arr, $hits_left_arr, $hits_right_arr)
    = $interval_tree->fetch($search_start, $search_end);

  print "Hits ($search_start <= x <= $search_end):\n";
  for my $hit (@$hits_arr) {
    print "$hit\n";
  }
  print "Closest Left (x < $search_start):\n";
  for my $hit (@$hits_left_arr) {
    print "$hit\n";
  }
  print "Closest Right (x > $search_end):\n";
  for my $hit (@$hits_right_arr) {
    print "$hit\n";
  }

  # == Output ==
  # Hits (5 <= x <= 25):
  #  First
  #  Second
  # Closest Left (x < 5):
  #  Zeroth
  # Closest Right (x > 25):
  #  Third

=head1 AUTHOR

Isaac Turner, C<< <turner.isaac at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests by emailing the author.  

=head1 LICENSE AND COPYRIGHT

Copyright 2011 Isaac Turner.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.

=cut

sub new
{
  my ($class, @intervals) = @_;

  # Sort by interval start position then by end position
  @intervals = sort _sort_intervals @intervals;

  my @element_boundaries = ();

  for my $interval (@intervals)
  {
    my ($start, $end) = @$interval;

    if(!defined($start) || !defined($end))
    {
      croak("Intervals must be an array of [start,end,object] entries");
    }
    elsif($start == $end)
    {
      croak("Interval start should not be equal to the end ($start <= x < $end)");
    }
    elsif($end < $start)
    {
      croak("Interval end should not be less than the start ($start, $end)");
    }
    
    push(@element_boundaries, $start, $end);
  }

  @element_boundaries = sort {$a <=> $b} @element_boundaries;

  my @elements = ();

  for(my $i = 0; $i < @element_boundaries; $i++)
  {
    my $start = $element_boundaries[$i];

    #print "Creating element $i starting at $start\n";

    push(@elements, [$start]);

    # Skip over matching entries
    while($i+1 < @element_boundaries && $element_boundaries[$i+1] == $start)
    {
      $i++;
    }
  }

  # Now loop through intervals, adding them to elements
  for(my $i = 0; $i < @intervals; $i++)
  {
    _add_interval_to_elements($intervals[$i], $i, \@elements);
  }

  my $self = {_elements => \@elements,
              _intervals => \@intervals};

  bless $self, $class;
  return $self;
}

sub fetch
{
  my ($self, $start, $end) = @_;

  return unless defined wantarray;  # void context, do nothing

  # Don't find nearest
  # Don't find contained
  my ($results_ref) = $self->_find($start, $end, 0, 0);

  return @$results_ref;
}

sub fetch_nearest
{
  my ($self, $start, $end) = @_;

  return unless defined wantarray;  # void context, do nothing

  # Find nearest
  # Don't find contained
  return $self->_find($start, $end, 1, 0);
}

sub fetch_contained
{
  my ($self, $start, $end) = @_;

  return unless defined wantarray;  # void context, do nothing

  # Don't find nearest
  # Find contained
  my ($results_ref) = $self->_find($start, $end, 0, 1);

  return @$results_ref;
}


## Interval Functions ##

sub _sort_intervals
{
  my $cmp = ($a->[0] <=> $b->[0]);

  if($cmp != 0)
  {
    return $cmp;
  }
  
  return $a->[1] <=> $b->[1];
}

sub _add_interval_to_elements
{
  my ($interval, $interval_i, $elements) = @_;

  my ($int_start, $int_end) = @$interval;

  # Binary search to find start
  my $left = 0;
  my $right = scalar(@$elements) - 1;
  
  my $middle = int(($left + $right) / 2);

  while($elements->[$middle]->[0] != $int_start)
  {
    if($left+1 == $right)
    {
      $middle = $right;
      last;
    }
    elsif($elements->[$middle]->[0] < $int_start)
    {
      $left = $middle;
    }
    else
    {
      $right = $middle;
    }

    $middle = int(($left + $right) / 2);
  }

  for(my $i = $middle; $i < @$elements && $elements->[$i]->[0] < $int_end; $i++)
  {
    #print "Adding interval $interval_i to element $i\n";
    push(@{$elements->[$i]}, $interval_i);
  }
}

sub _find
{
  my ($self, $start, $end, $find_nearest, $find_contained) = @_;

  if(!defined($end))
  {
    $end = $start;
  }
  elsif($start > $end)
  {
    croak("Start cannot be greater than end ($start, $end)");
  }

  my $elements = $self->{_elements};

  if(@$elements < 2)
  {
    # No intervals added
    # (min entry is (interval0,endpos))
    return ([],[],[]);
  }
  if($end < $elements->[0]->[0])
  {
    # Return those on the right (intervals in the first element)
    my @intervals_in_element = @{$elements->[0]};

    my @hits = map {$self->{_intervals}->[$_]->[2]}
               (sort {$a <=> $b} @intervals_in_element[1..$#intervals_in_element]);

    return ([],[],\@hits);
  }
  elsif($start >= $elements->[@$elements-1]->[0])
  {
    # Return those on the left (intervals in the last element)
    # Last element just remembers where the end position is...
    # need to use the second last element
    my @intervals_in_element = @{$elements->[@$elements - 2]};

    my @hits = map {$self->{_intervals}->[$_]->[2]}
               (sort {$a <=> $b} @intervals_in_element[1..$#intervals_in_element]);

    return ([],\@hits,[]);
  }

  # Start on the left, work right
  my $first_element = 0;

  if(@$elements > 2 && $start > $elements->[1]->[0])
  {
    # Binary search to find start
    my $left = 0;
    my $right = scalar(@$elements) - 1;
  
    my $middle = int(($left + $right) / 2);

    while(1)
    {
      if($elements->[$middle]->[0] <= $start &&
         $elements->[$middle+1]->[0] > $start)
      {
        last;
      }
      elsif($left+1 == $right)
      {
        $middle = $right;
        last;
      }
      elsif($elements->[$middle]->[0] > $start)
      {
        $right = $middle;
      }
      else
      {
        $left = $middle;
      }

      $middle = int(($left + $right) / 2);
    }
  
    $first_element = $middle;
  }

  #print "Start is $first_element\n";

  # Get intervals
  my %hit_interval_indices = ();

  my $element_i;

  for($element_i = $first_element;
      $element_i < @$elements-1 && $elements->[$element_i]->[0] <= $end;
      $element_i++)
  {
    my @intervals_in_element = @{$elements->[$element_i]};

    # Start at 1 since index 0 is the start position
    # positions 1.. are all interval indices
    for(my $i = 1; $i < @intervals_in_element; $i++)
    {
      $hit_interval_indices{$intervals_in_element[$i]} = 1;
      #print "Found $element_i => $intervals_in_element[$i]\n";
    }
  }
  
  my $next_element = $element_i;

  # Get indices of intervals in the search region
  my @indices = keys %hit_interval_indices;
  
  if($find_contained)
  {
    # Filter out intervals not contained by the search region
    @indices = grep {$self->{_intervals}->[$_]->[0] >= $start &&
                     $self->{_intervals}->[$_]->[1] <= $end} @indices;
  }

  my @hits = map {$self->{_intervals}->[$_]->[2]} (sort {$a <=> $b} @indices);

  my @left_hits = ();
  my @right_hits = ();

  if($find_nearest)
  {
    my $left_element = $first_element-1;

    # If we have hits, we may not have hits in the element directly to the left
    if(@hits > 0 && $left_element > 0 && @{$elements->[$left_element]} == 1)
    {
      # No elements in the bin to the left, shift left again
      $left_element--;
    }
  
    if($left_element >= 0)
    {
      # Get hits to the left
      my %interval_indices = ();

      my @intervals_in_element = @{$elements->[$left_element]};

      for(my $i = 1; $i < @intervals_in_element; $i++)
      {
        if(!defined($hit_interval_indices{$intervals_in_element[$i]}))
        {
          $interval_indices{$intervals_in_element[$i]} = 1;
        }
      }
      
      @left_hits = map {$self->{_intervals}->[$_]->[2]}
                   (sort {$a <=> $b} keys %interval_indices);
    }
    
    my $right_element = $next_element;
    
    # If we have hits, we may not have hits in the element directly to the right
    if(@hits > 0 && $right_element < @$elements-2 &&
       @{$elements->[$right_element]} == 1)
    {
      # No elements in the bin to the left, shift left again
      $left_element--;
    }
    
    if($right_element < @$elements-1)
    {
      # Get hits to the right
      my %interval_indices = ();

      my @intervals_in_element = @{$elements->[$right_element]};

      for(my $i = 1; $i < @intervals_in_element; $i++)
      {
        if(!defined($hit_interval_indices{$intervals_in_element[$i]}))
        {
          $interval_indices{$intervals_in_element[$i]} = 1;
        }
      }
      
      @right_hits = map {$self->{_intervals}->[$_]->[2]}
                   (sort {$a <=> $b} keys %interval_indices);
    }
  }

  return (\@hits, \@left_hits, \@right_hits);
}

1;
