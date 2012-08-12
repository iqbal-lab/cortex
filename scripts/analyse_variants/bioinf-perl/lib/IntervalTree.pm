package IntervalTree;

use 5.006;
use strict;
use warnings;

use Carp;
use List::Util qw(min max);

our $VERSION = '0.01';

=head1 NAME
 
InternalTree - Quickly find intersecting / overlapping intervals

=head1 SYNOPSIS

  use IntervalTree;

  my $interval_tree = new IntervalTree(@intervals);

  my @hits = $intervals->fetch($search_position);
  my @hits = $intervals->fetch($search_start, $search_end);
  my ($hits_arr, $hits_left_arr, $hits_right_arr)
    = $intervals->fetch($search_start ,$search_end);

=head1 DESCRIPTION

Given a set of intervals, quickly find which intervals overlap a specified region.  
Ideal for bioinformatics (e.g. finding which annotations (repeats, genes etc.)
overlap a given region).  

See: http://www.dgp.toronto.edu/~jstewart/378notes/22intervals/ for a good
explanation of IntervalTree data structures.  

For dense interval sets, using IntervalTree is most effecient.
For sparse interval sets using IntervalList is probably faster / uses lower memory.  

=head1 EXAMPLE

  use IntervalTree;
  
  my @intervals = ([-12,3,'Zeroth'],
                   [10,20,'First'],
                   [15,45,'Second']);
  
  my $interval_tree = new IntervalTree(@intervals);

  my $search_start = 5;
  my $search_end = 25;

  my @hits = $interval_tree->fetch($search_start, $search_end);
  
  print "Hits ($search_start <= x < $search_end):\n";
  for my $hit (@hits) {
    print "$hit\n";
  }
  
  my ($hits_arr, $hits_left_arr, $hits_right_arr) = $interval_tree->fetch(5,25);

  print "Hits ($search_start <= x < $search_end):\n";
  for my $hit (@$hits_arr) {
    print "$hit\n";
  }
  print "Closest Left (x < $search_start):\n";
  for my $hit (@$hits_left_arr) {
    print "$hit\n";
  }
  print "Closest Right (x >= $search_end):\n";
  for my $hit (@$hits_right_arr) {
    print "$hit\n";
  }

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

  my $root;
  my @nodes = ();

  if(@intervals > 0)
  {
    # Find element boundaries, check input
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
    
    # sort boundaries
    @element_boundaries = sort {$a <=> $b} @element_boundaries;
  
    # find unique boundaries
    my @uniq_boundaries = ();
    my $last_boundary = $element_boundaries[0]-1;

    for my $boundary (@element_boundaries)
    {
      if($boundary > $last_boundary)
      {
        push(@uniq_boundaries, $boundary);
        $last_boundary = $boundary;
      }
    }

    #print "Boundaries: " . join(";", @uniq_boundaries) . "\n";

    # Now build tree
    # nodes: (boundary, left_node, right_node, index, interval_indices..)
    # leaf nodes: (undef, interval_indices..)
    my @left_boundaries = ();
    my @right_boundaries = ();

    $root = _build_tree(0, $#uniq_boundaries, \@uniq_boundaries,
                        \@left_boundaries, \@right_boundaries,
                        $uniq_boundaries[1], \@nodes);

    # Now add interval indices to the tree
    for(my $i = 0; $i < @intervals; $i++)
    {
      _add_interval_to_tree($root, $intervals[$i], $i,
                            \@left_boundaries, \@right_boundaries);
    }

    #print _tree_to_string($root) . "\n";
  }

  my $self = {_intervals => \@intervals,
              _root => $root,
              _nodes => \@nodes};

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

sub _tree_to_string
{
  my ($tree) = @_;

  if(!defined($tree->[0]))
  {
    # Leaf node
    my ($ignore,@intervals) = @$tree;
    return "[leaf:(".join(",",@intervals).")]";
  }
  else
  {
    # Internal node
    my ($boundary, $left_node, $right_node, $index, @intervals) = @$tree;

    return "[node$boundary:$index:(" . join(",",@intervals) . "):" .
           _tree_to_string($left_node) . ":" .
           _tree_to_string($right_node) . "]";
  }
}

sub _build_tree
{
  my ($left, $right, $boundaries,
      $node_left_boundaries, $node_right_boundaries, $prev_boundary,
      $nodes_arr) = @_;

  # Shifts to the right (5 -> 2:1:2, 4 -> 2:1:1) 
  my $middle = $right - int(($right - $left) / 2);

  #print "_build_tree $left:$middle:$right ($prev_boundary)\n";

  my $left_node;
  my $right_node;

  my $left_boundary = min($boundaries->[max($left-1,0)], $prev_boundary);
  my $right_boundary = max($boundaries->[min($right+1,@$boundaries-1)],
                           $prev_boundary);

  # Store boundaries
  $node_left_boundaries->[$middle] = $left_boundary;
  $node_right_boundaries->[$middle] = $right_boundary;

  my $num_of_nodes = $right - $left + 1;

  if($num_of_nodes == 1)
  {
    $left_node = [undef];
    $right_node = [undef];
  }
  elsif($num_of_nodes == 2)
  {
    #print "  Creating half node $middle) $boundaries->[$middle] " .
    #      "($left_boundary, $right_boundary)\n";
    $left_node = _build_tree($left, $left, $boundaries,
                             $node_left_boundaries, $node_right_boundaries,
                             $boundaries->[$middle], $nodes_arr);
    $right_node = [undef];
  }
  else
  {
    #print "Creating node $middle) $boundaries->[$middle] " .
    #      ($left_boundary, $right_boundary)\n";
    $left_node = _build_tree($left, $middle-1, $boundaries,
                             $node_left_boundaries, $node_right_boundaries,
                             $boundaries->[$middle], $nodes_arr);

    $right_node = _build_tree($middle+1, $right, $boundaries,
                              $node_left_boundaries, $node_right_boundaries,
                              $boundaries->[$middle], $nodes_arr);
  }

  # Create, store and return node
  my $new_node = [$boundaries->[$middle], $left_node, $right_node, $middle];
  $nodes_arr->[$middle] = $new_node;
  return $new_node;
}

sub _add_interval_to_tree
{
  my ($tree, $interval, $i, $node_left_boundaries, $node_right_boundaries) = @_;

  my ($interval_start, $interval_end) = @$interval;

  if(!defined($tree->[0]))
  {
    # This is a leaf node - store here
    #print "Adding to tree leaf\n";
    push(@$tree, $i);
  }
  else
  {
    my ($boundary, $left_node, $right_node, $index) = @$tree;
    
    my $left_boundary = $node_left_boundaries->[$index];
    my $right_boundary = $node_right_boundaries->[$index];
    
    if($interval_start <= $left_boundary && $interval_end >= $right_boundary)
    {
      # This is an internal node - store here
      #print "Adding to internal node\n";
      push(@$tree, $i);
    }
    else
    {
      if($interval_start < $boundary)
      {
        _add_interval_to_tree($left_node, $interval, $i,
                              $node_left_boundaries, $node_right_boundaries);
      }

      if($interval_end > $boundary)
      {
        _add_interval_to_tree($right_node, $interval, $i,
                              $node_left_boundaries, $node_right_boundaries);
      }
    }
  }
}

sub _find_recursive
{
  my ($tree, $start, $end, $interval_indices_hash_ref) = @_;

  #print "find $start,$end in "._tree_to_string($tree)."\n";

  if(!defined($tree->[0]))
  {
    # Leaf node - ignore first entry, rest are indices
    my ($ignore, @interval_indices) = @$tree;

    for my $interval_index (@interval_indices)
    {
      $interval_indices_hash_ref->{$interval_index} = 1;
    }
  }
  else
  {
    # Internal (aka boundary) node
    my ($boundary, $left_node, $right_node, $index, @intervals) = @$tree;

    # Add the intervals held by this internal node
    for my $interval (@intervals)
    {
      $interval_indices_hash_ref->{$interval} = 1;
    }

    if($start < $boundary)
    {
      _find_recursive($left_node, $start, $end, $interval_indices_hash_ref);
    }
    
    if($end >= $boundary)
    {
      _find_recursive($right_node, $start, $end, $interval_indices_hash_ref);
    }
  }
}

sub _walk_down_side_of_graph
{
  my ($tree, $intervals_arr, $down_left_side) = @_;

  # work down the left side of the graph
  # Intervals can't occur more than once down a single path of the graph
  # So we don't need to use a hash
  my @side_intervals = ();

  my ($boundary, $left_node, $right_node, $index, @intervals) = @$tree;
  push(@side_intervals, @intervals);
  $tree = ($down_left_side ? $left_node : $right_node);
  
  while(defined($tree->[0]))
  {
    # Grab intervals on internal node
    ($boundary, $left_node, $right_node, $index, @intervals) = @$tree;
    push(@side_intervals, @intervals);
    $tree = ($down_left_side ? $left_node : $right_node);
  }
    
  # Grab intervals on leaf
  ($boundary, @intervals) = @$tree;
  push(@side_intervals, @intervals);

  return map {$intervals_arr->[$_]->[2]} (sort {$a <=> $b} @side_intervals);
}

sub _find_nearest
{
  my ($tree, $search, $go_left, $intervals_in_region, $hits_arr) = @_;

  #print "nearest $go_left " . _tree_to_string($tree) . "\n";

  if(!defined($tree->[0]))
  {
    # Leaf node
    my ($ignore, @intervals) = @$tree;
    push(@$hits_arr, grep {!defined($intervals_in_region->{$_})} @intervals);
    return (@$hits_arr > 0);
  }
  else
  {
    my ($boundary, $left_node, $right_node, $index, @intervals) = @$tree;

    push(@$hits_arr, grep {!defined($intervals_in_region->{$_})} @intervals);
    
    if($search < $boundary)
    {
      if(_find_nearest($left_node, $search, $go_left,
                       $intervals_in_region, $hits_arr) ||
         (!$go_left && _find_nearest($right_node, $search, $go_left,
                                     $intervals_in_region, $hits_arr)))
      {
        #print "left/right!\n";
        return 1;
      }
    }
    else
    {
      if(_find_nearest($right_node, $search, $go_left,
                       $intervals_in_region, $hits_arr) ||
         ($go_left && _find_nearest($left_node, $search, $go_left,
                                    $intervals_in_region, $hits_arr)))
      {
        #print "right/left!\n";
        return 1;
      }
    }
    
    return (@$hits_arr > 0);
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

  #print "find ($start,$end)\n";

  my $tree = $self->{_root};
  my $nodes_arr = $self->{_nodes};
  my $intervals_arr = $self->{_intervals};

  if(!defined($tree))
  {
    # No intervals added
    #print "No tree :-(\n";
    return ([],[],[]);
  }

  # Get edges of extreme intervals
  my $left_boundary = $nodes_arr->[0]->[0];
  my $right_boundary = $nodes_arr->[@$nodes_arr - 1]->[0];

  my ($boundary, $left_node, $right_node, $index) = @$tree;

  my @hits = ();
  my @hits_left = ();
  my @hits_right = ();

  if($end < $left_boundary)
  {
    # work down LEFT side of the graph to get hits on the RIGHT of the region
    @hits_right = _walk_down_side_of_graph($tree, $intervals_arr, 1);
  }
  elsif($start >= $right_boundary)
  {
    # Work down RIGHT side of the graph to get hits on the LEFT of the region
    @hits_left = _walk_down_side_of_graph($tree, $intervals_arr, 0);
  }
  else
  {
    my %intervals_in_region = ();

    _find_recursive($tree, $start, $end, \%intervals_in_region);
    
    my @indices = keys %intervals_in_region;
    
    if($find_contained)
    {
      # Filter out intervals not contained by the search region
      @indices = grep {$intervals_arr->[$_]->[0] >= $start &&
                       $intervals_arr->[$_]->[1] <= $end} @indices;
    }
    
    @hits = map {$intervals_arr->[$_]->[2]} (sort {$a <=> $b} @indices);

    if($find_nearest)
    {
      my @hits_indices = ();    
      _find_nearest($tree, $start, 1, \%intervals_in_region, \@hits_indices);
      @hits_left =  map {$intervals_arr->[$_]->[2]} (sort {$a <=> $b} @hits_indices);

      @hits_indices = (); 
      _find_nearest($tree, $end, 0, \%intervals_in_region, \@hits_indices);
      @hits_right =  map {$intervals_arr->[$_]->[2]} (sort {$a <=> $b} @hits_indices);
    }
  }

  return (\@hits, \@hits_left, \@hits_right);
}

1;
