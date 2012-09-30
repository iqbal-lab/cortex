package UsefulModule;

use strict;
use warnings;

use Carp;

# Inherit from the "Exporter" module which handles exporting functions.
# Most procedural modules make use of this.

use base 'Exporter';

# When the module is invoked, export, by default, the function "hello" into 
# the namespace of the using code.

our @EXPORT = qw(num2str mem2str round_int round_decimal
                 pretty_fraction binary_search_nearest trim
                 open_stdin);

=head1 NAME
 
num2str - Formats numbers so they are readable

=head1 SYNOPSIS

  use Num2Str;
  print num2str($num);
  print num2str($num, $separator);

=head1 DESCRIPTION

Format 

=head2 Functions

The following function is exported by default

=head3 hello

  print num2str(12443);
  print num2str(4523123.4521, " ");

Returns the number with a separator (',' by default) between every 3 digits.  
Example '12,443' and '4 523 123.4521'.  

=cut

sub num2str
{
  my $num = shift;
  my $sep = shift; # optional
  my $num_decimals = shift; # optional
  
  my @places = split(/\./, abs($num));
  my $intDigits = $places[0];

  # add commas or spaces to $intDigits

  my $initialPlaces = length($intDigits) % 3;
  
  if($initialPlaces == 0 && length($intDigits) >= 3)
  {
    $initialPlaces = 3;
  }
  
  my $strNum = substr($intDigits,0,$initialPlaces);
  
  for(my $i = $initialPlaces; $i <= length($intDigits)-3; $i += 3)
  {
    $strNum .= (defined($sep) ? $sep : ",") . substr($intDigits,$i,3);
  }
  
  if($num < 0)
  {
    $strNum = "-" . $strNum;
  }
  
  if(@places > 1)
  {
    if(!defined($num_decimals) || length($places[1]) <= $num_decimals)
    {
      $strNum .= "." . $places[1];
    }
    else
    {
      $strNum .= "." . substr($places[1], 0, $num_decimals);
    }
  }
  
  return $strNum;
}

sub mem2str
{
  my $num_bytes = shift;
  my $full_unit = shift;

  # Don't want more than double digits
  if($num_bytes > 2**38)
  {
    return "" . num2str($num_bytes / (2**40)) . " " .
           (defined($full_unit) ? "terabytes" : "TB");
  }
  elsif($num_bytes > 2**28)
  {
    return "" . num2str($num_bytes / (2**30)) . " " .
           (defined($full_unit) ? "gigabytes" : "GB");
  }
  elsif($num_bytes > 2**18)
  {
    return "" . num2str($num_bytes / (2**20)) . " " .
           (defined($full_unit) ? "megabytes" : "MB");
  }
  elsif($num_bytes > 2**8)
  {
    return "" . num2str($num_bytes / (2**10)) . " " .
           (defined($full_unit) ? "kilobytes" : "kB");
  }
}

sub pretty_fraction
{
  my ($nominator, $denominator, $places) = @_;

  if(!defined($places))
  {
    $places = 2;
  }

  if($denominator == 0)
  {
    return num2str($nominator) . " / " . num2str($denominator);
  }

  my $percent = sprintf("%.".$places."f", 100 * $nominator / $denominator);

  return num2str($nominator) . " / " . num2str($denominator) . " " .
         "(" . $percent . "%)";
}

sub round_int
{
  my ($num, $round_to) = @_;

  if(!defined($round_to))
  {
    $round_to = 1;
  }

  return int($num / $round_to + 0.5) * $round_to;
}

sub round_decimal
{
  my ($num, $decimal_places) = @_;

  my $multiply = 10**$decimal_places;

  return int($num * $multiply + 0.5) / $multiply;
}

# Returns index of nearest value
sub binary_search_nearest
{
  # Note rounds up: 1.5 is 'nearer' to 2 than 1

  my ($arr, $search_value, $lower_bound, $upper_bound) = @_;

  if(@$arr == 0)
  {
    die("binary_search_nearest cannot an search empty array for value " .
        "'$search_value'");
  }
  elsif(!defined($lower_bound))
  {
    $lower_bound = 0;
    $upper_bound = scalar(@$arr) - 1;
  }

  #print "$lower_bound,$upper_bound,".@$arr."\n";

  if($lower_bound + 1 == $upper_bound)
  {
    # This bit makes it nearest rather than pure search
    my $dist_lower = $search_value - $arr->[$lower_bound];
    my $dist_upper = $arr->[$upper_bound] - $search_value;

    return ($dist_lower < $dist_upper ? $lower_bound : $upper_bound);
  }

  my $middle = int(($lower_bound + $upper_bound) / 2);

  if($search_value > $arr->[$middle])
  {
    return binary_search_nearest($arr, $search_value, $middle, $upper_bound);
  }
  elsif($search_value < $arr->[$middle])
  {
    return binary_search_nearest($arr, $search_value, $lower_bound, $middle);
  }
  else
  {
    return $middle;
  }
}

sub trim
{
  my ($str) = @_;

  $str =~ s/^\s+//;
  $str =~ s/\s+$//;

  return $str;
}

sub open_stdin
{
  my ($error_msg) = @_;

  my $stdin_handle;

  # -p checks STDIN is connected to a pipe
  if(!(-p STDIN) || !open($stdin_handle, "<&=STDIN"))
  {
    croak(defined($error_msg) ? $error_msg : "Cannot open STDIN");
  }

  return $stdin_handle;
}

1;
