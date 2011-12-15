require 'Descriptive.pm';
use Benchmark;

print "1..14\n";

$testct = 1;

# test #1

$stat = Statistics::Descriptive::Full->new();
@result = $stat->least_squares_fit();
print ( (@result? 'not ': '' ) . 'ok ' . $testct++ . "\n" );

# test #2
# data are y = 2*x - 1

$stat->add_data( 1, 3, 5, 7 );
@result = $stat->least_squares_fit();
$ok = ( $result[0] == -1 ) && ( $result[1] == 2 );
print ( ($ok? '': 'not ' ) . 'ok ' . $testct++ . "\n" );

# test #3
# test error condition on harmonic mean : one element zero
$stat = Statistics::Descriptive::Full->new();
$stat->add_data( 1.1, 2.9, 4.9, 0.0 );
$result = $stat->harmonic_mean();
$ok = ! defined( $result );
print ( ($ok? '': 'not ' ) . 'ok ' . $testct++ . "\n" );

# test #4
# test error condition on harmonic mean : sum of elements zero
$stat = Statistics::Descriptive::Full->new();
$stat->add_data( 1.0, -1.0 );
$result = $stat->harmonic_mean();
$ok = ! defined( $result );
print ( ($ok? '': 'not ' ) . 'ok ' . $testct++ . "\n" );

# test #5
# test error condition on harmonic mean : sum of elements near zero
$stat = Statistics::Descriptive::Full->new();
$savetol = $Statistics::Descriptive::Tolerance;
$Statistics::Descriptive::Tolerance = 0.1;
$stat->add_data( 1.01, -1.0 );
$result = $stat->harmonic_mean();
$ok = ! defined( $result );
print ( ($ok? '': 'not ' ) . 'ok ' . $testct++ . "\n" );

$Statistics::Descriptive::Tolerance = $savetol;

# test #6
# test normal function of harmonic mean
$stat = Statistics::Descriptive::Full->new();
$stat->add_data( 1,2,3 );
$result = $stat->harmonic_mean();
$ok = defined( $result ) && abs( $result - 1.6363 ) < 0.001;
print ( ($ok? '': 'not ' ) . 'ok ' . $testct++ . "\n" );

# test #7
# test stringification of hash keys in frequency distribution
$stat = Statistics::Descriptive::Full->new();
$stat->add_data(0.1,
                0.15,
                0.16,
               1/3);
%f = $stat->frequency_distribution(2);

$ok = ($f{0.216666666666667} == 3) &&
      ($f{0.333333333333333} == 1);
print ( ($ok? '': 'not ' ) . 'ok ' . $testct++ . "\n" );

# test #8
##Test memorization of last frequency distribution
%g = $stat->frequency_distribution();
$ok = 1;
foreach $key (keys %f) {
  $ok = 0 if $f{$key} != $g{$key};
}
print ( ($ok? '': 'not ' ) . 'ok ' . $testct++ . "\n" );

# test #9
# test the frequency distribution with specified bins
$stat = Statistics::Descriptive::Full->new();
@freq_bins=(20,40,60,80,100);
$stat->add_data(23.92,
                32.30,
                15.27,
                39.89,
                8.96,
                40.71,
                16.20,
                34.61,
                27.98,
                74.40);
%f = $stat->frequency_distribution(\@freq_bins);

$ok = ($f{20} == 3) &&
      ($f{40} == 5) &&
      ($f{60} == 1) &&
      ($f{80} == 1) &&
      ($f{100} == 0);
print ( ($ok? '': 'not ' ) . 'ok ' . $testct++ . "\n" );

# test #10 and #11
# Test the percentile function and caching
$stat = Statistics::Descriptive::Full->new();
$stat->add_data(-5,-2,4,7,7,18);
##Check algorithm
$ok =  ( $stat->percentile(50) == 4 );
print ( ($ok? '': 'not ' ) . 'ok ' . $testct++ . "\n" );
$ok =  ( $stat->percentile(25) == -2 );
print ( ($ok? '': 'not ' ) . 'ok ' . $testct++ . "\n" );

# tests #12 and #13
# Check correct parsing of method parameters
$stat = Statistics::Descriptive::Full->new();
$stat->add_data(1,2,3,4,5,6,7,8,9,10);
$ok = ($stat->trimmed_mean(0.1,0.1) == $stat->trimmed_mean(0.1));
print ( ($ok? '': 'not ' ) . 'ok ' . $testct++ . "\n" );

$ok = ($stat->trimmed_mean(0.1,0) == 6);
print ( ($ok? '': 'not ' ) . 'ok ' . $testct++ . "\n" );

# tests #14
# Make sure that trimmed_mean caching works but checking execution times
# This test may fail on very fast machines but I'm not sure how to get better
# timing without requiring extra modules to be added.

$stat = Statistics::Descriptive::Full->new();
##Make this a really big array so that it takes some time to execute!
$stat->add_data((1,2,3,4,5,6,7,8,9,10,11,12,13) x 10000);

my ($t0,$t1,$td);
my @t = ();
foreach (0..1) {
  $t0 = new Benchmark;
  $stat->trimmed_mean(0.1,0.1);
  $t1 = new Benchmark;
  $td = timediff($t1,$t0);
  push @t, $td->cpu_p();
}

$ok = $t[1] < $t[0];
print ( ($ok? '': 'not ' ) . 'ok ' . $testct++ . "\n" );
