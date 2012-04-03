package GeneticsModule;

use strict;
use warnings;

use Carp;

use base 'Exporter';

our @EXPORT = qw(get_clean_chr_name
                 gc_content
                 rev_comp complement rev_comp_cyclic
                 dna_rev_comp_group dna_word_group
                 dna_pattern_group dna_weak_strong_group);


sub gc_content
{
  my ($seq) = @_;
  my $gc = 0;

  while($seq =~ /([gc]+)/gi)
  {
    $gc += length($1);
  }

  return $gc;
}

sub get_clean_chr_name
{
  #takes: chromosomeX chr2a 12_random etc
  #returns one of: chr1, chr2a, chr2b, chr3, ..., chr22 chrX, chrY, chrM
  # (plus some _suffix if present e.g. _random)
  my ($chr) = @_;
  
  # prioritise:
  # 1) chr... 2 digits or '2a/b'
  # 2) chr... 1 digit
  # 3) 2 digits or '2a/b'
  # 4) 1 digit
  if($chr =~ /(?:chromosome|chrom|chr)\s*(1[0-9]|2[012ab]?|[mxy]|Un)(_\w+)?(?:[^\d]|$)/i ||
     $chr =~ /(?:chromosome|chrom|chr)\s*(\d)(_\w+)?(?:[^\d]|$)/i || # look for 1 digit
     $chr =~ /(1[0-9]|2[012ab]?|[mxy]|Un)(_\w+)?(?:[^\d]|$)?/i || # favour 2 digits
     $chr =~ /(\d)(_\w+)?(?:[^\d]|$)/i) # look for 1 digit
  {
    my $code = $1;
    my $suffix = defined($2) ? $2 : '';
    
    if($code =~ /^[mxy]$/i)
    {
      return 'chr'.uc($code).lc($suffix);
    }
    elsif($code =~ /^\d+$/ && $code > 22)
    {
      carp "Warning: Cannot resolve chromosome '$chr' - " .
           "chr '$code' number too high\n";
      return $chr;
    }
    else
    {
      return 'chr'.lc($code.$suffix);
    }
  }
  
  print STDERR "get_clean_chr_name: Cannot resolve chromosome '$chr'\n";
  return $chr;
}

sub rev_comp
{
  my ($seq) = @_;

  # Reverse
  my $rev = reverse(complement($seq));

  return $rev;
}

sub complement
{
  my ($seq) = @_;

  my %mapping = ('a' => 't', 'A' => 'T',
                 'c' => 'g', 'C' => 'G',
                 'g' => 'c', 'G' => 'C',
                 't' => 'a', 'T' => 'A');
  
  my $complement = $seq;
  
  # Complement
  for(my $i = 0; $i < length($seq); $i++)
  {
    my $chr = $mapping{substr($seq,$i,1)};

    if(!defined($chr))
    {
      carp("rev_comp: Cannot complement base '".substr($seq,$i,1)."' " .
           "in string '$seq'");
    }

    substr($complement, $i, 1) = $chr;
  }

  return $complement;
}

sub dna_rev_comp_group
{
  my ($seq) = @_;

  my $rev = rev_comp($seq);

  return ($seq le $rev ? $seq : $rev);
}

sub rev_comp_cyclic
{
  my ($word) = @_;

  #print "cyclic: $word\n";

  $word = uc($word);
  # Get word reverse complement
  my $word_rev_comp = rev_comp($word);

  # Get all cyclic permutations of word and its rev. comp.
  my @words = ();
  
  for my $w ($word, $word_rev_comp)
  {
    for(my $i = 0; $i < length($word); $i++)
    {
      push(@words, substr($w,$i).substr($w,0,$i));
    }
  }

  @words = sort {$a cmp $b} @words;

  return @words;
}

sub dna_word_group
{
  my ($word) = @_;

  my @words = rev_comp_cyclic($word);

  return $words[0];
}

sub dna_pattern_group_internal
{
  my ($word) = @_;

  my @chars = split(//,$word);

  my %grammar = ();

  $grammar{$chars[0]} = 'a';
  $grammar{rev_comp($chars[0])} = 'A';

  # Look for the next undefined character
  for(my $i = 1; $i < @chars; $i++)
  {
    if(!defined($grammar{$chars[$i]}))
    {
      $grammar{$chars[$i]} = 'b';
      $grammar{rev_comp($chars[$i])} = 'B';
      last;
    }
  }
  
  my $pattern = join("", map {$grammar{$_}} @chars);

  return $pattern;
}

sub dna_pattern_group
{
  my ($word) = @_;

  my $pattern1 = dna_pattern_group_internal($word);
  my $pattern2 = dna_pattern_group_internal(rev_comp($word));
  
  return $pattern1 le $pattern2 ? $pattern1 : $pattern2;
}

sub dna_weak_strong_group
{
  my ($word) = @_;

  $word = uc($word);
  
  $word =~ s/A/W/g;
  $word =~ s/T/W/g;

  $word =~ s/C/S/g;
  $word =~ s/G/S/g;

  my $rev = reverse($word);

  return $word le $rev ? $word : $rev;
}

1;
