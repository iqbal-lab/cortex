package FASTNFile;

use strict;
use warnings;

use Carp;

use List::Util qw(min);

# Object methods:
# is_fastq read_next read_all peek_line read_line

use base 'Exporter';
our @EXPORT = qw(read_all_from_files estimate_fastq_size print_FASTA print_FASTQ);

sub new
{
  my ($class, $handle, $descriptor) = @_;
  my $next_line = <$handle>;
  my $is_fastq = -1;
  
  my $tmp_line = $next_line;
  chomp($tmp_line);
  
  if(defined($tmp_line) && length($tmp_line) > 0)
  {
    if($tmp_line =~ /^@/) {
      $is_fastq = 1;
    }
    elsif($tmp_line =~ /^>/) {
      $is_fastq = 0;
    }
    else {
      my $len = min(length($tmp_line),10);
      croak("Invalid first line '" . substr($tmp_line, 0, $len) . "...'");
    }
  }

  my $self = {
      _descriptor => $descriptor,
      _handle => $handle,
      _next_line => $next_line,
      _is_fastq => $is_fastq,
      _line_number => 1
  };

  bless $self, $class;
  return $self;
}

sub is_fastq
{
  my ($self) = @_;
  return $self->{_is_fastq};
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
  $self->{_line_number}++;
  
  return $temp_line;
}

sub read_next
{
  my ($self) = @_;
  
  return $self->{_is_fastq} ? $self->read_next_fastq()
                            : $self->read_next_fasta();
}

sub read_next_fasta
{
  my ($self) = @_;

  my $line;
  while(defined($line = $self->read_line()) && $line =~ /^\s*$/) {}

  if(!defined($line))
  {
    return (undef);
  }
  elsif($line !~ /^>/)
  {
    chomp($line);
    croak("Fasta file title line does not begin with '>': '$line'");
  }
  
  chomp($line);
  
  my $title = substr($line,1);
  my $seq = "";
  my $peek;

  while(defined($peek = $self->peek_line()) && $peek =~ /^[^>]/i)
  {
    $line = $self->read_line();
    chomp($line);
    $seq .= $line;
  }

  return ($title, $seq);
}

#
# Read next FASTQ entry from a file
sub read_next_fastq
{
  my ($self) = @_;

  my $line;
  
  # Ignore empty lines
  while(defined($line = $self->read_line()) && $line =~ /^\s*$/) {}

  if(!defined($line)) {
    return (undef, undef);
  }
  elsif($line !~ /^@/)
  {
    chomp($line);
    croak("Fastq file title line does not begin with @ ('".substr($line,0,50)."')");
  }
  
  chomp($line);

  my $title = substr($line,1);

  my $sequence = "";
  my $skip = "";
  my $quality = "";

  my $num_of_seq_lines = 0;

  while(defined($line = $self->read_line()) && $line !~ /^\+/i)
  {
    chomp($line);
    $sequence .= $line;
    $num_of_seq_lines++;
  }
  
  # $line is now FASTQ 'skip line' (e.g. "+readname" or just "+")
  chomp($line);
  $skip = $line;

  if($sequence eq "")
  {
    croak("FASTQ file ended early - expected sequence line");
  }
  elsif($sequence !~ /[acgtn]+/i)
  {
    croak("FASTQ file - invalid sequence line '".substr($sequence,0,50)."'");
  }

  for(my $i = 0;
      length($quality) < length($sequence) &&
      defined($line = $self->read_line());
      $i++)
  {
    chomp($line);
    $quality .= $line;
  }
  
  if($quality eq "")
  {
    croak("FASTQ file ended early - expected quality line");
  }
  elsif(length($sequence) != length($quality))
  {
    print STDERR '@'."$title\n";
    print STDERR "$sequence\n";
    print STDERR "$quality\n";

    croak("FASTQ file: length of sequence (".length($sequence).") does not match " .
          "length of quality line (".length($quality).") for read '$title' " .
          (defined($self->{_descriptor}) ? "in ".$self->{_descriptor}.":"
                                         : "on line ") . $self->{_line_number});
  }

  return ($title, $sequence, $skip, $quality);
}

#
# load all reads from this file
# returns hash of name -> (sequence,quality)
sub read_all
{
  my ($self, $reads_hashref, $fastq_qual_hashref) = @_;

  if(!defined($reads_hashref)) {
    $reads_hashref = {};
  }
  
  if(!defined($fastq_qual_hashref) && $self->{_is_fastq}) {
    $fastq_qual_hashref = {};
  }

  my ($name, $seq, $skip, $qual) = $self->read_next();

  while(defined($name))
  {
    $reads_hashref->{$name} = $seq;

    if(defined($self->{_is_fastq})) {
      $fastq_qual_hashref->{$name} = $qual;
    }
    
    ($name, $seq, $skip, $qual) = $self->read_next();
  }

  if($self->{_is_fastq}) {
    return ($reads_hashref, $fastq_qual_hashref);
  }
  else {
    return ($reads_hashref);
  }
}

#
# load all reads from FASTA/Q files passed as arguments
# returns hash of name -> (sequence,quality)
sub read_all_from_files
{
  my $reads_hashref = {}; # hashref
  my $qualities_hashref = {};

  for my $file (@_)
  {
    my $handle;

    if($file eq "-")
    {
      if(-p STDIN)
      {
        # STDIN is connected to a pipe
        open($handle, "<&=STDIN") or croak("Cannot read STDIN pipe");
      }
      else
      {
        croak("Cannot open STDIN to read fasta/fastq");
      }
    }
    else
    {
      open($handle, $file)
        or croak("Cannot open fasta/fastq file '$file'");
    }

    my $fastn_file = new FASTNFile($handle);

    $fastn_file->read_all($reads_hashref, $qualities_hashref);

    close($handle)
      or print STDERR "Cannot close file '$file'\n";
  }

  return ($reads_hashref, $qualities_hashref);
}

#
# Returns: (first read length, file size, est. num of reads, est. num of bases)
sub estimate_fastq_size
{
  my ($file) = @_;

  my $read_file = $file;

  $read_file =~ s/(.*\.gz)\s*$/gzip -dc < $1|/;

  # Read first two lines of file
  open(FASTQ, $read_file) or croak("Can't read file: '$read_file'");
  my $firstLine = <FASTQ>;
  my $secondLine = <FASTQ>;
  close(FASTQ);

  if(!defined($firstLine) || !defined($secondLine))
  {
    croak("Can't read first two lines of '$file'");
  }

  # Get file size
  my $file_size = -s $file;

  # Read length
  my $read_length = $secondLine;
  chomp($read_length);
  $read_length = length($read_length);

  # Number of bytes per read +1 for '+'
  my $singleReadBytes = length($firstLine) + 2*length($secondLine) + 1;
  my $num_of_reads = int($file_size / $singleReadBytes);
  my $num_of_bases = $num_of_reads * $read_length;

  return ($read_length, $file_size, $num_of_reads, $num_of_bases);
}

# Line wrap is optional (-1 means no wrap)
sub print_wrap
{
  my ($txt, $line_wrap, $out) = @_;

  my $txt_len = length($txt);

  for(my $i = 0; $i < $txt_len; $i += $line_wrap)
  {
    print $out substr($txt, $i, $line_wrap) . "\n";
  }
}

sub print_entry
{
  my ($self, $title, $seq, $qualities, $line_wrap, $out) = @_;

  if($self->is_fastq())
  {
    print_FASTQ($title, $seq, $qualities, $line_wrap, $out);
  }
  else
  {
    print_FASTA($title, $seq, $line_wrap, $out);
  }
}

sub print_FASTA
{
  my ($title, $seq, $line_wrap, $out) = @_;

  if(!defined($out))
  {
    open($out, ">-") or die("Couldn't open pipe to stdout");
  }

  print $out ">$title\n";

  if(defined($line_wrap) && $line_wrap > 0)
  {
    print_wrap($seq, $line_wrap, $out);
  }
  else
  {
    print $out "$seq\n";
  }
}

# Line wrap is optional (-1 means no wrap)
# if out not given, prints to STDOUT
sub print_FASTQ
{
  my ($title, $seq, $qualitites, $line_wrap, $out) = @_;

  if(!defined($out))
  {
    open($out, ">-") or die("Couldn't open pipe to stdout");
  }

  print $out '@'.$title."\n";

  if(defined($line_wrap) && $line_wrap > 0)
  {
    print_wrap($seq, $line_wrap, $out);
  }
  else
  {
    print $out "$seq\n";
  }

  print $out "+\n";

  if(!defined($qualitites))
  {
    $qualitites = "!" x length($seq);
  }

  if(defined($line_wrap) && $line_wrap > 0)
  {
    print_wrap($qualitites, $line_wrap, $out);
  }
  else
  {
    print $out "$qualitites\n";
  }
}

1;
