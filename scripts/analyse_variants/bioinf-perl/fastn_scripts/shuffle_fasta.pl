#!/usr/bin/perl

use strict;
use warnings;

use List::Util qw(shuffle);
use Fcntl qw{SEEK_SET};

# These two values are actually not needed...
#my $kmer_size = 55; # fixed for this project
#my $clean_binary
#  = "/data/zam/HLA/hlab/k55/B_nuc.fasta.fixed_readid.extended.k55.ctx";

# Fastq:
# /data/zam/projects/data/wtchg/na12878/fastq/fq_for_submission/

# HLA-B binary and info:
# /data/zam/HLA/hlab/k55/

# HLA-B binary built from fastq:
# /data/zam/HLA/hlab/B_nuc.fasta.fixed_readid.extended

#my $fastq_dir = "/ib/users/turner/cortexPaper/HLA/original_fastq/";
#my $save_dir = "/ib/users/turner/cortexPaper/HLA/fastq_files/";
#my $save_dir = "/ib/users/turner/cortexPaper/HLA/test_fastq/";

## Config
my $covg_per_file = 2;
my $genome_size = 3.1*10**9; # Genome 3.1 billion bases
my $num_of_output_files = 12;
##

## Testing
#$genome_size = 100; # only get 2 reads per output file, since reads are 101
#$num_of_output_files = 2;
##

if(@ARGV != 2)
{
  print "usage: ./shuffle_fasta.pl <fastq_dir> <save_dir>\n";
  print "  Sample without replacement from fastq files to create " .
        "$num_of_output_files new files\n";
  print "  Assume genome $genome_size bp, " .
        "covg per file needed = ".$covg_per_file."X\n";
  exit;
}

my $fastq_dir = shift;
my $save_dir = shift;

if($fastq_dir !~ /\/$/) { $fastq_dir .= "/"; }
if($save_dir !~ /\/$/) { $save_dir .= "/"; }

# Get fastq files
opendir(FASTQ_DIR, $fastq_dir) or die "Couldn't open dir '$fastq_dir': $!";
my @fastq_files = grep { /\.fastq$/i } readdir(FASTQ_DIR);
closedir(FASTQ_DIR);

if(@fastq_files == 0) {
  print STDERR "Couldn't find any fastq files!\n";
  exit;
}

my @shuffled_files = ();
my @read_counts = ();
my @handles = ();

for my $fastq_file (@fastq_files)
{
  my ($tmp_shuffled, $num_of_reads) = create_tmp_shuffled($fastq_dir,
                                                          $fastq_file);
  push(@shuffled_files, $tmp_shuffled);
  push(@read_counts, $num_of_reads);
}

for my $shuffled_file (@shuffled_files)
{
  my $handle;
  open($handle, $shuffled_file) or die("Cannot open file '$shuffled_file'");
  push(@handles, $handle);
}

# Write files
my $bases_per_file = $covg_per_file*$genome_size;
my $run_out_of_reads = 0;

my $file_num;

for($file_num = 0; $file_num < $num_of_output_files && !$run_out_of_reads;
    $file_num++)
{
  my $letter = chr(ord('a')+$file_num);
  my $sample_file = $save_dir.$covg_per_file."covg.$letter.fastq";

  print "Sampling reads into $sample_file\n";

  open(OUT, ">$sample_file") or die("Cannot open file for writing '$sample_file'");

  # Sample reads
  my $bases_read = 0;
  
  while($bases_read < $bases_per_file)
  {
    # sample a random read
    my $rand_index = int(rand() * @handles);
    my $rand_handle = $handles[$rand_index];
    
    # 4 lines per read
    my $header_line = <$rand_handle>;
    my $read_line = <$rand_handle>;
    my $spacer_line = <$rand_handle>;
    my $qual_line = <$rand_handle>;
    
    print OUT $header_line;
    print OUT $read_line;
    print OUT $spacer_line;
    print OUT $qual_line;
    
    chomp($read_line);
    $bases_read += length($read_line);
    
    $read_counts[$rand_index]--;
    
    if($read_counts[$rand_index] == 0)
    {
      # Remove file from @read_counts, @handles
      close($rand_handle) or print STDERR "Failed to close handle\n";
      
      print "Removing $rand_index\n";
      splice(@read_counts, $rand_index, 1);
      splice(@handles, $rand_index, 1);

      if(@read_counts == 0)
      {
        # Run out of reads
        print STDERR "Warning: ran out of reads\n";
        $run_out_of_reads = 1;
        last;
      }
    }
  }
  
  close(OUT) or print STDERR "Failure to close '$sample_file'\n";
}

# Close remaining handles
for my $handle (@handles)
{
  close($handle) or print STDERR "Failed to close handle ($!)\n";
}

print "Finished sampling reads\n";

# remove @shuffled_files
for my $shuffled_file (@shuffled_files)
{
  unlink($shuffled_file)
    or print STDERR "Cannot delete temporary file '$shuffled_file'\n";
}

## Now create fastalist files
print "Creating .fastqlist files\n";

# iterate through coverage levels required
for(my $i = 0; $i < $file_num; $i++)
{
  my $file = $save_dir.(($i+1)*$covg_per_file)."covg.fastqlist";

  open(OUT, ">$file") or die("Couldn't open '$file'");

  for(my $j = 0; $j <= $i; $j++) {
    print OUT $save_dir.$covg_per_file."covg.".chr(ord('a')+$j).".fastq\n";
  }
  
  close(OUT) or print STDERR "Couldn't close '$file'\n";
}

print "Done.\n";

# DONE!


sub create_tmp_shuffled
{
  my ($dir, $file) = @_;
  
  my $tmp_shuffled = $file.".tmp";
  
  if(-e $tmp_shuffled)
  {
    my $num = 1;
    while(-e $tmp_shuffled.$num) {
      $num++;
    }
    $tmp_shuffled .= $num;
  }
  
  my @read_positions = ();
  
  my $path = $dir.$file;
  
  print "$path: Loading positions...\n";
  
  open(FASTQ, $path) or die("Cannot read FASTQ file '$path'");
  
  my $line;
  my $position = 0;

  while(defined($line = <FASTQ>))
  {
    if($line =~ /^\@WTCHG/) {
      push(@read_positions, $position);
    }
    $position += length($line);
  }
  
  close(FASTQ);
  
  # Shuffle read positions...
  print "$path: Shuffling positions...\n";
  @read_positions = shuffle(@read_positions);
  
  print "$path: Saving positions...\n";
  
  open(TMP, ">$tmp_shuffled")
    or die("Cannot write temporary file '$tmp_shuffled'");
  
  open(FASTQ, $path) or die("Cannot read FASTQ file '$path'");
  
  for my $seek_pos (@read_positions)
  {
    if(!seek(FASTQ, $seek_pos, SEEK_SET))
    {
      die("Seek failure - seeking $seek_pos in file $path");
    }
    
    # 4 lines per read
    my $line;

    for(my $i = 0; $i < 4; $i++) {
      if(!defined($line = <FASTQ>)) {
        die("Unexpectedly ran out of lines ($path)");
      }
      print TMP $line;
    }
  }

  close(FASTQ);
  close(TMP);

  my $num_of_reads = @read_positions;

  return ($tmp_shuffled, $num_of_reads);
}
