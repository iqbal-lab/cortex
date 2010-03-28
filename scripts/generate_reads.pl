#this scripts generates "short-reads" from a fasta entry
#it generates paired reads (in two files _1 & _2)
use strict;
use lib "/homes/marioc/perl-modules/bioperl-live";
use Getopt::Long;
use  Bio::SeqIO;

my ($file,$read_length,$coverage,$insert_size,$delta, $error_rate,$length);

&GetOptions(
     	    'file|f:s'               => \$file,              #file contains one big fasta entry
	    'read|r:i'               => \$read_length,
	    'cov|c:i'                => \$coverage,
	    'insert|i:i'             => \$insert_size,
	    'delta|s:i'              => \$delta,             #range of insert size (ie real insert size will be between insert_size-delta .. insert_size+delta)
            'error|r:i'          => \$error_rate,         #1 error every error_rate bases in the final 20% of the read
	    'length|l:i'             => \$length,            #length of entry in file
           );


print "open file $file\n";
open(F,"$file") || die "cannot open file $file\n";

my $entry = <F>;
chomp($entry);
$entry =~s/>//;

my $erate_for_printing;
if ($error_rate>100000000000)
{
    $erate_for_printing=0;
}
else
{
    $erate_for_printing=$error_rate;
}

open(F1,">reads_${entry}_1_read_len_".$read_length."coverage_".$coverage."_error_rate_".$erate_for_printing_end) || die "cannot open file 1\n";
open(F2,">reads_${entry}_2_read_len_".$read_length."coverage_".$coverage."_error_rate_".$erate_for_printing_end) || die "cannot open file 2\n";


print "entry:$entry length:$length coverage:$coverage read length: $read_length\n";
  
my $number_reads = int(($coverage * $length)/$read_length);

print "generate $number_reads read (read length $read_length) for $coverage X coverage\n";

#generate positions for reads
my @reads;
my $i;
while($i<$number_reads){
  my $pos = int(rand($length-$read_length))+1;
  if ($pos+$insert_size+$delta<$length){
    $reads[$pos]++;
    $i+=2; #becuase we generate read pairs
  }
}
  

#print reads

#every how many reads to introduce errors
my $reads_rate = int($error_rate/(0.2*$read_length));
print "print reads - introduce errors in the end of a read every $reads_rate reads ($error_rate)\n";

my $count_reads=0;
my $index_error_1;
my $index_error_2;

my $pos_1;
my $pos_2;


#reads the input fasta entry en pair of buffers (of size 10000)
my $file_buf1 = read_file(10000);
my $file_buf2 = read_file(10000);
my $current_pos = 0;
my $file_buf = $file_buf1.$file_buf2;

for($pos_1=1;$pos_1<=@reads;$pos_1++){
  print "." if $pos_1 % 1000000 == 0;

  my $i;
  for($i=1;$i<=$reads[$pos_1];$i++){
    
    if ($pos_1>$current_pos+length($file_buf1)){
      $file_buf1 = $file_buf2;
      $file_buf2 = read_file(10000);
      if (length($file_buf2)>0){
	$file_buf = $file_buf1.$file_buf2;
	$current_pos+=length($file_buf1);
      }
    }

    #delta for insert size
    my $var = int(rand(2*$delta))-$delta;
    my $real_insert=$insert_size+$var;      
    my $pos_2 = $pos_1+$real_insert-$read_length;
    
    if ($count_reads==0){#define where - wich reads to introduce error
      $index_error_1 = int(rand($reads_rate+1));
      $index_error_2 = int(rand($reads_rate+1));
    }
    
    my $seq_1 = substr($file_buf,$pos_1-$current_pos-1,$read_length);
    #print "TEST ",$pos_1-$current_pos,"\n";
    my $suf1;
    if ($count_reads == $index_error_1){
      my $base;
      ($base,$seq_1) = change_base($seq_1);
      $suf1="*_$base";
    }
    
    print F1 ">read_${i}_${pos_1}_${real_insert}_1$suf1\n";
    print F1 $seq_1,"\n";
    
    my $seq_2 =  substr($file_buf,$pos_2-$current_pos-1,$read_length);
    #print "TEST ",$pos_2-$current_pos,"\n";

    my $suf2;
    if ($count_reads == $index_error_2){
      my $base;
      ($base,$seq_2) = change_base($seq_2);
      $suf2="*_$base";
    }
    
    print F2 ">read_${i}_${pos_1}_${real_insert}_2$suf2\n";
    print F2 $seq_2,"\n";
    
    $count_reads++;
    
    if ($count_reads>$reads_rate){
	$count_reads=0;
      }
  }
  
}
print "\n";


sub read_file{
  my ($n) = @_;

  my $seq;
  while(<F>){
    chomp;
    $seq.=$_;
    if (length($seq)>$n){
      last;
    }
  }
  
  return $seq;    
}


## Errors in the last 20% of a read
sub change_base{

  my ($str) = @_;

  my %num2base = (0 => 'A', 1 => 'C', 2 => 'G', 3 => 'T');

  my $pos = int(0.8*$read_length + int(rand(0.2*$read_length)));

  my @str = split //,$str;
  my $index;
  do {
    $index = int(rand(4));
  } until (uc $str[$pos] ne $num2base{$index});

  $str[$pos] = $num2base{$index};

  return ($pos+1,join("",@str));
}

