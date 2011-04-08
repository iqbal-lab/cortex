#this scripts generates "short-reads" from a fasta entry
#it generates paired reads (in two files _1 & _2)
use strict;
#use lib "/homes/marioc/perl-modules/bioperl-live";
use lib "/home/zam/dev/perl/installed_modules/bioperl-live";

use Getopt::Long;
use  Bio::SeqIO;

my ($file,$read_length,$coverage,$insert_size,$delta, $error_rate_main, $error_rate_end,$length);

&GetOptions(
     	    'file|f:s'               => \$file,              #file contains one big fasta entry
	    'read|r:i'               => \$read_length,
	    'cov|c:i'                => \$coverage,
	    'insert|i:i'             => \$insert_size,
	    'delta|s:i'              => \$delta,             #range of insert size (ie real insert size will be between insert_size-delta .. insert_size+delta)
            'error_main|m:i'          => \$error_rate_main,         #1 error every error_rate_main bases in the first 80% of the read
            'error_end|e:i'          => \$error_rate_end,         #1 error every error_rate_end bases in the final 20% of the read
            'length|l:i'             => \$length,            #length of entry in file
           );


print "open file $file\n";
open(F,"$file") || die "cannot open file $file\n";

my $entry = <F>;
chomp($entry);
$entry =~s/>//;
if ($entry =~ /^(\S+)/)
{
    $entry=$1;
}


my $erate_for_printing_main;
my $erate_for_printing_end;
if ($error_rate_main>100000000000)
{
    $erate_for_printing_main=0;
}
else
{
    $erate_for_printing_main=$error_rate_main;
}
if ($error_rate_end>100000000000)
{
    $erate_for_printing_end=0;
}
else
{
    $erate_for_printing_end=$error_rate_end;
}

open(F1,">reads_${entry}_1_read_len_".$read_length."coverage_".$coverage."_error_rate_main".$erate_for_printing_main."_error_rate_end".$erate_for_printing_end) || die "cannot open file 1\n";
open(F2,">reads_${entry}_2_read_len_".$read_length."coverage_".$coverage."_error_rate_main".$erate_for_printing_main."_error_rate_end".$erate_for_printing_end) || die "cannot open file 2\n";


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

#every how many reads to introduce errors, in the main body, and in the tail-end of the reads
my $reads_rate_main = int($error_rate_main/(0.8*$read_length));
my $reads_rate_end = int($error_rate_end/(0.2*$read_length));
print "print reads - introduce errors in the body (first 80%) of a read every $reads_rate_main reads ($error_rate_main)\n";
print "            - introduce errors in the end of a read every $reads_rate_end reads ($error_rate_end)\n";

if (($reads_rate_main==0) || ($reads_rate_end==0))
{
    die("Should not have 1 error every 0 reads");
}

my $count_reads_main = 0;#will keep count of reads, modulo the rate at which we want errors in main body
my $count_reads_end  = 0; # will keep count of reads, modulo the rate at which we want errors in the end
my $index_error_main_1;
my $index_error_main_2;
my $index_error_end_1;
my $index_error_end_2;

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
    
    if ($count_reads_main==0){#Specify which reads will put errors in in the main 80% of the read, in the next chunk of reads_rate_main reads
      $index_error_main_1 = int(rand($reads_rate_main+1));
      $index_error_main_2 = int(rand($reads_rate_main+1));
    }

    if ($count_reads_end==0){#Specify which reads will put errors in, at the tail end of the read,  in the next chunk of reads_rate_end reads
      $index_error_end_1 = int(rand($reads_rate_end+1));
      $index_error_end_2 = int(rand($reads_rate_end+1));
    }
    

    ## First sort out left hand mate read - 
    my $seq_1 = substr($file_buf,$pos_1-$current_pos-1,$read_length);
    #print "TEST ",$pos_1-$current_pos,"\n";
    my $suf1="";
    if ($count_reads_main == $index_error_main_1){
      my $base;
      ($base,$seq_1) = change_base_in_first_80percent($seq_1);
      $suf1="*_$base";
    }

    if ($count_reads_end == $index_error_end_1){
      my $base;
      ($base,$seq_1) = change_base_in_final_20percent($seq_1);
      $suf1= $suf1."*_$base";
    }
    
    print F1 ">read_${i}_${pos_1}_${real_insert}_1$suf1\n";
    print F1 $seq_1,"\n";
    
    my $seq_2 =  substr($file_buf,$pos_2-$current_pos-1,$read_length);
    #print "TEST ",$pos_2-$current_pos,"\n";

    my $suf2="";
    if ($count_reads_main == $index_error_main_2){
      my $base;
      ($base,$seq_2) = change_base_in_first_80percent($seq_2);
      $suf2="*_$base";
    }

    if ($count_reads_end == $index_error_end_2){
      my $base;
      ($base,$seq_2) = change_base_in_final_20percent($seq_2);
      $suf2= $suf2."*_$base";
    }
    
    print F2 ">read_${i}_${pos_1}_${real_insert}_2$suf2\n";
    print F2 $seq_2,"\n";
    
    $count_reads_main++;
    $count_reads_end++;
    
    if ($count_reads_main>$reads_rate_main){
	$count_reads_main=0;
      }

    if ($count_reads_end>$reads_rate_end){
	$count_reads_end=0;
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



sub change_base_in_final_20percent{

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


sub change_base_in_first_80percent{

  my ($str) = @_;

  my %num2base = (0 => 'A', 1 => 'C', 2 => 'G', 3 => 'T');

  my $pos = int(rand(0.8*$read_length));

  my @str = split //,$str;
  my $index;
  do {
    $index = int(rand(4));
  } until (uc $str[$pos] ne $num2base{$index});

  $str[$pos] = $num2base{$index};

  return ($pos+1,join("",@str));
}

