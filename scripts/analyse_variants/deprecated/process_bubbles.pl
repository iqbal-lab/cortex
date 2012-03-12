use strict;

use lib "Algorithm-NeedlemanWunsch-0.03/lib/";
use Algorithm::NeedlemanWunsch;


sub score_sub {
  if (!@_) {
    return -1; # gap penalty
  }
  
  return ($_[0] eq $_[1]) ? 1 : -1;
}

my $matcher = Algorithm::NeedlemanWunsch->new(\&score_sub);

my (@seq1,@seq2);

my @align1;
my @align2;


sub on_align {
  my ($a,$b) = @_;
  $align1[$a]=$b;
  $align2[$b]=$a;
  #print join(" ","align",$a,$b,$seq1[$a],$seq2[$b],$seq1[$a] ne $seq2[$b]?"*":""),"\n";
}

sub shift_a {
  my ($a) = @_;
  $align1[$a]=-1;
  #print join(" ","shift a",$a,$seq1[$a]),"\n";

}


sub shift_b {
  my ($b) = @_;
  $align2[$b]=-1;
  #print join(" ","shift b",$b,$seq2[$b]),"\n";

}


my $count=0;
my %seq;


my $file = shift;
my $prefix = shift; ## to add in front of var names to make them globally unique, otherwise all files contain var_1, var_2, etc
open(FILE, $file)||die("Cnnot open $file");



my $printed_at_start_of_var=0;
while(<FILE>){

    my $line = $_;
    my $var_name;
    if ($line =~ /branch\_(\d+)\_(1|2)/)
    {
	if ($printed_at_start_of_var==0)
	{
	    print "\n\nSTART NEW VAR\n";
	    $printed_at_start_of_var=1;
	}
	elsif ($printed_at_start_of_var==1)
        {
            $printed_at_start_of_var=0;
        }


	$var_name = $prefix."_var_".$1;
	my $which_branch=$2;
	print "$var_name branch $which_branch\n";
	my $a = <FILE>;
	chomp $a;
	$seq{$which_branch} = $a;
	$count++;
    }
    elsif ($line =~ /var\_(\d+)\_(\S+)_branch/)
    {
	$var_name = $prefix."_var_".$1;
	my $which = $2;

	if ($printed_at_start_of_var==0)
	{
	    print "\n\nSTART NEW VAR\n";
	    $printed_at_start_of_var=1;
	}
	elsif ($printed_at_start_of_var==1)
        {
            $printed_at_start_of_var=0;
        }




	my $which_num;
	if ($which eq "trusted")
	{
	    $which_num=1;
	}
	elsif ($which eq "variant")
	{
	    $which_num=2
	}
	else
	{
	    die("Unexpected.")
	}
	print "$var_name branch $which_num\n";
	my $a = <FILE>;
	chomp $a;
	

	$seq{$which_num} = $a;
	$count++;
    }
    else
    {

    }


  if ($count==2){
    #print $seq{1},"\n";
    #print $seq{2},"\n";

    @seq1 = split //,$seq{1};
    @seq2 = split //,$seq{2};
   

    my $score = $matcher->align(
				\@seq1,
				\@seq2,
				{
				 align => \&on_align,
				 shift_a => \&shift_a,
				 shift_b => \&shift_b,
				}
			       );
    print "FORWARD ALIGNMENT\n";
    my ($count_snps_f,$count_indels_f) = pretty_printing(\@seq1,\@seq2,\@align1,\@align2);
    print join(" ",$count_snps_f,$count_indels_f),"\n\n";

    if ( ($count_indels_f > @seq1/3) || ($count_indels_f > @seq2/3) ){
	print "REVERSE ALIGNMENT\n";
	@seq2 = split //, rev_comp($seq{2});
	
	my $score = $matcher->align(
	    \@seq1,
	    \@seq2,
	    {
		align => \&on_align,
		shift_a => \&shift_a,
		shift_b => \&shift_b,
				  }
	    );
      
	my ($count_snps_r,$count_indels_r) = pretty_printing(\@seq1,\@seq2,\@align1,\@align2);
      
	print join(" ",$count_snps_r,$count_indels_r),"\n";
	
	print "\n";
    }
    else
    {
	print "NO REVERSE ALIGNMENT\n";
    }
    $count=0;

  }
}
close(FILE);

sub rev_comp{
  my ($seq) = @_;

  my $r_seq = reverse($seq);
  $r_seq =~ tr/acgtACGT/tgcaTGCA/;
 # print join(" ",$seq,$r_seq),"\n";
  return $r_seq;
}


sub pretty_printing{

  my ($seq1,$seq2,$align1,$align2) = @_;

  #print join(" ","TEST",scalar(@{$seq1}),scalar(@{$seq2}),scalar(@{$align1}),scalar(@{$align2})),"\n";
  
  my ($count_snps,$count_indels)=(0,0);

  my ($count1,$count2) = (0,0);
    my @line1;
    my @line2;
    my @align;

    while($count1<scalar(@{$seq1}) || $count2<scalar(@{$seq2})){
      if ($align1->[$count1] == -1){
	 push @line1,$seq1->[$count1];
	 push @line2,"_";
	 push @align," ";
	 $count1++;
	 $count_indels++;
      }
      elsif ($align2->[$count2] == -1){
	push @line2,$seq2->[$count2];
	push @line1,"_";
	push @align," ";
	$count2++;
	$count_indels++;
      }
      elsif ($align1->[$count1] == $count2 && $align2->[$count2] == $count1){
	push @line1,$seq1->[$count1];
	push @line2,$seq2->[$count2];
	
	if ($seq1->[$count1] eq $seq2->[$count2]){ 
	  push @align,"|";
	}
	else{
	  push @align,"*";
	  $count_snps++;
	}
	$count1++;
	$count2++;
      }
      else{
	die;
      }
    }
  print "Br1:";
  print join("",@line1),"\n";
  print "    ";
  print join("",@align),"\n";
  print "Br2:";
  print join("",@line2),"\n";

  return ($count_snps,$count_indels);
}

