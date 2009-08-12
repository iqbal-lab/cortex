use strict;
use Getopt::Long;
use Bio::SeqIO;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

my $registry;
my $masked     = 0;
my $species;
my $start;
my $end;
my $id;

&GetOptions( 
	    'registry|r:s'   => \$registry,
	    'species|s:s'    => \$species,
	    'masked'         => \$masked,
	    'start:i'        => \$start,
	    'end:i'          => \$end,
	    'id:s'           => \$id,
	   );


Bio::EnsEMBL::Registry->load_all($registry);

#get dbg
my $db=  Bio::EnsEMBL::Registry->get_DBAdaptor($species,'core');
throw("database $species not found in registry\n") unless defined $db;

my $slice_adaptor = $db->get_SliceAdaptor;


my $seqout = Bio::SeqIO->new( '-format' => 'fasta',
                              '-fh'   => \*STDOUT);
    
print STDERR "looking up $id\n";

my $seq;

my $slice = $slice_adaptor->fetch_by_region('toplevel',$id,$start,$end);

if ($masked) {
  $seq = $slice->get_repeatmasked_seq([''],1);     
}    
else {
  $seq = $slice;
}    

die if not defined $seq;
$seqout->write_seq($seq);

