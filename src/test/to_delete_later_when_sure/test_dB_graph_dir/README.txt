
******  TEST 1.1 generates_graph_with_two_self_loops.fasta 

Contains this sequence:

>IL:blah blah 
ACGTAC

kmers therefore are:

ACG   CGT   GTA   TAC
TGC   GCA   CAT   ATG

Note those first two kmers are the same, and the last two kmers are also the same.
So in fact, the graph is a self loop out of kmer1 and back into itself, in the opposite orientation,
followed by an edge joining the two kmers, then a self loop out of the second kmer and back into itself
in the opposite orientation.

So if I ask, get me all the kmers (with orientation)  from GTA to the end of the supernode, I ought to get GTA/+ and GTA/-.
If I ask get me all the kmers (with orientation) from ACG to the end of the supernode, I should get ACG/+, ACG/-, GTA/+, GTA/-.

