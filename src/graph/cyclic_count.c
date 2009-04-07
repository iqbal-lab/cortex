#include <element.h>
#include <open_hash/hash_table.h>


//shifts the least sign 2k bits cyclically.
//on some compilers the following may work 
//struct {
//int val: 12;
//} mystruct;
// compiler reats it as a 12 bit integer, so you can do bit shifts etc. Does not get allocated at=s 12 bits, but can treat it as such.

unsigned long long  rotate_least_sig_2k_bits(unsigned long long value, short kmer_size) 
{
  int two_k_minus_two = 2* kmer_size-2;

  int mask = (1<<2*kmer_size)-1;
  int x = mask & value;

  x = ((x & 3 )<<two_k_minus_two) | x >> 2;
  
  return x;
}

//marks node and cyclic permutations (that are in the graph) as visited
int db_node_how_many_cyclic_perms_of_this_node_are_in_graph(dBNode* node, HashTable* hash_table)
{

  int i;

 //go through all cyclic permuations of of binary_kmer, and see if they are in the graph;

  int count_of_shifts_that_are_in_the_graph;
  BinaryKmer kmer = element_get_kmer(node);



  for (i=0; i< hash_table->kmer_size; i++)
    {
      //rotate - ie shift cyclically - by 2 bits. ie one base
      kmer = rotate_least_sig_2k_bits(kmer,hash_table->kmer_size);

      Key key = element_get_key(kmer, hash_table->kmer_size);
      dBNode* query_node = hash_table_find(key,hash_table);
      if (query_node && !db_node_check_status(query_node,visited))
	{
	  //mark node as visited
	  db_node_set_status(query_node, visited);
	  count_of_shifts_that_are_in_the_graph++;
	}
    }

  return count_of_shifts_that_are_in_the_graph;

  
}


