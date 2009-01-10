/*
  dB_graph_population.c - implementation
 */

#include <element.h>
#include <hash_table.h>
#include <dB_graph.h>
#include <dB_graph_population.h>
#include <pqueue_pop.h>


// wrapper for hash_table_find, which allows you to look in the hash table
// specifically for nodes related to a specific person or population
// person or population  specified by which edge array type
// which person or pop specified by index

dBNode *  db_graph_find_node_restricted_to_specific_person_or_population(Key key, dBGraph * hash_table, EdgeArrayType type, int index)
{

  dBNode *  e = hash_table_find(key, hash_table);

  //ASSUMING read length is strictly greater than kmer-length
  //then you should never see a kmer (node) which is unconnected to another (even if the other is itself)
  //ie to check if this kmer is seen in a given person/pop, it is enough to check if it has an edge

  Edges edge_for_this_person_or_pop = get_edge_copy(*e, type, index);

  if (edge_for_this_person_or_pop == 0)
    {
      return NULL;
    }
  else
    {
      return e;
    }
  
}

void db_graph_traverse_specific_person_or_pop(void (*f)(Element *, EdgeArrayType, int),HashTable * hash_table, EdgeArrayType type, int index){
  int i;
  for(i=0;i<hash_table->number_buckets;i++){
    pqueue_traverse_specific_person_or_pop(f,&(hash_table->table[i]), type,index);
  }
}
