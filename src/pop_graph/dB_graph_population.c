/*
  dB_graph_population.c - implementation
 */

#include <stdlib.h>
#include <stdio.h>

#include <element.h>
#include <hash_table.h>
#include <dB_graph.h>
#include <dB_graph_population.h>
#include <pqueue_pop.h>


boolean db_graph_is_this_node_in_this_person_or_populations_graph(dBNode* node, EdgeArrayType type, int index)
{
 
  Edges edge_for_this_person_or_pop = get_edge_copy(*node, type, index);

  if (edge_for_this_person_or_pop == 0)
    {
      return false;
    }
  else
    {
      return true;
    }
 
}


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

void db_graph_traverse_specific_person_or_pop(void (*f)(HashTable*, Element *, EdgeArrayType, int, char**, int*),HashTable * hash_table, EdgeArrayType type, int index, char** for_test, int* index_for_test){

  printf("*** db_graph_population *** start db_graph_traverse_specific_person_or_pop\n");


  int i;
  for(i=0;i<hash_table->number_buckets;i++){
    pqueue_traverse_specific_person_or_pop(f,hash_table, &(hash_table->table[i]), type,index, for_test, index_for_test);
  }

  printf("*** db_graph_population *** end db_graph_traverse_specific_person_or_pop\n");


}





void print_supernode_for_specific_person_or_pop(HashTable* db_graph, dBNode * node,EdgeArrayType type, int index, char** for_test, int* index_for_test)
{
  printf("*** db_graph_population *** start print_supernode_for_specific_person_or_pop\n");

  FILE * fout;
  long count=0;
  
  char filename [200];
  if (count % 100000000 == 0)
    {
      int num = count / 100000000;
      
      if (count !=0)
	{
	  fclose(fout);
	}

      if (type == individual_edge_array)
	{
	  sprintf(filename,"out_nodes_kmer_%i_person_%i_subset_%i",db_graph->kmer_size,index,num);
	}
      else
	{
	  sprintf(filename,"out_nodes_kmer_%i_population_%i_subset_%i",db_graph->kmer_size,index,num);
	}
      //fprintf(stderr,"opening file %s\n",filename);
      fout = fopen(filename,"w");
    }
  
  printf("in db_graph_pop index before print call is %d\n",*index_for_test);
  db_graph_print_supernode_for_specific_person_or_pop(fout,node,db_graph, type,index, for_test,index_for_test);
  printf("QQQQQQQQQ in db_graph_pop index after print call is %d, latest supernode is %s\n",*index_for_test, for_test[*index_for_test]);

  printf("*** db_graph_population *** end print_supernode_for_specific_person_or_pop\n");


}



void db_graph_set_all_visited_nodes_to_status_none(dBGraph* hash_table)
{
  hash_table_traverse(&db_node_set_status_to_none, hash_table);
}

void db_graph_set_all_visited_nodes_to_status_none_for_specific_person_or_population(dBGraph* hash_table, EdgeArrayType type, int index)
{
  printf("not implemented yet");
  exit(1);
}

