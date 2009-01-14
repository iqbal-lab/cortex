/*
  dB_graph_population.h

  wrapper for hash_table for people who want to think of the hash table as 
  many de bruijn graphs drawn on top of each other.

*/

#ifndef DB_GRAPH__POPULATION_H_
#define DB_GRAPH_POPULATION_H_

#include <element.h>
#include <hash_table.h>

boolean db_graph_is_this_node_in_this_person_or_populations_graph(dBNode* node, EdgeArrayType type, int index);
dBNode *  db_graph_find_node_restricted_to_specific_person_or_population(Key key, dBGraph * hash_table, EdgeArrayType type, int index);
void db_graph_traverse_specific_person_or_pop(void (*f)(HashTable*, Element *, EdgeArrayType, int),HashTable * hash_table, EdgeArrayType type, int index);
void print_supernode_for_specific_person_or_pop(HashTable* db_graph, dBNode * node,EdgeArrayType type, int index);

#endif
