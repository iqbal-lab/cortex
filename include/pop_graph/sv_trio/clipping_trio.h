/*
  clipping_trio.h

  Mario - Only reason this is not in pop_graph_core, is that when you clip something from someone's graph, you need to set its status
  You end up needing n enum which has values pruned_from_person_n_and_person_m etc.

*/

#ifndef DB_GRAPH__POPULATION_H_
#define DB_GRAPH_POPULATION_H_

#include <element.h>
#include <hash_table.h>

int db_graph_clip_tip_for_specific_person_or_pop(dBNode * node, int limit,dBGraph * db_graph, EdgeArrayType type, int index);



#endif
