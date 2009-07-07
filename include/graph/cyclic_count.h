
#ifndef CYCLIC_COUNT_H_
#define CYCLIC_COUNT_H_

#include <open_hash/hash_table.h>
#include <element.h>

unsigned long long  rotate_least_sig_2k_bits(unsigned long long value, short kmer_size);
int db_node_how_many_cyclic_perms_of_this_node_are_in_graph(dBNode* node, HashTable* hash_table);

#endif
