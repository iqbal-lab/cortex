/*
  supernode.h defines the interface for dealing with supernodes
*/

#ifndef SUPERNODE_H_
#define SUPERNODE_H_

#include <element.h>


// 1.  dBSupernode "inherits" the same orientation as its nodes. Given the hashtable, each node has
//     a well-defined forward direction, and so, so does the Sopernode
// 2.  To save memory, the dBSupernode object contains only one Node, and a length. You can always use the get_next_node function in dB_graph.h to work
//     your way along the supernode.
// 3.  Max tells you how much space is allocated
typedef struct
{
  dBNode* first_node;
  int length;
  int max;
} dBSupernode;



typedef struct
{
  dBNode** nodes;
  int length;
  int max;
} dBNodeArray;


#endif /* SUPERNODE_H_ */
