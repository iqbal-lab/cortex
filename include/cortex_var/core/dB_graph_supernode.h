/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo
 * 
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * **********************************************************************
 *
 * This file is part of CORTEX.
 *
 * CORTEX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CORTEX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CORTEX.  If not, see <http://www.gnu.org/licenses/>.
 *
 * **********************************************************************
 */
/*
  dB_graph_supernode.h defines the interface for dealing with supernodes
*/

#ifndef DB_GRAPH_SUPERNODE_H_
#define DB_GRAPH_SUPERNODE_H_

#include "global.h"
#include "element.h"


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


#endif /* DB_GRAPH_SUPERNODE_H_ */
