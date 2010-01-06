#include <db_genotyping.h>
#include <element.h>
#include <open_hash/hash_table.h>
#include <dB_graph.h>
#include <dB_graph_population.h>
#include <seq.h>
#include <string.h>
#include <limits.h>
#include <file_reader.h>
#include <stdlib.h>


//Utility function - only exported so I can test it.
//Given two branches (arrays of dBNode*'s), return (in arguments 5,6,7,8) two arrays for each branch. One gives, for each node in the branch, the 
//number of times that node is seen in that branch. The second gives, for each node in a branch, the number of times that node is seen in
//the OTHER branch.
// The argument only_count_nodes_with_edge_in_specified_colour_func allows you to specify if you care whether a node has an edge in some colour or not.
// e.g. if we are only interested in the blue subgraph, and we do not want to count nodes that are not in the blue subgraph, we set this to true,
// and use the following two arguments also:
//The last two arguments allow you to do this to subgraph of the de Bruijn graph defined by these two functions - eg the union of all colours,
// or just one colour. These arguments are ignored if only_count_nodes_with_edge_in_specified_colour_func==false
void get_node_multiplicities(dBNode** branch1, int len_branch1, dBNode** branch2, int len_branch2, 
			     int** br1_multiplicity_in_br1, int** br2_multiplicity_in_br2,
			     int** br1_multiplicity_in_br2, int** br2_multiplicity_in_br1,
			     boolean only_count_nodes_with_edge_in_specified_colour_func,
			     Edges (*get_colour)(const dBNode*), int (*get_covg)(const dBNode*) )
{
  int working_array1[len_branch1];
  int working_array1_num=0; //number of things we have added to this array
  int working_array2[len_branch2];
  int working_array2_num=0; //number of things we have added to this array

  int k;

  void init_array(int* arr, int len)
  {
    for (k=0; k<len; k++)
    {
      arr[k]=-1;
    }
  }
  init_array(working_array1, len_branch1);
  init_array(working_array2, len_branch2);

  //we will reuse this;
  void get_mult(dBNode** br_src, int len_br_src, dBNode** br_target, int len_br_target, int** mult_array)
  {
    int i,j;
    int count_occurrences=0; //will be number of things we have put in this array

    for (i=0; i<len_br_src; i++)
      {
	*mult_array[i]=0;
      }

    for (i=0; i<len_br_src ; i++)
      {
	count_occurrences=0;
	
	for (j=0 ; j<len_br_target; j++)
	  {
	    //if i-th and j-th elements are the same, AND they exist in the colour (or function of colours) we are interested in
	    if (db_node_addr_cmp(&br_src[i], &br_target[j])==0 )
	      {
		if ( (only_count_nodes_with_edge_in_specified_colour_func==true) && 
		     (!db_node_is_this_node_in_subgraph_defined_by_func_of_colours(br_src[i], get_colour)) )
		  {
		    //does not count if node does not exist in the specified subgraph
		  }
		else
		  {
		    count_occurrences++;
		  }
	      }
	  }
	
	*mult_array[i]=count_occurrences;
      }
  }
  get_mult(branch1, len_branch1, branch1, len_branch1,  br1_multiplicity_in_br1);
  get_mult(branch2, len_branch2, branch2, len_branch2,  br2_multiplicity_in_br2);
  get_mult(branch1, len_branch1, branch2, len_branch2,  br1_multiplicity_in_br2);
  get_mult(branch2, len_branch2, branch1, len_branch1,  br2_multiplicity_in_br1);
  
}


