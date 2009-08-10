#include <element.h>
#include <stdio.h>
#include <stdlib.h>
#include <file_reader.h>
#include <dB_graph.h>
#include <cyclic_count.h>

int main(int argc, char **argv){

  FILE *fp_fnames;
  char filename[1000];
  int hash_key_bits;
  dBGraph * db_graph = NULL; 
  short kmer_size;
  int bucket_size;
  int action; //0 dump graph - 1 call SNPs
  int ctg_length = 200000;
  
   //command line arguments 
  fp_fnames= fopen(argv[1], "r");    //open file of file names
  kmer_size        = atoi(argv[2]); 
  hash_key_bits    = atoi(argv[3]);  //number of buckets: 2^hash_key_bits
  bucket_size      = atoi(argv[4]);
  action           = atoi(argv[5]);
  DEBUG            = atoi(argv[6]);
  

  fprintf(stderr,"Input file of filenames: %s - action: %i\n",argv[1],action);
  fprintf(stderr,"Kmer size: %d hash_table_size (%d bits): %d - bucket size: %d - total size: %qd\n",kmer_size,hash_key_bits,1 << hash_key_bits, bucket_size, ((long long) 1<<hash_key_bits)*bucket_size);

  //Create the de Bruijn graph/hash table
  db_graph = hash_table_new(hash_key_bits,bucket_size, 10,kmer_size);
  fprintf(stderr,"table created: %d\n",1 << hash_key_bits);

  int count_file   = 0;
  long long total_length = 0; //total sequence length

  //Go through all the files, loading data into the graph
  while (!feof(fp_fnames)){

    fscanf(fp_fnames, "%s\n", filename);
    
    long long seq_length = 0;
    count_file++;

    seq_length += load_binary_from_filename_into_graph(filename,db_graph);

    total_length += seq_length;
    
    fprintf(stderr,"\n%i kmers: %qd file name:%s seq:%qd total seq:%qd\n\n",count_file,hash_table_get_unique_kmers(db_graph),filename,seq_length, total_length);

    hash_table_print_stats(db_graph);

    //print mem status
    FILE* fmem=fopen("/proc/self/status", "r");
    char line[500];
    while (fgets(line,500,fmem) !=NULL){
      if (line[0] == 'V' && line[1] == 'm'){
	fprintf(stderr,"%s",line);
      }
    }
    fclose(fmem);
    fprintf(stderr,"************\n");
  }  
    
  
 

  
  switch (action){
  case 0 :
    printf("dumping graph %s\n",argv[7]);
    db_graph_dump_binary(argv[7],&db_node_check_nothing,db_graph);
    break;

  case 1 :
    printf("call SNPs\n");
    //db_graph_detect_snps(db_graph);
    break;

  case 2:
    printf("clip tips\n");
    db_graph_clip_tips(db_graph);
    break;

  case 3:
    printf("print supernodes\n");
    db_graph_print_supernodes(argv[7],ctg_length,db_graph); 
    break;

  case 4:
    printf("clip tips\n");
    db_graph_clip_tips(db_graph);
    printf("print supernodes\n");
    db_graph_print_supernodes(argv[7],ctg_length,db_graph); 
    break;

  case 5:
    printf("detect SNPs\n");
    //db_graph_detect_snps(db_graph);
    break;
    
  case 6:
    printf("count kmers\n");
    db_graph_print_coverage(db_graph);
    break;
  
  case 7:
    printf("clip tips\n");
    db_graph_clip_tips(db_graph);
    printf("detect SNPs\n");
    //db_graph_detect_snps(db_graph);
    break;


  case 8:    
    printf("clip tips\n");
    db_graph_clip_tips(db_graph);

    //printf("remove low coverage nodes\n");
    //db_graph_remove_low_coverage_nodes(1,db_graph);

    printf("print supernodes\n");
    db_graph_print_supernodes(argv[8],ctg_length,db_graph); 
 
    printf("dumping graph %s\n",argv[7]);
    db_graph_dump_binary(argv[7],&db_node_check_status_not_pruned,db_graph);
 
    break;

  case 9:
    printf("clip tips\n");
    db_graph_clip_tips(db_graph);
    printf("detect SNPs\n");
    //db_graph_detect_snps(db_graph);
    break;


    
  case 10:
    
    printf("remove low coverage nodes (<=1) \n");
    db_graph_remove_low_coverage_nodes(1,db_graph);

    printf("clip tips\n");
    db_graph_clip_tips(db_graph);

    printf("print supernodes\n");
    db_graph_print_supernodes(argv[8],ctg_length,db_graph); 
 
    printf("dumping graph %s\n",argv[7]);
    db_graph_dump_binary(argv[7],&db_node_check_status_not_pruned,db_graph);

    break;


  case 11:
    
    printf("remove low coverage nodes (<=4) \n");
    db_graph_remove_low_coverage_nodes(4,db_graph);

    //printf("clip tips\n");
    //db_graph_clip_tips(db_graph);

    //printf("smooth bubbles\n");
    //db_graph_smooth_bubbles(10,kmer_size*2,50,db_graph);
    
    printf("print supernodes\n");
    db_graph_print_supernodes(argv[8],ctg_length,db_graph); 
 
    printf("dumping graph %s\n",argv[7]);
    db_graph_dump_binary(argv[7],&db_node_check_status_not_pruned,db_graph);

    break;
  

  case 12:
    
    //printf("remove low coverage nodes (<=4) \n");
    //db_graph_remove_low_coverage_nodes(4,db_graph);

    //printf("clip tips\n");
    //db_graph_clip_tips(db_graph);

    //printf("smooth bubbles\n");
    //db_graph_smooth_bubbles(10,kmer_size*2,50,db_graph);
    
    //printf("print supernodes\n");
    //db_graph_print_supernodes(argv[8],ctg_length,db_graph); 
 
    //printf("dumping graph %s\n",argv[7]);
    //db_graph_dump_binary(argv[7],&db_node_check_status_not_pruned,db_graph);

    db_graph_detect_vars(200,1000,db_graph);
    break;
  

  case 13:
    
    printf("smooth bubbles\n");
    db_graph_smooth_bubbles(10,kmer_size*2,db_graph);
    
    printf("print supernodes\n");
    db_graph_print_supernodes(argv[8],ctg_length,db_graph); 
 
    printf("dumping graph %s\n",argv[7]);
    db_graph_dump_binary(argv[7],&db_node_check_status_not_pruned,db_graph);

    break;

  case 14:
    printf("Graph loaded. Count how many nodes have cyclic shifts\n");
    int results[kmer_size];
    int i;
    for (i=0; i<kmer_size; i++)
      {
	results[i]=0;
      }
    int* results_ptrs[kmer_size];
    for (i=0; i<kmer_size; i++)
    {
      results_ptrs[i]=&results[i];
    }
    
    //db_graph_traverse_with_array(&db_node_count_shifts_that_are_in_graph_and_log_in_array, db_graph, results_ptrs, kmer_size);

    for (i=0; i<kmer_size; i++)
      {
	printf("Number of times where a node matches %d other nodes by shifting =  %d\n", i, results[i]);
      }

  
  case 15:
    printf("Print supernodes containing entirely novel(non-reference) sequence\n");
    read_all_ref_chromosomes_and_mark_graph(db_graph);
    int min_covg=2; //demand coverage of at least 2

    //todo - add condition to make sure ignore pruned
    printf("Start printing..\n");
    db_graph_print_supernodes_where_condition_is_true_for_all_nodes_in_supernode(db_graph, &db_node_check_status_is_not_exists_in_reference,min_covg, stdout, false, NULL, 0);

    hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph); //cleanup - the printing process has set all printed nodes to visited
                                                                                                            //        - or to visited_and_exists_in_reference

  case 16:
    printf("Print supernodes that match reference at start and end, but not the middle\n");
    read_all_ref_chromosomes_and_mark_graph(db_graph);
    int minimum_covg=5; //demand covg of at least 5
    printf("Start printing..\n");
    db_graph_print_supernodes_where_condition_is_true_at_start_and_end_but_not_all_nodes_in_supernode(db_graph, &db_node_check_status_exists_in_reference, minimum_covg,
												      2,25,40, stdout, false,NULL,0);//2,25 are required numbers of overlapping nodes with ref at start and end
                                                                                                                             //40 is min number of supernodes which are NOT in reference
    //cleanup - the printing process has set all printed nodes to visited
    //           - or to visited_and_exists_in_reference
    hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph); 


  case 17:
    min_covg=2; //demand coverage of at least 2
    read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("/nfs/1000g-work/G1K/work/zi/projects/marzam/output/zi_20090716_hla_drb1/hla_drb1.fasta" , db_graph);
    printf("Loaded chromosome HLA DRB1\n");

    //    printf("Start printing supernodes that match HLA DRB1  exactly\n");
    //db_graph_print_supernodes_where_condition_is_true_for_all_nodes_in_supernode(db_graph, &db_node_check_status_exists_in_reference,min_covg, stdout, false, NULL, 0);
    //printf("STart printing supernodes that match HLA_DRB1 at start and end but not middle\n");
    //db_graph_print_supernodes_where_condition_is_true_at_start_and_end_but_not_all_nodes_in_supernode(db_graph, &db_node_check_status_exists_in_reference, minimum_covg,
    //                                                                                                      2,25,1, stdout, false,NULL,0);//2,25 are required numbers of overlapping nodes with ref at start and end, 1 is min nuber of nodes not in the reference


    printf("Start printing all supernodes that match the referemce (HLA DRB1) in any way\n");
    db_graph_print_supernodes_where_condition_is_true_for_at_least_one_node_in_supernode(db_graph, &db_node_check_status_exists_in_reference, minimum_covg,
                                                                                         stdout, false,NULL,0);



    hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph); //cleanup - the printing process has set all printed nodes to visited
    printf("Finished\n");


  case 18:
    printf("We have now loaded NA12878, Take first 40Mb of chromosome 2, and mark all supernodes that touch that. Then dump a binary.This is for debugging a segfault - want a small binary I can valgrind");
    read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference("/nfs/1000g-work/G1K/work/zi/projects/marzam/output/zi_20090807_fix_segfault/chrom2_first_40mb.fasta" , db_graph);

    int length_of_arrays = 100000;
    dBNode**     path_nodes        = (dBNode**) malloc(sizeof(dBNode*)*length_of_arrays); //everything these dBNode*'s point to will be in the hash table - ie owned by the hash table
    Orientation* path_orientations = (Orientation*) malloc(sizeof(Orientation)*length_of_arrays); 
    Nucleotide*  path_labels         = (Nucleotide*) malloc(sizeof(Nucleotide)*length_of_arrays);
    char*        path_string            = (char*) malloc(sizeof(char)*length_of_arrays+1); //+1 for \0

    if ( (path_nodes==NULL) || (path_orientations==NULL) || (path_labels==NULL) || (path_string==NULL) )
      {
	printf("Cannot even start. OOM. Failed to alloc arrays\n");
	exit(1);
      }

    int n;
    for (n=0; n<length_of_arrays; n++)
      {
	path_nodes[n]=NULL;
	path_orientations[n]=forward;
	path_labels[n]=Undefined;
	path_string[n]='N';
      }
    path_string[0]='\0';
    

    void mark_supernodes_of_ref_nodes_as_to_be_dumped(dBNode* node)
      {
	if ((node !=NULL) && (db_node_check_status(node,exists_in_reference)) )
	  {
	    int limit = 100000;
	    int min_covg=0;
	    int max_covg=0;
	    double avg_covg=0;
	    boolean is_cycle=false;

	    db_graph_supernode(node, limit, &db_node_action_set_status_to_be_dumped, 
			       path_nodes, path_orientations, path_labels, path_string,
			       &avg_covg, &min_covg, &max_covg, &is_cycle, db_graph);

	  }

      }

    FILE* binary_of_subsection_fptr = fopen("dumped_bin.ctx", "w");
    if ( binary_of_subsection_fptr ==NULL)
      {
	printf("Cannot open file to dump binary\n");
	exit(1);
      }

    void print_nodes_marked_as_to_be_dumped_as_binary(dBNode* elem)
      {
	if (db_node_check_status(elem, to_be_dumped))
	  {
	    db_node_print_binary(binary_of_subsection_fptr, elem);
	  }
      }

    hash_table_traverse(&mark_supernodes_of_ref_nodes_as_to_be_dumped, db_graph);

    //then dump everything to be dumped!
    hash_table_traverse(&print_nodes_marked_as_to_be_dumped_as_binary, db_graph);

    fclose(binary_of_subsection_fptr);

  }

  return 0;
}
