#include <element.h>
#include <stdio.h>
#include <stdlib.h>
#include <file_reader.h>
#include <dB_graph.h>
#include <cyclic_count.h>
#include <string.h>

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
    {
      printf("smooth bubbles\n");
      db_graph_smooth_bubbles(10,kmer_size*2,db_graph);
      
      printf("print supernodes\n");
      db_graph_print_supernodes(argv[8],ctg_length,db_graph); 
      
      printf("dumping graph %s\n",argv[7]);
      db_graph_dump_binary(argv[7],&db_node_check_status_not_pruned,db_graph);
      
      break;
    }
  case 14:
    {
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
      break;
    }
  case 15:
    {
      printf("Print supernodes containing entirely novel(non-reference) sequence\n");
      read_all_ref_chromosomes_and_mark_graph(db_graph);
      int min_covg=2; //demand coverage of at least 2
      
      //todo - add condition to make sure ignore pruned
      printf("Start printing..\n");
      db_graph_print_supernodes_where_condition_is_true_for_all_nodes_in_supernode(db_graph, &db_node_check_status_is_not_exists_in_reference,min_covg, stdout, false, NULL, 0);
      
      hash_table_traverse(&db_node_action_unset_status_visited_or_visited_and_exists_in_reference, db_graph); //cleanup - the printing process has set all printed nodes to visited
      //        - or to visited_and_exists_in_reference
      break;
    }
  case 16:
    {
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
      break;
    }
  case 17:
    {
      int minimum_covg=2; //demand coverage of at least 2
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
      break;
    }
  case 18:
    {
      //printf("print supernodes\n");
      //db_graph_print_supernodes(stdout,ctg_length,db_graph);

      printf("Given this list of fastq files %s, read them each and dump a fasta file for each one containing only those reads lying inside our (cleaned) graph\n", argv[9]);


      int file_reader(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
	* full_entry = true;
	
	if (new_entry!= true){
	  puts("new_entry has to be true for fastq\n");
	  exit(1);
	}
	
	return read_sequence_from_fastq(fp,seq,max_read_length);
      }
      

      FILE* list_fptr = fopen(argv[9], "r");
      if (list_fptr==NULL)
	{
	  printf("Cannot open %s\n", argv[9]);
	  exit(1);
	}
      
      char line[MAX_FILENAME_LENGTH+1];
      

      
      while(fgets(line,MAX_FILENAME_LENGTH, list_fptr) !=NULL)
	{
	  
	  //remove newline from end of line - replace with \0
	  char* p;
	  if ((p = strchr(line, '\n')) != NULL)
	    {
	      *p = '\0';
	    }

	  printf("Looking at %s\n", line);

	  char outputfile[200];
	  sprintf(outputfile,"%s.cleaned.fasta",line);
	  FILE* out = fopen(outputfile, "w");
	  if (out ==NULL)
	    {
	      printf("Cannot open %s, exiting", outputfile);
	      exit(1);
	    }
	  printf("We will output to %s\n", outputfile);


	  FILE* read_fptr = fopen(line, "r");
	  if (read_fptr==NULL)
	    {
	      printf("Cannot open %s", line);
	      exit(1);
	    }

	  long long b_reads=0; //will ignore base reads, but needed by API
	  int max_read_length=10000;
	  
	  read_fastq_and_print_reads_that_lie_in_graph(read_fptr, out, &file_reader, 
						       &b_reads, max_read_length, db_graph,
						       false, NULL, NULL);//last 3 arguments show is not for testing
	  fclose(out);
	  fclose(read_fptr);
	  
	}
      fclose(list_fptr);
      break;
    }
  case 19:
    {
      // We suspect that a couple of bubbles found by detect vars when Mario mixed ref and person in one graph are there
      // wrongly, because the edge between two specific nodes is only there in the reference, but the algorithm
      // acts as if it is there in the individual. Test these specific nodes and edges

      BinaryKmer tmp_kmer1;
      BinaryKmer tmp_kmer2;

      dBNode* node_of_fst_kmer_branch2_bubble_1434000 = hash_table_find(element_get_key(seq_to_binary_kmer("GTGTGGTGTGAATAGGCTCCTGCATGTGGCA", 
													   db_graph->kmer_size,&tmp_kmer1), db_graph->kmer_size, &tmp_kmer2),db_graph);
      Orientation fst_orientation = forward;
      //      dBNode* node_of_lst_kmer_branch2_bubble_1434000 = hash_table_find(element_get_key(seq_to_binary_kmer("TGTGGTGTGAATAGGCTCCTGCATGTGGCAT", 
      //		     										   db_graph->kmer_size,&tmp_kmer1), db_graph->kmer_size, &tmp_kmer2),db_graph);
      //Orientation lst_orientation = reverse;
      Nucleotide expected_nuc_joining_nodes = Thymine;

      if (node_of_fst_kmer_branch2_bubble_1434000==NULL)
	{
	  printf("Problem - get NULL when try to find nodes for bubble_1434000\n");
	}
      else
	{
	  if (db_node_edge_exist(node_of_fst_kmer_branch2_bubble_1434000, expected_nuc_joining_nodes , fst_orientation)==true)
	    {
	      printf("Damn, that edge does exist in NA12878 for bubble_1434000\n");
	    }
	  else
	    {
	      printf("Good. That edge for bubble_1434000 is NOT there is NA12878");
	    }
	}

      
      dBNode* node_of_fst_kmer_branch2_bubble_2293897 = hash_table_find(element_get_key(seq_to_binary_kmer("TCCGGTCTCCACGCCCCGCGGAGCAGGCCCA",
													   db_graph->kmer_size,&tmp_kmer1), db_graph->kmer_size, &tmp_kmer2),db_graph);
      fst_orientation = forward;
      expected_nuc_joining_nodes = Thymine;

      if  (node_of_fst_kmer_branch2_bubble_2293897==NULL)
        {
          printf("Problem - get NULL when try to find nodes for bubble_2293897\n");
        }
      else
        {
          if (db_node_edge_exist(node_of_fst_kmer_branch2_bubble_2293897, expected_nuc_joining_nodes , fst_orientation)==true)
            {
              printf("Damn, that edge does exist in NA12878 for bubble_2293897\n");
            }
          else
            {
              printf("Good. That edge for bubble_2293897 is NOT there is NA12878");
            }
        }


      dBNode* node_of_fst_kmer_branch2_bubble_1427650 = hash_table_find(element_get_key(seq_to_binary_kmer("AGTTAAAATGACAGACATGAAGCAGAATCAT",
													   db_graph->kmer_size,&tmp_kmer1), db_graph->kmer_size, &tmp_kmer2),db_graph);
      fst_orientation = forward;
      expected_nuc_joining_nodes = Thymine;

      if (node_of_fst_kmer_branch2_bubble_1427650==NULL) 
        {
          printf("Problem - get NULL when try to find nodes for bubble_1427650\n");
        }
      else
        {
          if (db_node_edge_exist(node_of_fst_kmer_branch2_bubble_1427650, expected_nuc_joining_nodes , fst_orientation)==true)
            {
              printf("Damn, that edge does exist in NA12878 for bubble_1427650\n");
            }
          else
            {
              printf("Good. That edge for bubble_1427650 is NOT there is NA12878");
            }
        }

      dBNode* node_of_lst_kmer_branch2_bubble_1427650 = hash_table_find(element_get_key(seq_to_binary_kmer("ATTTGGTATAACCTACTGATAGTCTTAGGTT",
													   db_graph->kmer_size,&tmp_kmer1), db_graph->kmer_size, &tmp_kmer2),db_graph);
      Orientation lst_orientation = forward;
      expected_nuc_joining_nodes = Guanine;

      if  (node_of_lst_kmer_branch2_bubble_1427650==NULL) 
        {
          printf("Problem - get NULL when try to find nodes for bubble_1427650\n");
        }
      else
        {
          if (db_node_edge_exist(node_of_lst_kmer_branch2_bubble_1427650, expected_nuc_joining_nodes , lst_orientation)==true)
            {
              printf("Damn, that edge does exist in NA12878 for bubble_1427650\n");
            }
          else
            {
              printf("Good. That edge for bubble_1427650 is NOT there is NA12878");
            }
        }

      

      
      break;
    }


  }

  return 0;
}
