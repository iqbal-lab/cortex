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

#include <db_differentiation.h>
#include <string.h>


void align_list_of_fastaq_to_graph_and_print_coverages_in_all_colours(FileFormat format, char* list_of_fastaq, int max_read_length, int* array_of_colours, char** array_of_names_of_colours,
								      int num_of_colours, dBGraph* db_graph,int fastq_ascii_offset,
								      boolean is_for_testing, char** for_test_array_of_strings, int* for_test_index)
{

  if ( (format != FASTA) && (format !=FASTQ) )
    {
      printf("Calling align_list_of_fastaq_to_graph_and_print_coverages_in_all_colours with file format not set to fasta or fastq\n");
      exit(1);
    }

  //For each file in list_of_fasta, go through the reads, and for each read, print one  "coverage read" per colour (space separated)
  // e.g. for a read print
  //    >read_id colour 0
  //    coverages of each of the nodes in the ref (colour 0)
  //   >read_id colour 1
  //     ... covgs in colour 1
  //   >read_id colour 2
  //     ... 
	
  
  //----------------------------------
  // allocate the memory used to read the sequences
  //----------------------------------
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  alloc_sequence(seq,max_read_length,LINE_MAX);
  
  //We are going to load all the bases into a single sliding window 
  KmerSlidingWindow* kmer_window = malloc(sizeof(KmerSlidingWindow));
  if (kmer_window==NULL)
    {
      printf("Failed to malloc kmer sliding window in align_list_of_fastaq_to_graph_and_print_coverages_in_all_colours. Exit.\n");
      exit(1);
    }
  

  kmer_window->kmer = (BinaryKmer*) malloc(sizeof(BinaryKmer)*(max_read_length-db_graph->kmer_size-1));
  if (kmer_window->kmer==NULL)
    {
      printf("Failed to malloc kmer_window->kmer in align_list_of_fastaq_to_graph_and_print_coverages_in_all_colours. Exit.\n");
      exit(1);
    }
  kmer_window->nkmers=0;
  
  
  //end of intialisation 
	  
	  
  //create file readers
  int file_reader_fasta(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
    long long ret;
    int offset = 0;
    if (new_entry == false){
      printf("new_entry must be true in hsi test function");
      exit(1);
    }
    ret =  read_sequence_from_fasta(fp,seq,max_read_length,new_entry,full_entry,offset);
    
    return ret;
  }

  int file_reader_fastq(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
    * full_entry = true;

    if (new_entry!= true){
      puts("new_entry has to be true for fastq\n");
      exit(1);
    }

    return read_sequence_from_fastq(fp,seq,max_read_length,fastq_ascii_offset);
  }

  
	  
  dBNode* array_nodes[max_read_length+db_graph->kmer_size+1];
  Orientation array_or[max_read_length+db_graph->kmer_size+1];
  
  //loop through the fasta/q in list_of_fasta/q, and for each, print out a new coverage fasta
  
  FILE* list_fptr = fopen(list_of_fastaq, "r");
  if (list_fptr==NULL)
    {
      printf("Cannot open %s\n", list_of_fastaq);
      exit(1);
    }
  //printf("List of fasta is %s\n", list_of_fastaq);

  char line[MAX_FILENAME_LENGTH+1];
  
  while(fgets(line,MAX_FILENAME_LENGTH, list_fptr) !=NULL)
    {
      
      //remove newline from end of line - replace with \0
      char* p;
      if ((p = strchr(line, '\n')) != NULL)
	{
	  *p = '\0';
	}
      
      //printf("Looking at %s\n", line);

      char outputfile[200];
      sprintf(outputfile,"%s.colour_covgs",line);
      FILE* out = fopen(outputfile, "w");
      if (out ==NULL)
	{
	  printf("Cannot open %s, exiting", outputfile);
	  exit(1);
	}
      //printf("Output to %s\n", outputfile);
      
      FILE* fp = fopen(line, "r");
      if (fp==NULL)
	{
	  printf("Cannot open %s. Exit.\n", line);
	  exit(1);
	}
      
      int dummy_colour_ignored=0;//this will be ignored, as setting to false - we don't want to demand the read all lies in any colour
      int num_kmers;
      do
	{
	  if (format==FASTA)
	    {
	      num_kmers = align_next_read_to_graph_and_return_node_array(fp, max_read_length, array_nodes, array_or, false, file_reader_fasta,
									 seq, kmer_window, db_graph, dummy_colour_ignored);
	    }
	  else if (format==FASTQ)
	    {
	      num_kmers = align_next_read_to_graph_and_return_node_array(fp, max_read_length, array_nodes, array_or, false, file_reader_fastq,
									 seq, kmer_window, db_graph, dummy_colour_ignored);
	    }

	  
	  if (num_kmers>0)
	    {
	      int k;
	      //print out the original read
	      if (is_for_testing==false)
		{
		  fprintf(out, ">%s\n%s\n", seq->name, seq->seq);
		}
	      else//for testing only
		{
		  for_test_array_of_strings[*for_test_index][0]='\0';
		  strcat(for_test_array_of_strings[*for_test_index], seq->name);
		  *for_test_index = *for_test_index+1;

		  for_test_array_of_strings[*for_test_index][0]='\0';
		  strcat(for_test_array_of_strings[*for_test_index], seq->seq);
		  *for_test_index = *for_test_index+1;
		}

	      //print out multiplicities of nodes in each of the colours in turn
	      int j;
	      for (j=0; j<num_of_colours ; j++)//for each colour
		{
		  char read_id[200];
		  sprintf(read_id, "%s_%s_kmer_coverages", seq->name, array_of_names_of_colours[j]);

		  if (is_for_testing==false)//print out read_id
		    {
		      fprintf(out, ">%s\n", read_id);
		    }
		  else
		    {
		      for_test_array_of_strings[*for_test_index][0]='\0';
		      strcat(for_test_array_of_strings[*for_test_index], read_id);
		      *for_test_index = *for_test_index+1;
		    }

		  for (k=0; k<num_kmers; k++)//for each kmer in the read
		    {
		      if (is_for_testing==false)//print covg
			{
			  if (array_nodes[k] !=NULL)
			    {
			      fprintf(out, "%d ", array_nodes[k]->coverage[array_of_colours[j]]);
			    }
			  else
			    {
			      fprintf(out, "0 ");
			    }
			}
		      else
			{
			  for_test_array_of_strings[*for_test_index][0]='\0';
			  char tmp[100];
			  sprintf(tmp, "%d", array_nodes[k]->coverage[array_of_colours[j]]);
			  strcat(for_test_array_of_strings[*for_test_index],  tmp);
			  *for_test_index = *for_test_index+1;
			}
		    }
		  if (is_for_testing==false)
		    {
		      fprintf(out, "\n");
		    }
		  
		}

	    }
	}while((num_kmers>0)||(!feof(fp)) );
      
      
      fclose(out);
      
    }
  fclose(list_fptr);
  
  
}
