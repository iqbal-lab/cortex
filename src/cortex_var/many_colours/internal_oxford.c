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
  internal_oxford.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// cortex_var headers
#include "element.h"
#include "file_reader.h"
#include "dB_graph.h"
#include "dB_graph_population.h"


void set_ref_chromosome_file_pointers(char** reference_chromosome_file_ptrs, int num_chromosomes)
{
  //  if ( (num_chromosomes != 25) && (num_chromosomes != 24) )
  //  {
  //    die("Expected 24 or 25 ref chromosomes. Exit");
  //  }
  

  reference_chromosome_file_ptrs[0] = "/data/zam/ref/hs/fasta/Homo_sapiens.NCBI36.52.dna.chromosome.MT.fa";

  int i;
  for (i=1; i<23; i++)
    {
      sprintf(reference_chromosome_file_ptrs[i], "/data/zam/ref/hs/fasta/Homo_sapiens.NCBI36.52.dna.chromosome.%i.fa", i);
    }
  reference_chromosome_file_ptrs[23]="/data/zam/ref/hs/fasta/Homo_sapiens.NCBI36.52.dna.chromosome.X.fa";

  //if (num_chromosomes==24)
  //  {
  //    reference_chromosome_file_ptrs[24]="/data/zam/ref/hs/fasta/Homo_sapiens.NCBI36.52.dna.chromosome.Y.fa";
  //  }

}

void create_uniqueness_file(char** ref_chroms, dBGraph* db_graph)
{
  printf("For each reference chromosome, print a fasta file which has a 1 if the 31-mer starting at that base exists precisely once in the reference, and 0 otherwise\n");
	
  char** uniq_files = malloc( sizeof(char*) * 25); //one for each chromosome, ignoring haplotypes like MHC
  if (uniq_files==NULL)
    {
      die("Out of memory. Give up can't even allocate space for the names of the uniq files");
    }
  int i;
  for (i=0; i< 25; i++)
    {
      uniq_files[i] = malloc(sizeof(char)*100); 
      if (uniq_files[i]==NULL)
	{
	  die("Out of memory. Giveup can't even allocate space for the names of the uniqueness files for  i = %d",i);
	}
    }
  
  uniq_files[0] = "kmer_uniqueness_in_chrom_MT";
  for (i=1; i<23; i++)
    {
      sprintf(uniq_files[i], "kmer_uniqueness_in_chrom_%i", i);
    }
  uniq_files[23]="kmer_uniqueness_in_chrom_X";
  uniq_files[24]="kmer_uniqueness_in_chrom_Y";
  
  
  for (i=1; i<25; i++) 
    {
      printf("Print kmer uniqueness-in-entire-ref file for  %s\n", ref_chroms[i]);
      
      FILE* chrom_fptr = fopen(ref_chroms[i], "r");
      if (chrom_fptr==NULL)
	{
	  die("Cannot open %s", ref_chroms[i]);
	}
      
      FILE* out_fptr = fopen(uniq_files[i], "w");
      if (out_fptr==NULL)
	{
	  die("Cannot open %s for output", uniq_files[i]);
	}
      
      BinaryKmer marked_kmer; //will have all longlongs in array being ~0
      int i;
      for (i=0; i<NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
	{
	  marked_kmer[i]=~0;
	}
      
      void print_uniqueness(dBNode* node)
      {
	if ((node==NULL) || (binary_kmer_comparison_operator(node->kmer, marked_kmer)) )
	  {
	    fprintf(out_fptr, "0");
	  }
	else
	  {
	    // we assume that the reference is individual at index 0.
	    if (db_node_get_coverage(node,0)==1)
	      {
		//this kmer exists precisely once in the ref
		fprintf(out_fptr, "1");
	      }
	    else
	      {
		fprintf(out_fptr,"0");
	      }
	  }
      }
      
      apply_to_all_nodes_in_path_defined_by_fasta(&print_uniqueness, chrom_fptr, 10000, db_graph);//work through fasta in chunks of size 10kb
      fclose(out_fptr);
      fclose(chrom_fptr);
      
      
    }
  
}



void print_coverage_and_ref_multiplicities_for_list_of_fasta(char* list_of_fasta, dBGraph* db_graph)
{
  	printf("Take the file %s, and for each fasta file listed within, print out its own output file of coverage/number_of_occurrences_in_reference for each edge, tab separated. We expect no N's in the fasta", list_of_fasta);

	FILE* fptr = fopen(list_of_fasta, "r");
	if (fptr==NULL)
	  {
	    die("Cannot open %s", list_of_fasta);
	  }

	//file contains a list of fasta file names
	char line[MAX_FILENAME_LENGTH+1];
	
	while(fgets(line,MAX_FILENAME_LENGTH, fptr) !=NULL)
	  {
	    
	    //remove newline from endof line- replace with \0
	    char* p;
	    if ((p = strchr(line, '\n')) != NULL)
	      *p = '\0';
	    
	    //open fasta file
	    FILE* path_fptr = fopen(line, "r");
	    if (path_fptr==NULL)
	      {
		die("Cannot open %s", line);
	      }

	    //create output file 
	    char name[300];
	    sprintf(name,"%s.covg",line);

	    FILE* out_fptr = fopen(name, "w");
	    if (out_fptr==NULL)
	      {
		die("Cannot open %s", name);
	      }


	    BinaryKmer marked_kmer; //will have all longlongs in array being ~0
	    int i;
	    for (i=0; i<NUMBER_OF_BITFIELDS_IN_BINARY_KMER; i++)
	      {
		marked_kmer[i]=~0;
	      }

	    
	    void print_covg_and_ref_covg(dBNode* node)
	    {
		if (node==NULL) 
		  {
		    fprintf(out_fptr, "node is NULL\t");
		  }
		else if ( binary_kmer_comparison_operator(node->kmer, marked_kmer)==true ) 
		  {
		    fprintf(out_fptr, "N in fasta - unexpected\t");
		  }
		else
		  {
		      // we assume that the reference is individual at index 0.
		    int ref_covg =  db_node_get_coverage(node,0);
		    int covg     =  db_node_get_coverage(node,1);
		    fprintf(out_fptr, "%d/%d\t", covg, ref_covg);
		  }
	      }
	    
	    apply_to_all_nodes_in_path_defined_by_fasta(&print_covg_and_ref_covg, path_fptr, 2, db_graph);//work through fasta in chunks of size 50 bases - fasta not too big
	    fclose(out_fptr);
	    fclose(path_fptr);
	    
	    
	  }
	
	fclose(fptr);

}
