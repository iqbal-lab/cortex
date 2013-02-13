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
  error_correction.c 
*/

//library headers
#include "string_buffer.h"
#include "seq_file.h"

// cortex_var headers
#include "error_correction.h"
#include "global.h"


inline void error_correct_file_against_graph(char* fastq_file, char quality_cutoff, char ascii_qual_offset,
					     dBGraph *db_graph, char* outfile,
					     int *bases_modified_count_array,//distribution across reads; how many of the read_length bases are fixed
					     int *posn_modified_count_array,//where in the read are we making corrections?
					     int bases_modified_count_array_size,
					     HandleLowQualUncorrectable policy)
{
  quality_cutoff+=ascii_qual_offset;
  short kmer_size = db_graph->kmer_size;
  int count_corrected_bases=0;

  FILE* out_fp = fopen(outfile, "w");
  if (out_fp==NULL)
    {
      die("Unable to open output file %s\n", outfile);
    }

  SeqFile *sf = seq_file_open(fastq_file);
  if(sf == NULL)
    {
      // Error opening file
      fprintf(stderr, "Error: cannot read seq file '%s'\n", fastq_file);
      exit(EXIT_FAILURE);
    }
  char is_fastq = seq_has_quality_scores(sf);
  if (is_fastq==0)
    {
      die("Error correction is only meant to work on FASTQ and this file: %s is not\n", fastq_file);
    }
  
  StrBuf* buf_seq  = strbuf_new();
  StrBuf* buf_qual = strbuf_new();
  StrBuf* working_buf=strbuf_new();

  while(seq_next_read(sf))
    {
      seq_read_all_bases(sf, buf_seq);
      seq_read_all_quals(sf, buf_qual);
      int read_len = seq_get_length(sf);
      int num_kmers = read_len-kmer_size+1;
      int quality_good[read_len];
      set_int_array(quality_good, read_len, 1);

      int first_good=0;//index of first kmer in graph

      //populate the qual array showing which bases have qual >threshold
      //if all quals are high, will Print uncorrected
      //else, if all kmers NOT in graph, will discard or print uncorrected depending on policy
      //else print corrected.
      ReadCorrectionDecison dec = 
	get_first_good_kmer_and_populate_qual_array(buf_seq, buf_qual, num_kmers, read_len,
						    quality_good, quality_cutoff, 
						    &first_good, db_graph, policy);
      

      //*** start of local functions

      //if going right, keep going to right hand end. if going left, keep going to left hand end
      boolean condition(WhichEndOfKmer direction, int pos)
      {
	if ((direction==Right) && (pos<num_kmers))
	  {
	    return true;
	  }
	if ((direction==Left) && (pos>=0))
	  {
	    return true;
	  }
	return false;
      }
      boolean kmer_is_in_graph(char* kmer, dBGraph* db_g)
      {
	BinaryKmer curr_kmer;
	if (seq_to_binary_kmer(kmer, kmer_size, &curr_kmer)==NULL)
	  {
	    //is an N
	    return false;
	  }

	BinaryKmer temp_key;
	element_get_key(&curr_kmer, kmer_size, &temp_key);
	dBNode* node = hash_table_find(&temp_key, db_g);
	if (node==NULL)
	  {
	    return false;
	  }
	else
	  {
	    return true;
	  }
      }
      int increment(int i, WhichEndOfKmer direction)
      {
	if (direction==Right)
	  {
	    return i+1;
	  }
	else
	  {
	    return i-1;
	  }
      }
      char working_str[kmer_size+1];

      // start_pos is in kmer units
      void check_bases_to_end_of_read(int start_pos, ReadCorrectionDecison* decision, WhichEndOfKmer direction)

      {
	if ((start_pos<0) || (start_pos>=num_kmers))
	  {
	    return;
	  }
	int pos=start_pos;
	int offset=0;
	if (direction==Right)
	  {
	    offset= kmer_size-1;
	  }
	char local_kmer[kmer_size+1];
	local_kmer[kmer_size]='\0';	
	BinaryKmer curr_kmer;

	while ( (*decision==PrintCorrected) && (condition(direction,pos)==true) )
	  {
	    strncpy(local_kmer, buf_seq->buff+pos, kmer_size);	  

	    if (quality_good[pos+offset]==1) 
	      {
		//nothing to do
	      }
	    else if (kmer_is_in_graph(local_kmer, db_graph)==true)
	      {
		//nothing to do - don't correct if kmer is in graph
	      }
	    else//kmer not in graph and quality bad
	      {
		boolean fixed = fix_end_if_unambiguous(direction, buf_seq, pos, 
						       working_buf, working_str, db_graph);
		if ( (policy==DiscardReadIfLowQualBaseUnCorrectable) 
		     &&  
		     (fixed==false) )
		  {
		    *decision=Discard;
		  }
		if (fixed==true)
		  {
		    count_corrected_bases++;
		    posn_modified_count_array[offset+pos]++;
		  }
		}
	    pos = increment(pos, direction);
	  }
	
      }			      
      //end of local functions

      if (dec==PrintCorrected)
	{
	  check_bases_to_end_of_read(first_good+1, &dec, Right);
	  check_bases_to_end_of_read(first_good-1, &dec, Left);
	}
      if (dec!=Discard)
	{
	  fprintf(out_fp, ">%s\n%s\n", seq_get_read_name(sf), strbuf_as_str(buf_seq));
	  bases_modified_count_array[count_corrected_bases]++;
	}
    }

    seq_file_close(sf);
    fclose(out_fp);
    strbuf_free(buf_seq);
    strbuf_free(buf_qual);
    strbuf_free(working_buf);
}

//populate two arrays of integers, represent the kmers in the read
//one says whether each kmer is in the graph
//one says whether the final base has quality above the threshold.
//if return value is PrintUncorrected, then kmers_in_graph values may NOT be set
//also first_good_kmer
ReadCorrectionDecison get_first_good_kmer_and_populate_qual_array(StrBuf* seq, StrBuf* qual, 
								  int num_kmers, int read_len,
								  int* quals_good,
								  char quality_cutoff, int* first_good_kmer,
								  dBGraph* dbg, HandleLowQualUncorrectable policy)
{
  int i;
  BinaryKmer curr_kmer;
  BinaryKmer tmp_key;
  dBNode* curr_node;
  char local_kmer[dbg->kmer_size+1];
  local_kmer[dbg->kmer_size]='\0';
  boolean all_qualities_good=true;
  boolean found_first_good_kmer=false;

  for (i=0; i<read_len; i++)
    {
      //check quality good.
      if (strbuf_get_char(qual, i) < quality_cutoff)
	{
	  quals_good[i]=0;
	  all_qualities_good=false;
	}
    }

  if (all_qualities_good==true)
    {
      //performance choice - exit the function immediately,
      //as we know we will print uncorrected.
      //saves us the bother of making hash queries
      return PrintUncorrected;
    }

  for (i=0; (i<num_kmers) && (found_first_good_kmer==false); i++)
    {
      //check kmer in graph
      strncpy(local_kmer, seq->buff+i, dbg->kmer_size);

      if (seq_to_binary_kmer(local_kmer, dbg->kmer_size, &curr_kmer)==NULL)
	{
	  //there is an N
	}
      else
	{
	  element_get_key(&curr_kmer, dbg->kmer_size, &tmp_key);
	  curr_node = hash_table_find(&tmp_key, dbg);
	  if (curr_node==NULL)
	    {
	    }
	  else
	    {
	      found_first_good_kmer=true;
	      *first_good_kmer=i;
	    }
	}
    }

  if (found_first_good_kmer==false)
    {
      //no good kmers, and not all bases are high quality.
      if (policy==DiscardReadIfLowQualBaseUnCorrectable)
	{
	  return Discard;
	}
      else
	{
	  return PrintUncorrected;
	}
    }
  else
    {
      return PrintCorrected;
    }
}
				       

//This function fixes a kmer to match the graph if there is a unique way of doing so
boolean fix_end_if_unambiguous(WhichEndOfKmer which_end, StrBuf* read_buffer, int pos, 
			       StrBuf* kmer_buf, char* kmer_str,//kmer_buf and kmer_str both working variablesm prealloced by caller
			       dBGraph* dbg)
{

  strbuf_substr_prealloced(read_buffer, pos, dbg->kmer_size, kmer_str);
  strbuf_set(kmer_buf, kmer_str);

  BinaryKmer curr_kmer;
  BinaryKmer tmp_key;
  dBNode* curr_node;

  int i;
  int num_ways_of_fixing=0;
  int which_mutant_fixes=-1;
  for (i=0; (i<4) && (num_ways_of_fixing<=1); i++)
    {
      boolean mutation_worth_pursuing=true;
      if (which_end==Left)
	{
	  mutation_worth_pursuing = mutate_base(kmer_buf, 0, i);
	}
      else
	{
	  mutation_worth_pursuing = mutate_base(kmer_buf, dbg->kmer_size -1, i);
	}
      
      if (mutation_worth_pursuing==true)
	{
	  //check if kmer in graph
	  seq_to_binary_kmer(kmer_buf->buff, dbg->kmer_size, &curr_kmer);
	  element_get_key(&curr_kmer, dbg->kmer_size, &tmp_key);
	  curr_node = hash_table_find(&tmp_key, dbg);
	  if (curr_node!=NULL)
	    {
	      num_ways_of_fixing++;
	      which_mutant_fixes=i;	  
	    }
	}
      strbuf_set(kmer_buf, kmer_str);//reset
    }
  if (num_ways_of_fixing==1)
    {
      if (which_end==Left)
	{
	  mutate_base(read_buffer, pos, which_mutant_fixes);
	}
      else
	{
	  mutate_base(read_buffer, pos+dbg->kmer_size -1, which_mutant_fixes);
	}
      return true;
    }
  else
    {
      return false;
    }

}

//if the specified base is ACGT then there are 3 possible mutants
// and which mutant must be between 0 and 2
// else, there are 4 (else means it is N), and which mutant must be between 0 and 3
// returns a boolean saying yes, I made the modification and is worth parsing the modified string
// or no, becase you specified which_mutant=3, with a non-N base.
boolean mutate_base(StrBuf* strbuf, int which_base, int which_mutant)
{
  char mutants[4];
  mutants[3]='N';

  boolean base_is_N=false;
  if (strbuf_get_char(strbuf,which_base)=='A')
    {
      mutants[0]='C';
      mutants[1]='G';
      mutants[2]='T';
    }
  else if (strbuf_get_char(strbuf,which_base)=='C')
    {
      mutants[0]='A';
      mutants[1]='G';
      mutants[2]='T';
    }
  else if (strbuf_get_char(strbuf,which_base)=='G')
    {
      mutants[0]='A';
      mutants[1]='C';
      mutants[2]='T';
    }
  else if (strbuf_get_char(strbuf,which_base)=='T')
    {
      mutants[0]='A';
      mutants[1]='C';
      mutants[2]='G';
    }
  else
    {
      mutants[0]='A';
      mutants[1]='C';
      mutants[2]='G';
      mutants[3]='T';
      base_is_N=true;
    }

  if ( (base_is_N==false) && (which_mutant==3) )
    {
      return false;
    }
  else
    {
      strbuf_set_char(strbuf,which_base,mutants[which_mutant]);
      return true;
    }

}








