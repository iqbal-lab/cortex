/*
 * Copyright 2009-2013 Zamin Iqbal and Mario Caccamo
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
#include <libgen.h> // basename

// cortex_var headers
#include "error_correction.h"
#include "global.h"
#include "dB_graph_population.h"


void error_correct_list_of_files(StrBuf* list_fastq,char quality_cutoff, char ascii_qual_offset,
				 dBGraph *db_graph, HandleLowQualUncorrectable policy,
				 int max_read_len, StrBuf* suffix, char* outdir,
				 boolean add_greedy_bases_for_better_bwt_compression,
				 int num_greedy_bases, boolean rev_comp_read_if_on_reverse_strand)

{

  int len = max_read_len+2;
  uint64_t* distrib_num_bases_corrected     =(uint64_t*) malloc(sizeof(uint64_t)*len);
  uint64_t* distrib_position_bases_corrected=(uint64_t*) malloc(sizeof(uint64_t)*len);
  if ( (distrib_num_bases_corrected==NULL)|| (distrib_position_bases_corrected==NULL))
    {
      die("Unable to alloc arrays for keeping stats. Your machine must have hardly any spare memory\n");
    }
  set_uint64_t_array(distrib_num_bases_corrected,      len, (uint64_t) 0);
  set_uint64_t_array(distrib_position_bases_corrected, len, (uint64_t) 0);

  FILE* list_fastq_fp = fopen(list_fastq->buff, "r");
  if (list_fastq_fp==NULL)
    {
      printf("Cannot open file %s\n", list_fastq->buff);
    }
  StrBuf *next_fastq     = strbuf_new();
  StrBuf* corrected_file = strbuf_new();
  StrBuf* corrected_file_newpath = strbuf_new();

  while(strbuf_reset_readline(next_fastq, list_fastq_fp))
    {
      strbuf_chomp(next_fastq);
      if(strbuf_len(next_fastq) > 0)
	 {
	   strbuf_reset(corrected_file);
	   strbuf_reset(corrected_file_newpath);
	   strbuf_copy(corrected_file, 0,//dest
		       next_fastq,0,strbuf_len(next_fastq));
	   strbuf_append_str(corrected_file, suffix->buff);
	   char* corrected_file_basename = basename(corrected_file->buff);
	   strbuf_append_str(corrected_file_newpath, outdir);
	   strbuf_append_str(corrected_file_newpath, corrected_file_basename);

	   error_correct_file_against_graph(next_fastq->buff, quality_cutoff, ascii_qual_offset,
					    db_graph, corrected_file_newpath->buff,
					    distrib_num_bases_corrected,
					    distrib_position_bases_corrected,
					    len,
					    policy,
					    add_greedy_bases_for_better_bwt_compression,
					    num_greedy_bases, rev_comp_read_if_on_reverse_strand);
	     

	 }
    }
  fclose(list_fastq_fp);
  strbuf_free(next_fastq);
  strbuf_free(corrected_file);
  strbuf_free(corrected_file_newpath);
  free(distrib_num_bases_corrected);
  free(distrib_position_bases_corrected);
}

//outputs fastQ unless add_greedy_bases_for_better_bwt_compression==true, in which case is for 1000genomes, and they want fastA
inline void error_correct_file_against_graph(char* fastq_file, char quality_cutoff, char ascii_qual_offset,
					     dBGraph *db_graph, char* outfile,
					     uint64_t *bases_modified_count_array,//distribution across reads; how many of the read_length bases are fixed
					     uint64_t *posn_modified_count_array,//where in the read are we making corrections?
					     int bases_modified_count_array_size,
					     HandleLowQualUncorrectable policy,
					     boolean add_greedy_bases_for_better_bwt_compression,
					     int num_greedy_bases,
					     boolean rev_comp_read_if_on_reverse_strand)
{
  //reset the stats arrays, we get stats per input file
  set_uint64_t_array(bases_modified_count_array,bases_modified_count_array_size, (uint64_t) 0);
  set_uint64_t_array(posn_modified_count_array, bases_modified_count_array_size, (uint64_t) 0);


  //set some variables, quality etc
  quality_cutoff+=ascii_qual_offset;
  short kmer_size = db_graph->kmer_size;

  //setup output file for corrected reads, plus two stats files
  FILE* out_fp = fopen(outfile, "w");
  if (out_fp==NULL)
    {
      die("Unable to open output file %s\n", outfile);
    }
  char* suff1 = ".distrib_num_modified_bases";
  char* suff2 = ".distrib_posn_modified_bases";
  char* suff3 =".read_stats";
  char* stat1 = (char*) malloc(sizeof(char)*(strlen(outfile)+strlen(suff1)+1));
  char* stat2 = (char*) malloc(sizeof(char)*(strlen(outfile)+strlen(suff2)+1));
  char* stat3 = (char*) malloc(sizeof(char)*(strlen(outfile)+strlen(suff3)+1));

  if ( (stat1==NULL) || (stat2==NULL) || (stat3==NULL) )
    {
      die("Unable to malloc FILENAME strings. Something badly wrong with your server\n");
    }
  set_string_to_null(stat1, strlen(outfile)+strlen(suff1)+1);
  set_string_to_null(stat2, strlen(outfile)+strlen(suff2)+1);
  set_string_to_null(stat3, strlen(outfile)+strlen(suff3)+1);
  strcpy(stat1, outfile);
  strcat(stat1, suff1);
  strcat(stat2, outfile);
  strcat(stat2, suff2);
  strcat(stat3, outfile);
  strcat(stat3, suff3);

  FILE* out_stat1 = fopen(stat1, "w");
  FILE* out_stat2 = fopen(stat2, "w");
  FILE* out_stat3 = fopen(stat3, "w");
  if ( (out_stat1==NULL)|| (out_stat2==NULL) || (out_stat3==NULL) )
    {
      die("Unable to open %s or %s or %s to write to - permissions issue?\n", stat1, stat2, stat3);
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
  dBNode* last_node_in_read;
  Orientation last_or_in_read;
  int num_original_reads=0, num_final_reads=0, num_corrected_reads=0, num_discarded_reads=0;
  while(seq_next_read(sf))
    {
      int count_corrected_bases=0;
      //NOTE - uses modified version fo Isaacs code - new func
      seq_read_all_bases_and_quals(sf, buf_seq, buf_qual);
      StrBuf* buf_seq_debug  = strbuf_clone(buf_seq);
      
      int read_len = seq_get_length(sf);
      int num_kmers = read_len-kmer_size+1;
      int quality_good[read_len];
      set_int_array(quality_good, read_len, 1);

      int first_good=0;//index of first kmer in graph

      //populate the qual array showing which bases have qual >threshold
      //if all quals are high, will Print uncorrected
      //else, if all kmers NOT in graph, will discard or print uncorrected depending on policy
      //else print corrected.
      Orientation strand_first_good_kmer;
      ReadCorrectionDecison dec = 
	get_first_good_kmer_and_populate_qual_array(seq_get_read_name(sf), buf_seq, buf_qual, num_kmers, read_len,
						    quality_good, quality_cutoff, 
						    &first_good, &strand_first_good_kmer,
						    db_graph, policy, rev_comp_read_if_on_reverse_strand);
      

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
      boolean check_bases_to_end_of_read(int start_pos, ReadCorrectionDecison* decision, 
					 WhichEndOfKmer direction, 
					 int* num_corrected_bases_in_this_read_debug)

      {
	boolean any_correction_done=false;
	if ((start_pos<0) || (start_pos>=num_kmers))
	  {
	    return any_correction_done;
	  }
	int pos=start_pos;
	int offset=0;
	if (direction==Right)
	  {
	    offset= kmer_size-1;
	  }
	char local_kmer[kmer_size+1];
	local_kmer[kmer_size]='\0';	

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
		boolean fixed = fix_end_if_unambiguous(direction, buf_seq, buf_qual, quality_cutoff, pos, 
						       working_buf, working_str, db_graph);
		if ( (policy==DiscardReadIfLowQualBaseUnCorrectable) 
		     &&  
		     (fixed==false) )
		  {
		    *decision=Discard;
		  }
		else if (fixed==true)
		  {
		    any_correction_done=true;
		    count_corrected_bases++;
		    *num_corrected_bases_in_this_read_debug=*num_corrected_bases_in_this_read_debug+1;
		    if (offset+pos<bases_modified_count_array_size)
		      {
			posn_modified_count_array[offset+pos]++;
		      }
		    else
		      {
			posn_modified_count_array[bases_modified_count_array_size-1]++;
		      }
		  }
	      }
	    pos = increment(pos, direction);
	  }
	return any_correction_done;
      }			      
      //end of local functions


      num_original_reads++;

      boolean any_fixing_done=false;//remember you can do some fixing but then decide later to discard

      if (dec==PrintCorrected)//this means will try and correct, but might not be able to
	{
	  int num_corrections_to_this_read_for_debug=0;
	  boolean any_fix_right = check_bases_to_end_of_read(first_good+1, &dec, Right, &num_corrections_to_this_read_for_debug);
	  boolean any_fix_left  = check_bases_to_end_of_read(first_good-1, &dec, Left, &num_corrections_to_this_read_for_debug);
	  if (any_fix_right||any_fix_left)
	    {
	      //then it really was able to correct at least one base somewhere
	      any_fixing_done=true;
	    }
	}

      if (dec==Discard)
	{
	  num_discarded_reads++;
	}
      else
	{
	  if (add_greedy_bases_for_better_bwt_compression==false)
	    {
	      fprintf(out_fp, "@");
	    }
	  else
	    {
	      fprintf(out_fp, ">");
	    }
	  fprintf(out_fp, "%s\n", seq_get_read_name(sf));
	  char last_kmer_str[db_graph->kmer_size];
	  if ( (rev_comp_read_if_on_reverse_strand==true)
	       && (strand_first_good_kmer==reverse) )
	    {
	      strbuf_rev_comp(buf_seq);
	      strbuf_reverse(buf_qual);
	    }
	  //get the final kmer (after reverse complementing if that is happening)
	  strbuf_substr_prealloced(buf_seq, buf_seq->len-db_graph->kmer_size-1, db_graph->kmer_size, last_kmer_str);

	  if (add_greedy_bases_for_better_bwt_compression==true)
	    {
	      BinaryKmer last_kmer;
	      if (seq_to_binary_kmer(last_kmer_str, kmer_size, &last_kmer)==NULL)
		{
		  //contains an N
		  pad_to_N_with_Adenine(buf_seq, buf_seq->len+num_greedy_bases);
		}
	
	      BinaryKmer temp_key;
	      element_get_key(&last_kmer, kmer_size, &temp_key);
	      last_node_in_read = hash_table_find(&temp_key, db_graph);
	      if (last_node_in_read==NULL)
		{
		  pad_to_N_with_Adenine(buf_seq, buf_seq->len+num_greedy_bases);
		}
	      else
		{
		  last_or_in_read = db_node_get_orientation(&last_kmer,last_node_in_read, db_graph->kmer_size);
		  take_n_greedy_random_steps(last_node_in_read, last_or_in_read,
					     db_graph, num_greedy_bases, buf_seq);
		}
	    }
	  
	  fprintf(out_fp, "%s\n", buf_seq->buff );
	  if (add_greedy_bases_for_better_bwt_compression==false)
	    {
	      fprintf(out_fp, "+\n%s\n", buf_qual->buff);
	    }
	  if (count_corrected_bases<bases_modified_count_array_size)
	    {
	      bases_modified_count_array[count_corrected_bases]++;
	    }
	  else
	    {
	      bases_modified_count_array[bases_modified_count_array_size-1]++;
	    }
	  num_final_reads++;
	  if (any_fixing_done==true)
	    {
	      num_corrected_reads++;
	    }
	}
      strbuf_free(buf_seq_debug);
    }

    seq_file_close(sf);
    fclose(out_fp);

    //write out the stats files
    int i;
    for (i=0; i<bases_modified_count_array_size; i++)
      {
	fprintf(out_stat1, "%d\t%" PRIu64 "\n", i, bases_modified_count_array[i]);
	fprintf(out_stat2, "%d\t%" PRIu64 "\n", i, posn_modified_count_array[i]);
      }
    fclose(out_stat1);
    fclose(out_stat2);
    fprintf(out_stat3, "Original reads:\t%d\n", num_original_reads);
    fprintf(out_stat3, "Final reads:\t%d\n", num_final_reads);
    fprintf(out_stat3, "Corrected reads:\t%d\n", num_corrected_reads);
    fprintf(out_stat3, "Discarded reads:\t%d\n", num_discarded_reads);
    fclose(out_stat3);
    strbuf_free(buf_seq);
    strbuf_free(buf_qual);
    strbuf_free(working_buf);
    free(stat1);
    free(stat2);
    free(stat3);

}

//populate two arrays of integers, represent the kmers in the read
//one says whether each kmer is in the graph
//one says whether the final base has quality above the threshold.
//if return value is PrintUncorrected, then kmers_in_graph values may NOT be set
//also first_good_kmer
ReadCorrectionDecison get_first_good_kmer_and_populate_qual_array(const char* debug_read_id, StrBuf* seq, StrBuf* qual, 
								  int num_kmers, int read_len,
								  int* quals_good,
								  char quality_cutoff, int* first_good_kmer,
								  Orientation* strand_first_good_kmer,
								  dBGraph* dbg, HandleLowQualUncorrectable policy,
								  boolean we_will_want_to_revcomp_reads_on_rev_strand)
{
  *strand_first_good_kmer=forward;
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

      //debug:
      if (strlen(qual->buff)==0)
	{
	  printf("This read: id %s, seq: %s has empty qual: %s\n", debug_read_id, seq->buff, qual->buff);
	  die("Debug this!");
	}
      if (strbuf_get_char(qual, i) < quality_cutoff)
	{
	  quals_good[i]=0;
	  all_qualities_good=false;
	}
    }

 
  if ( (we_will_want_to_revcomp_reads_on_rev_strand==false) &&(all_qualities_good==true))
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
	      Orientation or_first_good_kmer  = db_node_get_orientation(&curr_kmer,curr_node, dbg->kmer_size); 
	      if ( (db_node_check_status(curr_node, rv_strand)==true)
		   && (or_first_good_kmer==forward) )
		{
		  *strand_first_good_kmer = reverse;
		}
	      else if ( (db_node_check_status(curr_node, fw_strand)==true)
		   && (or_first_good_kmer==reverse) )
		{
		  *strand_first_good_kmer = reverse;
		}
	      if (all_qualities_good==true)
		{
		  return PrintUncorrected;
		}
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
//the quality scores are here ONLY to allow fixing of the quality score at any fixed base.
boolean fix_end_if_unambiguous(WhichEndOfKmer which_end, StrBuf* read_buffer, StrBuf* qual_buffer, char quality_cutoff, int pos, 
			       StrBuf* kmer_buf, char* kmer_str,//kmer_buf and kmer_str both working variablesm prealloced by caller
			       dBGraph* dbg)
{
  //strncpy(kmer_str, read_buffer->buff+pos, dbg->kmer_size);
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
	  set_qual_to_just_above_cutoff(qual_buffer, pos, quality_cutoff);
	}
      else
	{
	  mutate_base(read_buffer, pos+dbg->kmer_size -1, which_mutant_fixes);
	  set_qual_to_just_above_cutoff(qual_buffer, pos+dbg->kmer_size-1, quality_cutoff);
	}
      return true;
    }
  else
    {
      return false;
    }

}

//if you do modify a base, you should also fix it's quality score.
//for the moment, do something ugly and fix it to be 1 more than the quality cutoff
void set_qual_to_just_above_cutoff(StrBuf* qualbuf, int pos, char cutoff)
{
  strbuf_set_char(qualbuf,pos,cutoff+1);
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


void read_ref_fasta_and_mark_strand(char* ref_fa, dBGraph * db_graph)
{

  int max_read_length = 20000;//chunk length for reading ref.

  FILE* fp = fopen(ref_fa, "r");
  if (fp==NULL)
    {
      printf("Cannot open reference fasta %s\n", ref_fa);
      die("Exit");
    }

  //file reader
  int file_reader(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry){
    long long ret;
    int offset = 0;
    if (new_entry == false){
      offset = db_graph->kmer_size;
      //die("new_entry must be true in hsi test function");
    }
    ret =  read_sequence_from_fasta(fp,seq,max_read_length,new_entry,full_entry,offset);
    
    return ret;
  }


    //----------------------------------
  // preallocate the memory used to read the sequences
  //----------------------------------
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  alloc_sequence(seq,max_read_length,LINE_MAX);
  
  
  int seq_length=0;
  short kmer_size = db_graph->kmer_size;

  //max_read_length/(kmer_size+1) is the worst case for the number of sliding windows, ie a kmer follow by a low-quality/bad base
  int max_windows = max_read_length/(kmer_size+1);
 
  //number of possible kmers in a 'perfect' read
  int max_kmers   = max_read_length-kmer_size+1;

  //we will be calling functions that apply to fastq as well as fasta. Since we are loading a ref fasta, we do not want to cut homopolymers,
  // and since we have fastq, the qualit-cutfoff os redundant
  boolean break_homopolymers=false;
  int homopolymer_cutoff=0;
  char quality_cut_off=0;
  

  //----------------------------------
  //preallocate the space of memory used to keep the sliding_windows. NB: this space of memory is reused for every call -- with the view 
  //to avoid memory fragmentation
  //NB: this space needs to preallocate memory for orthogonal situations: 
  //    * a good read -> few windows, many kmers per window
  //    * a bad read  -> many windows, few kmers per window    
  //----------------------------------
  KmerSlidingWindowSet * windows = malloc(sizeof(KmerSlidingWindowSet));  
  if (windows == NULL){
    fputs("Out of memory trying to allocate a KmerArraySet",stderr);
    exit(1);
  }  
  //allocate memory for the sliding windows 
  binary_kmer_alloc_kmers_set(windows, max_windows, max_kmers);

  BinaryKmer tmp_kmer;

  boolean full_entry = true;
  boolean prev_full_entry = true;
  
  int entry_length;

  while ((entry_length = file_reader(fp,seq,max_read_length, full_entry, &full_entry))){

    if (DEBUG){
      printf ("\nsequence %s\n",seq->seq);
    }
    
    int i,j;
    seq_length += (long long) (entry_length - (prev_full_entry==false ? db_graph->kmer_size : 0));

    
    //last two arguments mean not to break at homopolymers
    int nkmers = get_sliding_windows_from_sequence(seq->seq,seq->qual,entry_length,quality_cut_off,db_graph->kmer_size,windows,max_windows, 
						   max_kmers,break_homopolymers, homopolymer_cutoff);
    
    if (nkmers == 0) {
    }
    else {
      Element * current_node  = NULL;
      Orientation current_orientation;
      
      for(i=0;i<windows->nwindows;i++){ //for each window
	KmerSlidingWindow * current_window = &(windows->window[i]);
	
	for(j=0;j<current_window->nkmers;j++){ //for each kmer in window
	  
	  current_node = hash_table_find(element_get_key(&(current_window->kmer[j]),db_graph->kmer_size, &tmp_kmer),db_graph);	  	 

	  if (current_node){
	    current_orientation = db_node_get_orientation(&(current_window->kmer[j]),current_node, db_graph->kmer_size);
	    if (current_orientation==forward)
	      {
		db_node_set_status(current_node, fw_strand);
	      }
	    else
	      {
		db_node_set_status(current_node, rv_strand);
	      }
	  }
	  
	}
	
      }
    }
    
    if (full_entry == false){
      shift_last_kmer_to_start_of_sequence(seq,entry_length,db_graph->kmer_size);
    }

    prev_full_entry = full_entry;

  }
  
  free_sequence(&seq);
  binary_kmer_free_kmers_set(&windows);

}


void pick_a_random_edge(dBNode* node, Orientation or, Nucleotide* random_nuc)
{
  Edges edges = get_union_of_edges(*node);
  Edges edges_saved = edges;
  short edges_count = 0;

  if (or == reverse){
    edges >>= 4;
    edges_saved >>=4;
  }
  int n;

  int which_edges_present[4]={0,0,0,0};

  for(n=0;n<4;n++){
    
    if ((edges & 1) == 1){
      edges_count++;
      which_edges_present[n]=1;
    }
    edges >>= 1;    
  }

  if (edges_count==0)
    {
      *random_nuc=Undefined;
    }
  else
    {
      boolean found=false;
      while (found==false)
	{

	  //not that it matters, but see http://stackoverflow.com/questions/2509679/how-to-generate-a-random-number-from-within-a-range-c
	  int which_edge = rand() % 4;
	  if (which_edges_present[which_edge]==1)
	    {
	      *random_nuc =  (Nucleotide) which_edge;
	      found=true;
	    }
	}
    }
}

void pad_to_N_with_Adenine(StrBuf* strbuf, int n)
{
  int i;
  for (i=strbuf->len; i<n; i++)
    {
      strbuf_append_char(strbuf, 'A');
    }
}


void take_n_greedy_random_steps(dBNode* start_node, Orientation start_or, dBGraph* db_graph,
				int num_steps, StrBuf* greedy_seq)

{
  
  Orientation  current_orientation,next_orientation;
  dBNode * current_node = NULL;
  dBNode * next_node = NULL;
  Nucleotide nucleotide,rev_nucleotide,nucleotide2;
  int length =0;
  char tmp_seq[db_graph->kmer_size+1];
  tmp_seq[db_graph->kmer_size]='\0';

  //sanity checks
  if (start_node == NULL)
    {
      pad_to_N_with_Adenine(greedy_seq, num_steps);
      return;
    }

  current_node        = start_node;
  current_orientation = start_or;

  //first edge defined
  pick_a_random_edge(start_node, start_or, &nucleotide);

  if (nucleotide==Undefined)
    {
      pad_to_N_with_Adenine(greedy_seq, num_steps);
      return;
    }
  
  do{ 
    length++;
    next_node =  db_graph_get_next_node_in_subgraph_defined_by_func_of_colours(current_node,current_orientation,
									       &next_orientation,nucleotide,
									       &rev_nucleotide,
									       db_graph, 
									       &element_get_colour_union_of_all_colours);
    
    
    strbuf_append_char(greedy_seq, binary_nucleotide_to_char(nucleotide));
    current_node        = next_node;
    current_orientation = next_orientation;
    
    pick_a_random_edge(current_node, current_orientation, &nucleotide);
    if (nucleotide==Undefined)
      {
	pad_to_N_with_Adenine(greedy_seq, num_steps);
	return;
      }
    
  } while (length<num_steps);
  
}






