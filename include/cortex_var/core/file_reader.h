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

#ifndef FILE_READER_H_
#define FILE_READER_H_

#include <seq.h>
#include <dB_graph.h>
#include <file_format.h>
#include <graph_info.h>

extern int MAX_FILENAME_LENGTH;
extern int MAX_READ_LENGTH;


//pass in bases_read to track amount of sequence read in, and bases_pass_filters_and_loaded to see how much passed filters and got into the graph
void load_se_and_pe_filelists_into_graph_of_specific_person_or_pop(boolean se, boolean pe, char* se_f, char* pe_f1, char* pe_f2,
								   long long* bases_read, long long* bases_pass_filters_and_loaded,long long** readlen_count_array,
								   int qual_thresh, boolean remv_dups_se, int remv_dups_pe, 
								   boolean break_homopolymers, int homopol_limit, int ascii_fq_offset, FileFormat format, 
								   int max_read_length, int colour, dBGraph* db_graph);

void load_fastq_data_from_filename_into_graph_of_specific_person_or_pop(char* filename, long long* bases_read, long long* bases_pass_filters_and_loaded, long long** readlen_count_array,
									     long long * bad_reads,  char quality_cut_off, long long* dup_reads, int max_read_length, 
									     boolean remove_duplicates_single_endedly, boolean break_homopolymers, int homopolymer_cutoff,
									     int fastq_ascii_offset,
									     dBGraph* db_graph, EdgeArrayType type, int index);

void  load_paired_end_data_from_filenames_into_graph_of_specific_person_or_pop(char* filename1, char* filename2, FileFormat format,
									       long long* bases_read, long long* bases_pass_filters_and_loaded,long long** readlen_count_array,
									       long long * bad_reads,  char quality_cut_off, int max_read_length, 
									       long long* dup_reads, boolean remove_duplicates, boolean break_homopolymers, int homopolymer_cutoff, 
									       int fastq_ascii_offset,
									       dBGraph* db_graph, EdgeArrayType type, int index );

void load_fasta_data_from_filename_into_graph_of_specific_person_or_pop(char* filename, long long* bases_read, long long* bases_pass_filters_and_loaded,
									long long * bad_reads, long long* dup_reads, int max_chunk_length, 
									boolean remove_duplicates_single_endedly, boolean break_homopolymers, int homopolymer_cutoff, 
									dBGraph* db_graph, EdgeArrayType type, int index);

//pass in a single kmer sliding window and the Sequence* it was derived from. Will find the nodes correspinding to this seqeunce
//and put them in array. Also will check that edges exist as expected from the Sequence*
void load_kmers_from_sliding_window_into_array(KmerSlidingWindow* kmer_window, Sequence* seq, dBGraph* db_graph, dBNode** array_nodes, Orientation* array_orientations, 
					       int max_array_size, 
					       boolean require_nodes_to_lie_in_given_colour, int colour);

//penultimate argument is to specify whether to discard potential PCR duplicate reads single-endedly - ie if read starts at same kmer
// as a previous read, then discard it. This is a pretty harsh filter, and if ppssible, prefer to use paired end info.
// So in general when calling this function, would expect that boolean remove_dups_single_endedly to be set to false, unless you know you have low coverage, so have
// low probability of two reads starting at the same point.
void load_seq_data_into_graph_of_specific_person_or_pop(FILE* fp, long long* bases_read, long long* bases_pass_filters_and_loaded,long long** readlen_count_array,
							int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length,boolean new_entry, boolean * full_entry), 
							long long * bad_reads, char quality_cut_off, long long* dup_reads, 
							int max_read_length, boolean remove_dups_single_endedly, 
							boolean break_homopolymers, int homopolymer_cutoff, dBGraph * db_graph, EdgeArrayType type, int index);


//assume we have two lists of equal length, of mate files in the same order in each file
void load_list_of_paired_end_files_into_graph_of_specific_person_or_pop(char* list_of_left_mates, char* list_of_right_mates, FileFormat format,
									char quality_cut_off, int max_read_length, 
									long long* bases_read, long long* bases_loaded,long long** readlen_count_array,
									long long* bad_reads, long long* num_dups, int* count_file_pairs, boolean remove_dups, 
									boolean break_homopolymers, int homopolymer_cutoff, int fq_ascii_cutoff, 
									dBGraph* db_graph, EdgeArrayType type, int index);

//index tells you which person within a population it is
//bases_read is passed in to find out how much sequence there was in the files read-in.
//bases_loaded is passed in to find out how much sequence passed filters (qual, PCR dup, homopol) and was loaded into the graph
void load_all_fasta_for_given_person_given_filename_of_file_listing_their_fasta_files(char* f_name, long long* bases_read, long long* bases_loaded,
										      long long* bad_reads, dBGraph* db_graph, int index);


void load_population_as_fasta(char* filename, long long* bases_read, long long* bases_loaded, long long* bad_reads, dBGraph* db_graph);


// gets the next number_of_bases_to_load bases from fasta file, and returns them in the array of nodes.
// assumes this file has already been loaded into the graph.
// returns the number of nodes loaded. If this is less than what you asked for, you know it has hit the end of the file.
// We expect this to be used as follows:
// repeated calls of this function load etc into the LAST number_of_bases_to_load places of the relevant arrays

int load_seq_into_array(FILE* chrom_fptr, int number_of_nodes_to_load, int length_of_arrays, 
			dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels, char* path_string,
			Sequence* seq, KmerSlidingWindow* kmer_window, boolean expecting_new_fasta_entry, dBGraph * db_graph);





//functions for loading multicolour graphs
long long load_multicolour_binary_from_filename_into_graph(char* filename,  dBGraph* db_graph, 
							   int* num_cols_in_loaded_binary, int** array_mean_readlens, long long** array_total_seqs);


long long load_single_colour_binary_data_from_filename_into_graph(char* filename,  dBGraph* db_graph, 
								  int* mean_readlen, long long* total_seq,
								  boolean all_entries_are_unique, EdgeArrayType type, int index,
								  boolean only_load_kmers_already_in_hash, int colour_clean);

long long load_all_binaries_for_given_person_given_filename_of_file_listing_their_binaries(char* filename,  dBGraph* db_graph, GraphInfo* db_graph_info,
											   boolean all_entries_are_unique, EdgeArrayType type, int index,
											   boolean only_load_kmers_already_in_hash, int colour_clean);

long long load_population_as_binaries_from_graph(char* filename, int first_colour,boolean about_to_load_first_binary_into_empty_graph, 
						 dBGraph* db_graph, GraphInfo* db_graph_info,
						 boolean only_load_kmers_already_in_hash, int colour_clean);

void dump_successive_cleaned_binaries(char* filename, int in_colour, int clean_colour, char* suffix, dBGraph* db_graph, GraphInfo* db_graph_info );



//functions for comparing graph with reference, or comparing reads with the graph
void read_ref_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference(FILE* fp, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry),
                                                                            int max_read_length, dBGraph * db_graph);
void read_fastq_and_print_reads_that_lie_in_graph(FILE* fp, FILE* fout, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry),
                                                  long long * bad_reads, int max_read_length, dBGraph * db_graph,
                                                  boolean is_for_testing, char** for_test_array_of_clean_reads, int* for_test_index);
void read_fastq_and_print_subreads_that_lie_in_graph_breaking_at_edges_or_kmers_not_in_graph(FILE* fp, FILE* fout,
											     int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, 
														 boolean * full_entry), 
											     long long * bad_reads, int max_read_length, dBGraph * db_graph, 
											     EdgeArrayType type, int index,
											     boolean is_for_testing, char** for_test_array_of_clean_reads, int* for_test_index);

void read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference(char* f_name, dBGraph* db_graph);
void read_all_ref_chromosomes_and_mark_graph(dBGraph* db_graph);

int get_sliding_windows_from_sequence_breaking_windows_when_sequence_not_in_graph(char * seq,  char * qualities, int length, char quality_cut_off, 
										  KmerSlidingWindowSet * windows, int max_windows, int max_kmers, dBGraph* db_graph);

int get_sliding_windows_from_sequence_requiring_entire_seq_and_edges_to_lie_in_graph(char * seq,  char * qualities, int length, char quality_cut_off, 
										     KmerSlidingWindowSet * windows, int max_windows, int max_kmers, dBGraph* db_graph,
										     EdgeArrayType type, int index); 

void read_fastq_and_print_reads_that_lie_in_graph(FILE* fp, FILE* fout, int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry, boolean * full_entry), 
						  long long * bad_reads, int max_read_length, dBGraph * db_graph,
						  boolean is_for_testing, char** for_test_array_of_clean_reads, int* for_test_index);


void read_fastq_and_print_subreads_that_lie_in_graph_breaking_at_edges_or_kmers_not_in_graph(FILE* fp, FILE* fout,
                                                                                             int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length, boolean new_entry,
                                                                                                                 boolean * full_entry),
                                                                                             long long * bad_reads, int max_read_length, dBGraph * db_graph,
                                                                                             EdgeArrayType type, int index,
                                                                                             boolean is_for_testing, char** for_test_array_of_clean_reads, int* for_test_index);

void read_chromosome_fasta_and_mark_status_of_graph_nodes_as_existing_in_reference(char* f_name, dBGraph* db_graph);
void read_all_ref_chromosomes_and_mark_graph(dBGraph* db_graph);





//gets number_of_bases_to_load's worth of kmers, and returns the corresponding nodes, orientations etc in he array passed in.

int load_seq_into_array(FILE* chrom_fptr, int number_of_bases_to_load, int length_of_arrays,
			    dBNode * * path_nodes, Orientation * path_orientations, Nucleotide * path_labels, char* path_string,
			    Sequence* seq, KmerSlidingWindow* kmer_window, boolean expecting_new_fasta_entry,  dBGraph * db_graph);


int align_next_read_to_graph_and_return_node_array(FILE* fp, int max_read_length, dBNode** array_nodes, Orientation* array_orientations, 
						   boolean require_nodes_to_lie_in_given_colour,
						   int (* file_reader)(FILE * fp, Sequence * seq, int max_read_length,boolean new_entry, boolean * full_entry), 
						   Sequence* seq, KmerSlidingWindow* kmer_window,dBGraph * db_graph, int colour);

int read_next_variant_from_full_flank_file(FILE* fptr, int max_read_length,
                                           dBNode** flank5p,    Orientation* flank5p_or,    int* len_flank5p,
                                           dBNode** ref_allele, Orientation* ref_allele_or, int* len_ref_allele,
                                           dBNode** alt_allele, Orientation* alt_allele_or, int* len_alt_allele,
                                           dBNode** flank3p,    Orientation* flank3p_or,    int* len_flank3p,
                                           dBGraph* db_graph, int colour);


void print_binary_signature(FILE * fp,int kmer_size, int num_cols, int* array_mean_readlens, long long* array_total_seq);
boolean check_binary_signature(FILE * fp,int kmer_size, int bin_version, int* number_of_colours_in_binary, int** array_mean_readlens, long long** array_total_seqs);

boolean query_binary(FILE * fp,int* binary_version, int* kmer_size, int* num_bitfields, int* number_of_colours_in_binary); //return true if binary header readable and has magic number

int get_number_of_files_and_check_existence_from_filelist(char* filelist);
void get_filenames_from_list(char* filelist, char** array, int len);
boolean check_colour_list(char* filename, int kmer);
#endif /* FILE_READER_H_ */
