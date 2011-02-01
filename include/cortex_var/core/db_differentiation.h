#include <element.h>
#include <stdio.h>
#include <stdlib.h>
#include <dB_graph.h>
#include <seq.h>
#include <limits.h>
#include <file_reader.h>
#include <db_variants.h>


void align_list_of_fastaq_to_graph_and_print_coverages_in_all_colours(FileFormat format, char* list_of_fastaq, int max_read_length, int* array_of_colours, char** array_of_names_of_colours,
								      int num_of_colours, dBGraph* db_graph,int fastq_ascii_offset,
								      boolean is_for_testing, char** for_test_array_of_strings, int* for_test_index);
