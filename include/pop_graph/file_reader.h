#ifndef FILE_READER_H_
#define FILE_READER_H_

#include <dB_graph.h>

//for fasta

//new functions
int load_population_as_fasta(char* filename, dBGraph* db_graph);
int load_person_as_fasta_from_filename(char* filename, dBGraph* db_graph, int index);
int load_fasta_data_into_graph_for_specific_person_or_population(FILE* fp, dBGraph * db_graph, EdgeArrayType type, int index);

//unchanged functions
//int load_fasta_data_into_graph(FILE* fp, dBGraph * db_graph);
//int load_fasta_data_from_filename_into_graph(char* filename, dBGraph* db_graph);






//for fastq
//int load_fastq_data_into_graph(FILE* fp, dBGraph * db_graph);
//int load_fastq_data_from_filename_into_graph(char* filename, dBGraph* db_graph);


#endif /* FILE_READER_H_ */
