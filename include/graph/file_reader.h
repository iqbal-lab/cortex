#ifndef FILE_READER_H_
#define FILE_READER_H_

#include <dB_graph.h>

int load_fasta_data_into_graph(FILE* fp, dBGraph * db_graph);
int load_fasta_data_from_filename_into_graph(char* filename, dBGraph* db_graph);

#endif /* FILE_READER_H_ */
