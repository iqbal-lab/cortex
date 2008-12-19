include ../Makefile.inc

IDIR_BASIC =../include/basic
IDIR_COUNT =../include/count
IDIR_GRAPH =../include/graph
SRC_COUNT = ./count_src
SRC_GRAPH = ./src_graph
OBJ = ./obj
OBJ_COUNT = ./obj/count
OBJ_GRAPH = ./obj/graph

BIN = ../bin

ifdef MAC
 MACFLAG = -fnested-functions
endif 

CFLAGS_COUNT = -m64 -Wall  -I$(IDIR_BASIC) -I$(IDIR_COUNT) $(MACFLAG)
CFLAGS_GRAPH = -m64 -Wall  -I$(IDIR_BASIC) -I$(IDIR_GRAPH) $(MACFLAG)

GRAPH_OBJ = obj/graph/binary_kmer.o obj/graph/element.o obj/graph/seq.o obj/graph/hash_value.o obj/graph/hash_table.o obj/graph/priority_queue.o obj/graph/dB_graph.o obj/graph/file_reader.o obj/graph/graph.o 

GRAPH_OBJ_CL = obj/graph/binary_kmer.o obj/graph/element.o obj/graph/seq.o obj/graph/hash_value.o obj/graph/hash_table.o obj/graph/priority_queue.o obj/graph/dB_graph.o obj/graph/file_reader.o obj/graph/graph_clipping.o 

COUNT_OBJ = obj/count/binary_kmer.o obj/count/element.o obj/count/seq.o obj/count/hash_value.o obj/count/hash_table.o obj/count/priority_queue.o obj/count/count_kmers.o 


graph : $(GRAPH_OBJ)
	mkdir -p $(BIN); $(CC) $(CFLAGS_GRAPH) $(OPT) $(LDFLAGS) -o $(BIN)/graph $(GRAPH_OBJ)

graph_cl : $(GRAPH_OBJ_CL)
	mkdir -p $(BIN); $(CC) $(CFLAGS_GRAPH) $(OPT) $(LDFLAGS) -o $(BIN)/graph_cp $(GRAPH_OBJ_CL)


count : $(COUNT_OBJ)
	mkdir -p $(BIN); $(CC) $(CFLAGS_COUNT) $(OPT) $(LDFLAGS) -o $(BIN)/count_kmers $(COUNT_OBJ)

obj/graph/%.o :  graph/%.c
	mkdir -p obj/graph;$(CC) $(CFLAGS_GRAPH) $(OPT) $(LDFLAGS) -c $? -o $@

obj/graph/%.o :  basic/%.c
	mkdir -p obj/graph/; $(CC) $(CFLAGS_GRAPH) $(OPT) $(LDFLAGS) -c $? -o $@

obj/count/%.o :  count/%.c
	mkdir -p obj/count;$(CC) $(CFLAGS_COUNT) $(OPT) $(LDFLAGS) -c $? -o $@

obj/count/%.o :  basic/%.c
	mkdir -p obj/count/; $(CC) $(CFLAGS_COUNT) $(OPT) $(LDFLAGS) -c $? -o $@

