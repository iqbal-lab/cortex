CC = gcc

IDIR_TEST  =../include/test
IDIR_BASIC =../include/basic
IDIR_HASH  =../include/basic/hash_table 
IDIR_GRAPH =../include/graph

IDIR_POP_GRAPH = ../include/pop_graph
IDIR_CUNIT = /homes/zi/dev/hg/CUnit/CUnit-2.1-0/CUnit/Headers

BIN = ../bin


ifdef MAC
 MACFLAG    = -fnested-functions
 IDIR_CUNIT = /Users/zi/dev/hg/laptop_local/repos/CUnit/CUnit-2.1-0/CUnit/Headers
endif 

OPT  		      = -m64 -Wall -O0 -g
CFLAGS_GRAPH          = -I$(IDIR_BASIC) -I$(IDIR_HASH) -I$(IDIR_GRAPH) $(MACFLAG) 
CFLAGS_TEST	      = -I$(IDIR_TEST)  -I$(IDIR_CUNIT) 
CFLAGS_CUNIT          = -L/homes/zi/bin/lib

GRAPH_OBJ_DUMP = obj/basic/binary_kmer.o obj/graph/element.o obj/basic/seq.o obj/basic/hash_table/hash_key/bob_jenkins/hash_value.o obj/basic/hash_table/open_hash/hash_table.o obj/graph/dB_graph.o obj/graph/file_reader.o obj/graph/graph_dump.o 

GRAPH_OBJ      = obj/basic/binary_kmer.o obj/graph/element.o obj/basic/seq.o obj/basic/hash_table/hash_key/bob_jenkins/hash_value.o obj/basic/hash_table/open_hash/hash_table.o obj/graph/dB_graph.o obj/graph/file_reader.o obj/graph/graph.o 

BASIC_TESTS_OBJ = obj/basic/binary_kmer.o obj/basic/seq.o obj/test/test_binary_kmer.o obj/test/test_seq.o obj/test/run_basic_tests.o

DB_GRAPH_TESTS_OBJ = obj/basic/binary_kmer.o obj/graph/element.o obj/basic/seq.o obj/basic/hash_table/hash_key/bob_jenkins/hash_value.o obj/basic/hash_table/open_hash/hash_table.o obj/graph/dB_graph.o obj/graph/file_reader.o obj/test/test_dB_graph.o  obj/test/test_file_reader.o obj/test/run_dB_graph_tests.o


graph_dump : $(GRAPH_OBJ_DUMP)
	mkdir -p $(BIN); $(CC) $(OPT) $(LDFLAGS) -o $(BIN)/graph_dump $(GRAPH_OBJ_DUMP)

graph : $(GRAPH_OBJ)
	mkdir -p $(BIN); $(CC) $(OPT) $(LDFLAGS) -o $(BIN)/graph $(GRAPH_OBJ)

run_basic_tests : $(BASIC_TESTS_OBJ)
	mkdir -p $(BIN); $(CC) $(OPT) $(CFLAGS_CUNIT) $(LDFLAGS) -o $(BIN)/run_basic_tests $(BASIC_TESTS_OBJ) -lcunit

run_dB_graph_tests : $(DB_GRAPH_TESTS_OBJ)
	mkdir -p $(BIN); $(CC) $(OPT) $(CFLAGS_CUNIT) $(LDFLAGS) -o $(BIN)/run_dB_graph_tests $(DB_GRAPH_TESTS_OBJ) -lcunit 


#pattern rules

obj/basic/hash_table/open_hash/%.o :  basic/hash_table/open_hash/%.c 
	mkdir -p obj/basic/hash_table/open_hash; $(CC) $(CFLAGS_GRAPH) $(OPT) $(LDFLAGS) -c $? -o $@

obj/basic/hash_table/hash_key/bob_jenkins/%.o :  basic/hash_table/hash_key/bob_jenkins/%.c
	mkdir -p obj/basic/hash_table/hash_key/bob_jenkins; $(CC) $(CFLAGS_GRAPH) $(OPT) $(LDFLAGS) -c $? -o $@

obj/basic/hash_table/hash_key/crc32/%.o :  basic/hash_table/hash_key/crc32/%.c
	mkdir -p obj/basic/hash_table/hash_key/crc32; $(CC) $(CFLAGS_GRAPH) $(OPT) $(LDFLAGS) -c $? -o $@

obj/basic/%.o : basic/%.c 
	mkdir -p obj/basic/; $(CC) $(CFLAGS_GRAPH) $(OPT) $(LDFLAGS) -c $? -o $@

obj/graph/%.o :  graph/%.c
	mkdir -p obj/graph/; $(CC) $(CFLAGS_GRAPH) $(OPT) $(LDFLAGS) -c $? -o $@

obj/test/%.o :  test/%.c 
	mkdir -p obj/test; $(CC) $(CFLAGS_GRAPH) $(CFLAGS_TEST) $(OPT) $(LDFLAGS) -c $? -o $@

