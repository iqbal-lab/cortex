IDIR_BASIC =../include/basic
IDIR_GRAPH =../include/graph
IDIR_POP_GRAPH =../include/pop_graph

IDIR_CUNIT =/homes/zi/dev/hg/CUnit/CUnit-2.1-0/CUnit/Headers

CC = gcc
CFLAGS = -Wall -I$(IDIR_BASIC) -I$(IDIR_GRAPH) -I$(IDIR_POP_GRAPH) -I$(IDIR_CUNIT) -g -O0
CFLAGS_GRAPH = -Wall -I$(IDIR_BASIC) -I$(IDIR_GRAPH)  -I$(IDIR_CUNIT) -g -O0
CFLAGS_POP_GRAPH = -Wall -I$(IDIR_BASIC) -I$(IDIR_POP_GRAPH) -I$(IDIR_CUNIT) -g -O0

TEST_OBJ_DIR=obj
POP_GRAPH_TEST_OBJ_DIR=obj/pop_graph_test
GRAPH_TEST_OBJ_DIR=obj/graph_test

GRAPH_OBJ_DIR=../src/obj/graph
POP_GRAPH_OBJ_DIR=../src/obj/pop_graph



GRAPH_TEST_OBJ = $(GRAPH_TEST_OBJ_DIR)/test_binary_kmer.o $(GRAPH_TEST_OBJ_DIR)/test_seq.o $(GRAPH_TEST_OBJ_DIR)/test_dB_graph.o $(GRAPH_TEST_OBJ_DIR)/test_dB_graph_node.o  $(GRAPH_TEST_OBJ_DIR)/run_graph_tests.o $(GRAPH_TEST_OBJ_DIR)/supernode_cmp.o $(GRAPH_TEST_OBJ_DIR)/test_hash_table.o

POP_GRAPH_TEST_OBJ = $(POP_GRAPH_TEST_OBJ_DIR)/test_binary_kmer.o $(POP_GRAPH_TEST_OBJ_DIR)/test_seq.o $(POP_GRAPH_TEST_OBJ_DIR)/test_pop_load_and_print.o $(POP_GRAPH_TEST_OBJ_DIR)/run_pop_tests.o $(POP_GRAPH_TEST_OBJ_DIR)/supernode_cmp.o $(POP_GRAPH_TEST_OBJ_DIR)/test_pop_supernode_consensus.o $(POP_GRAPH_TEST_OBJ_DIR)/test_pop_element.o


GRAPH_OBJ = $(GRAPH_OBJ_DIR)/element.o $(GRAPH_OBJ_DIR)/dB_graph.o $(GRAPH_OBJ_DIR)/seq.o $(GRAPH_OBJ_DIR)/binary_kmer.o $(GRAPH_OBJ_DIR)/priority_queue.o $(GRAPH_OBJ_DIR)/file_reader.o  $(GRAPH_OBJ_DIR)/hash_table.o  $(GRAPH_OBJ_DIR)/hash_value.o 

POP_GRAPH_OBJ = $(POP_GRAPH_OBJ_DIR)/binary_kmer.o $(POP_GRAPH_OBJ_DIR)/element.o $(POP_GRAPH_OBJ_DIR)/seq.o $(POP_GRAPH_OBJ_DIR)/hash_value.o $(POP_GRAPH_OBJ_DIR)/hash_table.o $(POP_GRAPH_OBJ_DIR)/priority_queue.o $(POP_GRAPH_OBJ_DIR)/dB_graph.o $(POP_GRAPH_OBJ_DIR)/file_reader.o  $(POP_GRAPH_OBJ_DIR)/dB_graph_population.o $(POP_GRAPH_OBJ_DIR)/pqueue_pop.o



Z_LIB_DIR= ../src/third-party/zlib-1.2.3
Z_LIB_FILES= $(Z_LIB_DIR)/*.o
ZLIB_OBJ = $(Z_LIB_DIR)/adler32.o $(Z_LIB_DIR)/compress.o $(Z_LIB_DIR)/crc32.o $(Z_LIB_DIR)/deflate.o $(Z_LIB_DIR)/gzio.o $(Z_LIB_DIR)/infback.o $(Z_LIB_DIR)/inffast.o $(Z_LIB_DIR)/inflate.o $(Z_LIB_DIR)/inftrees.o $(Z_LIB_DIR)/trees.o $(Z_LIB_DIR)/uncompr.o $(Z_LIB_DIR)/zutil.o

BIN = ../bin


## You only need the -L when running on EBI network. On the laptop you can make install CUNit so it's in /usr/local/lib or whatever
## Note you also need to buils CUNit so installs to say /homes/zi/bin, and then export LD_LIBRARY_PATH=/homes/zi/bin/lib:$LD_LIBRARY_PATH

default : graph run_pop_tests run_graph_tests zlib

clean :
	-rm $(POP_GRAPH_TEST_OBJ_DIR)/*.o  -rm $(GRAPH_TEST_OBJ_DIR)/*.o $(BIN)/run_graph_tests  $(BIN)/run_pop_tests

zlib : 
	cd $(Z_LIB_DIR); ./configure; make; rm minigzip.o; rm example.o

run_graph_tests : $(GRAPH_TEST_OBJ) $(GRAPH_OBJ) $(ZLIB_OBJ)
	$(CC) $(CFLAGS_GRAPH) $(OPT) $(LDFLAGS) -L/homes/zi/bin/lib -o $(BIN)/run_graph_tests $(GRAPH_OBJ) $(ZLIB_OBJ) $(GRAPH_TEST_OBJ) -lcunit 

run_pop_tests : $(POP_GRAPH_TEST_OBJ) $(POP_GRAPH_OBJ) $(ZLIB_OBJ)
	$(CC) $(CFLAGS_POP_GRAPH) $(OPT) $(LDFLAGS) -L/homes/zi/bin/lib -o $(BIN)/run_pop_tests $(POP_GRAPH_OBJ) $(ZLIB_OBJ) $(POP_GRAPH_TEST_OBJ) -lcunit

obj/graph_test/%.o : %.c
	$(CC) $(CFLAGS_GRAPH) $(OPT) $(LDFLAGS) -c $? -o $@

obj/pop_graph_test/%.o : %.c
	$(CC) $(CFLAGS_POP_GRAPH) $(OPT) $(LDFLAGS) -c $? -o $@

