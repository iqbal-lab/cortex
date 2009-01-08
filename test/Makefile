include ../Makefile.inc

IDIR_BASIC =/homes/zi/dev/hg/marzam/marzam-harmonised/include/basic
IDIR_GRAPH =/homes/zi/dev/hg/marzam/marzam-harmonised/include/graph
IDIR_CUNIT =/homes/zi/dev/hg/CUnit/CUnit-2.1-0/CUnit/Headers

CC = gcc
CFLAGS = -Wall -I$(IDIR_BASIC) -I$(IDIR_GRAPH) -I$(IDIR_CUNIT) -g -O0


TEST_OBJ_DIR=obj
TEST_OBJ = $(TEST_OBJ_DIR)/test_binary_kmer.o $(TEST_OBJ_DIR)/test_dB_graph.o $(TEST_OBJ_DIR)/test_dB_graph_node.o $(TEST_OBJ_DIR)/test_seq.o  $(TEST_OBJ_DIR)/run_tests.o

GRAPH_OBJ_DIR=../src/obj/graph
GRAPH_OBJ = $(GRAPH_OBJ_DIR)/element.o $(GRAPH_OBJ_DIR)/dB_graph.o $(GRAPH_OBJ_DIR)/seq.o $(GRAPH_OBJ_DIR)/binary_kmer.o $(GRAPH_OBJ_DIR)/priority_queue.o $(GRAPH_OBJ_DIR)/file_reader.o  $(GRAPH_OBJ_DIR)/hash_table.o  $(GRAPH_OBJ_DIR)/hash_value.o

Z_LIB_DIR= ../src/third-party/zlib-1.2.3
Z_LIB_FILES= $(Z_LIB_DIR)/*.o
ZLIB_OBJ = $(Z_LIB_DIR)/adler32.o $(Z_LIB_DIR)/compress.o $(Z_LIB_DIR)/crc32.o $(Z_LIB_DIR)/deflate.o $(Z_LIB_DIR)/gzio.o $(Z_LIB_DIR)/infback.o $(Z_LIB_DIR)/inffast.o $(Z_LIB_DIR)/inflate.o $(Z_LIB_DIR)/inftrees.o $(Z_LIB_DIR)/trees.o $(Z_LIB_DIR)/uncompr.o $(Z_LIB_DIR)/zutil.o

BIN = ../bin


## You only need the -L when running on EBI network. On the laptop you can make install CUNit so it's in /usr/local/lib or whatever
## Note you also need to buils CUNit so installs to say /homes/zi/bin, and then export LD_LIBRARY_PATH=/homes/zi/bin/lib:$LD_LIBRARY_PATH

default : graph run_tests zlib

clean :
	-rm $(TEST_OBJ_DIR)/*.o $(BIN)/run_tests

zlib : 
	cd $(Z_LIB_DIR); ./configure; make; rm minigzip.o; rm example.o

run_tests : $(TEST_OBJ) $(GRAPH_OBJ) $(ZLIB_OBJ)
	$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -L/homes/zi/bin/lib -o $(BIN)/run_tests $(GRAPH_OBJ) $(ZLIB_OBJ) $(TEST_OBJ) -lcunit 

obj/%.o : %.c
	$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -c $? -o $@

