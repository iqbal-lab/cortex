CC = gcc

IDIR_TEST =../include/test
IDIR_BASIC =../include/basic
IDIR_COUNT =../include/count
IDIR_GRAPH =../include/graph
IDIR_POP_GRAPH = ../include/pop_graph

OBJ = ./obj
OBJ_COUNT = ./obj/count
OBJ_GRAPH = ./obj/graph
OBJ_POP_GRAPH = ./obj/pop_graph
OBJ_TEST = ./obj/test

BIN = ../bin

IDIR_CUNIT = /homes/zi/dev/hg/CUnit/CUnit-2.1-0/CUnit/Headers

ifdef MAC
 MACFLAG = -fnested-functions
 IDIR_CUNIT = /Users/zi/dev/hg/laptop_local/repos/CUnit/CUnit-2.1-0/CUnit/Headers
endif 


CFLAGS_COUNT = -m64 -Wall  -I$(IDIR_BASIC) -I$(IDIR_COUNT) $(MACFLAG) -O0 -g
CFLAGS_GRAPH = -m64 -Wall  -I$(IDIR_BASIC) -I$(IDIR_GRAPH) $(MACFLAG) -O0 -g
CFLAGS_POP_GRAPH = -m64 -Wall  -I$(IDIR_BASIC) -I$(IDIR_POP_GRAPH) $(MACFLAG) -O0 -g

CFLAGS_BASIC_TESTS = -m64 -Wall -I$(IDIR_TEST)  -I$(IDIR_BASIC) -I$(IDIR_CUNIT) -L/homes/zi/bin/lib  $(MACFLAG) -O0 -g

CFLAGS_DB_GRAPH_TESTS = -m64 -Wall -I$(IDIR_TEST)  -I$(IDIR_BASIC) -I$(IDIR_GRAPH) -I$(IDIR_CUNIT) -L/homes/zi/bin/lib  $(MACFLAG) -O0 -g


GRAPH_OBJ = obj/graph/binary_kmer.o obj/graph/element.o obj/graph/seq.o obj/graph/hash_value.o obj/graph/hash_table.o obj/graph/priority_queue.o obj/graph/dB_graph.o obj/graph/file_reader.o obj/graph/graph.o 

GRAPH_OBJ_CL = obj/graph/binary_kmer.o obj/graph/element.o obj/graph/seq.o obj/graph/hash_value.o obj/graph/hash_table.o obj/graph/priority_queue.o obj/graph/dB_graph.o obj/graph/file_reader.o obj/graph/graph_clipping.o 

GRAPH_OBJ_DUMP = obj/graph/binary_kmer.o obj/graph/element.o obj/graph/seq.o obj/graph/hash_value.o obj/graph/hash_table.o obj/graph/priority_queue.o obj/graph/dB_graph.o obj/graph/file_reader.o obj/graph/graph_dump.o 

GRAPH_OBJ_READ = obj/graph/binary_kmer.o obj/graph/element.o obj/graph/seq.o obj/graph/hash_value.o obj/graph/hash_table.o obj/graph/priority_queue.o obj/graph/dB_graph.o obj/graph/file_reader.o obj/graph/graph_read.o 

POP_GRAPH_OBJ = obj/pop_graph/graph.o obj/pop_graph/binary_kmer.o obj/pop_graph/element.o obj/pop_graph/seq.o obj/pop_graph/hash_value.o obj/pop_graph/hash_table.o obj/pop_graph/priority_queue.o obj/pop_graph/dB_graph.o obj/pop_graph/file_reader.o  obj/pop_graph/dB_graph_population.o obj/pop_graph/pqueue_pop.o

COUNT_OBJ = obj/count/binary_kmer.o obj/count/element.o obj/count/seq.o obj/count/hash_value.o obj/count/hash_table.o obj/count/priority_queue.o obj/count/count_kmers.o 

BASIC_TESTS_OBJ = obj/test/basic/binary_kmer.o obj/test/basic/seq.o obj/test/basic/test_binary_kmer.o obj/test/basic/test_seq.o obj/test/basic/run_basic_tests.o

DB_GRAPH_TESTS_OBJ = obj/test/graph/binary_kmer.o obj/test/graph/element.o obj/test/graph/seq.o obj/test/graph/hash_value.o obj/test/graph/hash_table.o obj/test/graph/priority_queue.o obj/test/graph/dB_graph.o obj/test/graph/file_reader.o obj/test/graph/test_dB_graph.o obj/test/graph/run_dB_graph_tests.o

pop_graph : $(POP_GRAPH_OBJ)
	mkdir -p $(BIN); $(CC) $(CFLAGS_POP_GRAPH) $(OPT) $(LDFLAGS) -o $(BIN)/pop_graph $(POP_GRAPH_OBJ)

graph : $(GRAPH_OBJ)
	mkdir -p $(BIN); $(CC) $(CFLAGS_GRAPH) $(OPT) $(LDFLAGS) -o $(BIN)/graph $(GRAPH_OBJ)

graph_cl : $(GRAPH_OBJ_CL)
	mkdir -p $(BIN); $(CC) $(CFLAGS_GRAPH) $(OPT) $(LDFLAGS) -o $(BIN)/graph_cp $(GRAPH_OBJ_CL)

graph_dump : $(GRAPH_OBJ_DUMP)
	mkdir -p $(BIN); $(CC) $(CFLAGS_GRAPH) $(OPT) $(LDFLAGS) -o $(BIN)/graph_dump $(GRAPH_OBJ_DUMP)

graph_read : $(GRAPH_OBJ_READ)
	mkdir -p $(BIN); $(CC) $(CFLAGS_GRAPH) $(OPT) $(LDFLAGS) -o $(BIN)/graph_read $(GRAPH_OBJ_READ)


count : $(COUNT_OBJ)
	mkdir -p $(BIN); $(CC) $(CFLAGS_COUNT) $(OPT) $(LDFLAGS) -o $(BIN)/count_kmers $(COUNT_OBJ)

run_basic_tests : $(BASIC_TESTS_OBJ)
	mkdir -p $(BIN); $(CC) $(CFLAGS_BASIC_TESTS) $(OPT) $(LDFLAGS) -o $(BIN)/run_basic_tests $(BASIC_TESTS_OBJ) -lcunit

run_dB_graph_tests : $(DB_GRAPH_TESTS_OBJ)
	mkdir -p $(BIN); $(CC) $(CFLAGS_DB_GRAPH_TESTS) $(OPT) $(LDFLAGS) -o $(BIN)/run_dB_graph_tests $(DB_GRAPH_TESTS_OBJ) -lcunit 

obj/pop_graph/%.o :  pop_graph/%.c
	mkdir -p obj/pop_graph;$(CC) $(CFLAGS_POP_GRAPH) $(OPT) $(LDFLAGS) -c $? -o $@

obj/pop_graph/%.o :  basic/%.c
	mkdir -p obj/pop_graph/; $(CC) $(CFLAGS_POP_GRAPH) $(OPT) $(LDFLAGS) -c $? -o $@

obj/graph/%.o :  basic/%.c
	mkdir -p obj/graph/; $(CC) $(CFLAGS_GRAPH) $(OPT) $(LDFLAGS) -c $? -o $@

obj/graph/%.o :  graph/%.c
	mkdir -p obj/graph/; $(CC) $(CFLAGS_GRAPH) $(OPT) $(LDFLAGS) -c $? -o $@

obj/count/%.o :  basic/%.c
	mkdir -p obj/count/; $(CC) $(CFLAGS_COUNT) $(OPT) $(LDFLAGS) -c $? -o $@

obj/test/basic/%.o :  basic/%.c 
	mkdir -p obj/test/basic; $(CC) $(CFLAGS_BASIC_TESTS) $(OPT) $(LDFLAGS) -c $? -o $@

obj/test/basic/%.o :  test/%.c 
	mkdir -p obj/test/basic; $(CC) $(CFLAGS_BASIC_TESTS) $(OPT) $(LDFLAGS) -c $? -o $@

obj/test/graph/%.o :  basic/%.c 
	mkdir -p obj/test/graph; $(CC) $(CFLAGS_DB_GRAPH_TESTS) $(OPT) $(LDFLAGS) -c $? -o $@

obj/test/graph/%.o :  graph/%.c 
	mkdir -p obj/test/graph; $(CC) $(CFLAGS_DB_GRAPH_TESTS) $(OPT) $(LDFLAGS) -c $? -o $@

obj/test/graph/%.o :  test/%.c 
	mkdir -p obj/test/graph; $(CC) $(CFLAGS_DB_GRAPH_TESTS) $(OPT) $(LDFLAGS) -c $? -o $@

