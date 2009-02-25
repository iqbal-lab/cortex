CC = gcc

IDIR_BASIC_TESTS =../include/test/basic
IDIR_DB_GRAPH_TESTS=../include/test/graph
IDIR_SV_TRIO_TESTS=../include/test/pop_graph/sv_trio

IDIR_BASIC =../include/basic
IDIR_COUNT =../include/count
IDIR_GRAPH =../include/graph
IDIR_POP_GRAPH_CORE = ../include/pop_graph/core
IDIR_SV_TRIO = ../include/pop_graph/sv_trio
IDIR_MANY_POP_GRAPH = ../include/pop_graph/many_pops

OBJ = ./obj
OBJ_COUNT = ./obj/count
OBJ_GRAPH = ./obj/graph
OBJ_SV_TRIO = ./obj/pop_graph/sv_trio
OBJ_MANY_POP_GRAPH = ./obj/pop_graph/many_pops
OBJ_TEST = ./obj/test

BIN = ../bin

IDIR_CUNIT = /homes/zi/dev/hg/CUnit/CUnit-2.1-0/CUnit/Headers

ifdef MAC
 MACFLAG = -fnested-functions
 IDIR_CUNIT = /Users/zi/dev/hg/laptop_local/repos/CUnit/CUnit-2.1-0/CUnit/Headers
endif 


CFLAGS_COUNT = -m64 -Wall  -I$(IDIR_BASIC) -I$(IDIR_COUNT) $(MACFLAG) -O0 -g
CFLAGS_GRAPH = -m64 -Wall  -I$(IDIR_BASIC) -I$(IDIR_GRAPH) $(MACFLAG) -O0 -g
CFLAGS_SV_TRIO = -m64 -Wall  -I$(IDIR_BASIC) -I$(IDIR_POP_GRAPH_CORE) -I$(IDIR_SV_TRIO)  $(MACFLAG) -O0 -g
CFLAGS_MANY_POP_GRAPH = -m64 -Wall  -I$(IDIR_BASIC) -I$(IDIR_POP_GRAPH_CORE) -I$(IDIR_MANY_POP_GRAPH)  $(MACFLAG) -O0 -g

CFLAGS_BASIC_TESTS = -m64 -Wall -I$(IDIR_BASIC_TESTS)  -I$(IDIR_BASIC) -I$(IDIR_CUNIT) -L/homes/zi/bin/lib  $(MACFLAG) -O0 -g
CFLAGS_DB_GRAPH_TESTS = -m64 -Wall -I$(IDIR_DB_GRAPH_TESTS)  -I$(IDIR_BASIC) -I$(IDIR_GRAPH) -I$(IDIR_CUNIT) -L/homes/zi/bin/lib  $(MACFLAG) -O0 -g
CFLAGS_SV_TRIO_TESTS = -m64 -Wall -I$(IDIR_SV_TRIO_TESTS)  -I$(IDIR_BASIC) -I$(IDIR_POP_GRAPH_CORE) -I$(IDIR_SV_TRIO) -I$(IDIR_CUNIT) -L/homes/zi/bin/lib  $(MACFLAG) -O0 -g


GRAPH_OBJ = obj/graph/binary_kmer.o obj/graph/element.o obj/graph/seq.o obj/graph/hash_value.o obj/graph/hash_table.o obj/graph/priority_queue.o obj/graph/dB_graph.o obj/graph/file_reader.o obj/graph/graph.o 

GRAPH_OBJ_CL = obj/graph/binary_kmer.o obj/graph/element.o obj/graph/seq.o obj/graph/hash_value.o obj/graph/hash_table.o obj/graph/priority_queue.o obj/graph/dB_graph.o obj/graph/file_reader.o obj/graph/graph_clipping.o obj/graph/cyclic_count.o

GRAPH_OBJ_DUMP = obj/graph/binary_kmer.o obj/graph/element.o obj/graph/seq.o obj/graph/hash_value.o obj/graph/hash_table.o obj/graph/priority_queue.o obj/graph/dB_graph.o obj/graph/file_reader.o obj/graph/graph_dump.o 

GRAPH_OBJ_READ = obj/graph/binary_kmer.o obj/graph/element.o obj/graph/seq.o obj/graph/hash_value.o obj/graph/hash_table.o obj/graph/priority_queue.o obj/graph/dB_graph.o obj/graph/file_reader.o obj/graph/graph_read.o 

SV_TRIO_OBJ = obj/pop_graph/sv_trio/graph.o obj/pop_graph/sv_trio/binary_kmer.o obj/pop_graph/sv_trio/element.o obj/pop_graph/sv_trio/seq.o obj/pop_graph/sv_trio/hash_value.o obj/pop_graph/sv_trio/hash_table.o obj/pop_graph/sv_trio/priority_queue.o obj/pop_graph/sv_trio/dB_graph.o obj/pop_graph/sv_trio/file_reader.o  obj/pop_graph/sv_trio/dB_graph_population.o obj/pop_graph/sv_trio/pqueue_pop.o obj/pop_graph/sv_trio/clipping_trio.o

MANY_POP_GRAPH_OBJ = obj/pop_graph/many_pops/graph.o obj/pop_graph/many_pops/binary_kmer.o obj/pop_graph/many_pops/element.o obj/pop_graph/many_pops/seq.o obj/pop_graph/many_pops/hash_value.o obj/pop_graph/many_pops/hash_table.o obj/pop_graph/many_pops/priority_queue.o obj/pop_graph/many_pops/dB_graph.o obj/pop_graph/many_pops/file_reader.o  obj/pop_graph/many_pops/dB_graph_population.o obj/pop_graph/many_pops/pqueue_pop.o obj/pop_graph/many_pops/clipping_trio.o

COUNT_OBJ = obj/count/binary_kmer.o obj/count/element.o obj/count/seq.o obj/count/hash_value.o obj/count/hash_table.o obj/count/priority_queue.o obj/count/count_kmers.o 

BASIC_TESTS_OBJ = obj/test/basic/binary_kmer.o obj/test/basic/seq.o obj/test/basic/test_binary_kmer.o obj/test/basic/test_seq.o obj/test/basic/run_basic_tests.o

DB_GRAPH_TESTS_OBJ = obj/test/graph/binary_kmer.o obj/test/graph/element.o obj/test/graph/seq.o obj/test/graph/hash_value.o obj/test/graph/hash_table.o obj/test/graph/priority_queue.o obj/test/graph/dB_graph.o obj/test/graph/file_reader.o obj/test/graph/test_dB_graph.o obj/test/graph/run_dB_graph_tests.o obj/test/graph/test_file_reader.o obj/test/graph/test_cyclic_count.o obj/test/graph/cyclic_count.o

SV_TRIO_TESTS_OBJ = obj/test/pop_graph/sv_trio/test_pop_load_and_print.o obj/test/pop_graph/sv_trio/run_sv_trio_tests.o obj/test/pop_graph/sv_trio/supernode_cmp.o obj/test/pop_graph/sv_trio/test_pop_supernode_consensus.o obj/test/pop_graph/sv_trio/test_pop_element.o  obj/test/pop_graph/sv_trio/test_dB_graph_population.o obj/test/pop_graph/sv_trio/binary_kmer.o obj/test/pop_graph/sv_trio/element.o obj/test/pop_graph/sv_trio/seq.o obj/test/pop_graph/sv_trio/hash_value.o obj/test/pop_graph/sv_trio/hash_table.o obj/test/pop_graph/sv_trio/priority_queue.o obj/test/pop_graph/sv_trio/dB_graph.o obj/test/pop_graph/sv_trio/file_reader.o  obj/pop_graph/sv_trio/dB_graph_population.o obj/pop_graph/sv_trio/pqueue_pop.o 


sv_trio : $(SV_TRIO_OBJ)
	mkdir -p $(BIN); $(CC) $(CFLAGS_SV_TRIO) $(OPT) $(LDFLAGS) -o $(BIN)/sv_trio $(SV_TRIO_OBJ)

many_pop_graph : $(MANY_POP_GRAPH_OBJ)
	mkdir -p $(BIN); $(CC) $(CFLAGS_MANY_POP_GRAPH) $(OPT) $(LDFLAGS) -o $(BIN)/many_pop_graph $(MANY_POP_GRAPH_OBJ)

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

run_sv_trio_tests : $(SV_TRIO_TESTS_OBJ)
	mkdir -p $(BIN); $(CC) $(CFLAGS_SV_TRIO_TESTS) $(OPT) $(LDFLAGS) -o $(BIN)/run_sv_trio_tests $(SV_TRIO_TESTS_OBJ) -lcunit 


obj/pop_graph/sv_trio/%.o :  pop_graph/sv_trio/%.c
	mkdir -p obj/pop_graph/sv_trio;$(CC) $(CFLAGS_SV_TRIO) $(OPT) $(LDFLAGS) -c $? -o $@

obj/pop_graph/sv_trio/%.o :  pop_graph/core/%.c
	mkdir -p obj/pop_graph/sv_trio;$(CC) $(CFLAGS_SV_TRIO) $(OPT) $(LDFLAGS) -c $? -o $@

obj/pop_graph/sv_trio/%.o :  basic/%.c
	mkdir -p obj/pop_graph/sv_trio; $(CC) $(CFLAGS_SV_TRIO) $(OPT) $(LDFLAGS) -c $? -o $@

obj/pop_graph/many_pops/%.o :  pop_graph/many_pops/%.c
	mkdir -p obj/pop_graph/many_pops;$(CC) $(CFLAGS_MANY_POP_GRAPH) $(OPT) $(LDFLAGS) -c $? -o $@

obj/pop_graph/many_pops/%.o :  pop_graph/core/%.c
	mkdir -p obj/pop_graph/many_pops;$(CC) $(CFLAGS_MANY_POP_GRAPH) $(OPT) $(LDFLAGS) -c $? -o $@

obj/pop_graph/many_pops/%.o :  basic/%.c
	mkdir -p obj/pop_graph/many_pops;$(CC) $(CFLAGS_MANY_POP_GRAPH) $(OPT) $(LDFLAGS) -c $? -o $@

obj/graph/%.o :  basic/%.c
	mkdir -p obj/graph/; $(CC) $(CFLAGS_GRAPH) $(OPT) $(LDFLAGS) -c $? -o $@

obj/graph/%.o :  graph/%.c
	mkdir -p obj/graph/; $(CC) $(CFLAGS_GRAPH) $(OPT) $(LDFLAGS) -c $? -o $@

obj/count/%.o :  basic/%.c
	mkdir -p obj/count/; $(CC) $(CFLAGS_COUNT) $(OPT) $(LDFLAGS) -c $? -o $@

obj/test/basic/%.o :  basic/%.c 
	mkdir -p obj/test/basic; $(CC) $(CFLAGS_BASIC_TESTS) $(OPT) $(LDFLAGS) -c $? -o $@

obj/test/basic/%.o :  test/basic/%.c 
	mkdir -p obj/test/basic; $(CC) $(CFLAGS_BASIC_TESTS) $(OPT) $(LDFLAGS) -c $? -o $@

obj/test/graph/%.o :  basic/%.c 
	mkdir -p obj/test/graph; $(CC) $(CFLAGS_DB_GRAPH_TESTS) $(OPT) $(LDFLAGS) -c $? -o $@

obj/test/graph/%.o :  graph/%.c 
	mkdir -p obj/test/graph; $(CC) $(CFLAGS_DB_GRAPH_TESTS) $(OPT) $(LDFLAGS) -c $? -o $@

obj/test/graph/%.o :  test/graph/%.c 
	mkdir -p obj/test/graph; $(CC) $(CFLAGS_DB_GRAPH_TESTS) $(OPT) $(LDFLAGS) -c $? -o $@

obj/test/pop_graph/sv_trio/%.o :  basic/%.c
	mkdir -p obj/test/pop_graph/sv_trio; $(CC) $(CFLAGS_SV_TRIO_TESTS) $(OPT) $(LDFLAGS) -c $? -o $@

obj/test/pop_graph/sv_trio/%.o :  pop_graph/core/%.c 
	mkdir -p obj/test/pop_graph/sv_trio; $(CC) $(CFLAGS_SV_TRIO_TESTS) $(OPT) $(LDFLAGS) -c $? -o $@

obj/test/pop_graph/sv_trio/%.o :  pop_graph/sv_trio/%.c 
	mkdir -p obj/test/pop_graph/sv_trio; $(CC) $(CFLAGS_SV_TRIO_TESTS) $(OPT) $(LDFLAGS) -c $? -o $@

obj/test/pop_graph/sv_trio/%.o :  test/pop_graph/sv_trio/%.c 
	mkdir -p obj/test/pop_graph/sv_trio; $(CC) $(CFLAGS_SV_TRIO_TESTS) $(OPT) $(LDFLAGS) -c $? -o $@


