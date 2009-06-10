CC = gcc

IDIR_BASIC_TESTS =../include/test/basic
IDIR_COUNT_KMERS_TESTS =../include/test/count_kmers
IDIR_GRAPH_TESTS=../include/test/graph
IDIR_SV_TRIO_TESTS=../include/test/pop_graph/sv_trio
IDIR_MANY_POP_TESTS=../include/test/pop_graph/many_pops

IDIR_BASIC =../include/basic
IDIR_HASH  =../include/hash_table 
IDIR_GRAPH =../include/graph
IDIR_COUNT_KMERS = ../include/count_kmers
IDIR_SV_TRIO = ../include/pop_graph/sv_trio 
IDIR_POP_GRAPH_CORE = ../include/pop_graph/core
IDIR_MANY_POP_GRAPH = ../include/pop_graph/many_pops ../include/pop_graph/core


IDIR_CUNIT = /homes/zi/dev/hg/CUnit/CUnit-2.1-0/CUnit/Headers

BIN = ../bin


ifdef MAC
 MACFLAG    = -fnested-functions
 IDIR_CUNIT = /Users/zi/dev/hg/laptop_local/repos/CUnit/CUnit-2.1-0/CUnit/Headers
endif 

OPT  		      = -m64 -Wall -O3 $(MACFLAG)
ifdef DEBUG
 OPT                  = -m64 -Wall -O0 -g $(MACFLAG)
endif

CFLAGS_CUNIT          = -L/homes/zi/bin/lib

CFLAGS_GRAPH          = -I$(IDIR_BASIC) -I$(IDIR_HASH) -I$(IDIR_GRAPH) 
CFLAGS_COUNT_KMERS    = -I$(IDIR_BASIC)  -I$(IDIR_HASH) -I$(IDIR_COUNT_KMERS) 
CFLAGS_SV_TRIO =   -I$(IDIR_BASIC)  -I$(IDIR_HASH)  -I$(IDIR_SV_TRIO) -I$(IDIR_POP_GRAPH_CORE)
CFLAGS_MANY_POP_GRAPH =   -I$(IDIR_BASIC)  -I$(IDIR_HASH)  -I$(IDIR_MANY_POP_GRAPH) -I$(IDIR_POP_GRAPH_CORE)
CFLAGS_BASIC_TESTS = -I$(IDIR_BASIC_TESTS)  -I$(IDIR_BASIC)  -I$(IDIR_CUNIT) -L/homes/zi/bin/lib
CFLAGS_COUNT_KMERS_TESTS =  -I$(IDIR_COUNT_KMERS_TESTS)  -I$(IDIR_COUNT_KMERS) -I$(IDIR_BASIC)  -I$(IDIR_HASH) -I$(IDIR_CUNIT) -L/homes/zi/bin/lib
CFLAGS_GRAPH_TESTS = -I$(IDIR_GRAPH_TESTS)  -I$(IDIR_BASIC) -I$(IDIR_HASH) -I$(IDIR_GRAPH) -I$(IDIR_CUNIT) -L/homes/zi/bin/lib 
CFLAGS_SV_TRIO_TESTS =  -I$(IDIR_SV_TRIO_TESTS)  -I$(IDIR_BASIC) -I$(IDIR_HASH)  -I$(IDIR_SV_TRIO) -I$(IDIR_POP_GRAPH_CORE) -I$(IDIR_CUNIT) -L/homes/zi/bin/lib
CFLAGS_MANY_POP_TESTS =  -I$(IDIR_MANY_POP_TESTS)  -I$(IDIR_BASIC) -I$(IDIR_HASH)  -I$(IDIR_MANY_POP_GRAPH) -I$(IDIR_POP_GRAPH_CORE) -I$(IDIR_CUNIT) -L/homes/zi/bin/lib


# Shared objects for all binaries using single colour graphs defined by element in src/graph/. Objects in obj/basic can be shared by all targets. Other objects must go in target specific obj-hierarchy

GRAPH_OBJ      = obj/basic/binary_kmer.o obj/basic/seq.o obj/graph/element.o obj/graph/hash_table/hash_key/bob_jenkins/hash_value.o obj/graph/hash_table/open_hash/hash_table.o obj/graph/dB_graph.o obj/graph/file_reader.o 
SV_TRIO_OBJ = obj/pop_graph/sv_trio/graph.o obj/pop_graph/sv_trio/binary_kmer.o obj/pop_graph/sv_trio/element.o obj/pop_graph/sv_trio/seq.o obj/pop_graph/sv_trio/hash_value.o obj/pop_graph/sv_trio/hash_table.o obj/pop_graph/sv_trio/dB_graph.o obj/pop_graph/sv_trio/file_reader.o  obj/pop_graph/sv_trio/dB_graph_population.o obj/pop_graph/sv_trio/pqueue_pop.o obj/pop_graph/sv_trio/clipping_trio.o


BASIC_TESTS_OBJ = obj/basic/binary_kmer.o obj/basic/seq.o obj/test/basic/test_binary_kmer.o obj/test/basic/test_seq.o obj/test/basic/run_basic_tests.o
GRAPH_TESTS_OBJ = obj/basic/binary_kmer.o obj/basic/seq.o obj/graph/element.o obj/graph/hash_table/hash_key/bob_jenkins/hash_value.o obj/graph/hash_table/open_hash/hash_table.o obj/graph/dB_graph.o obj/graph/file_reader.o obj/test/graph/test_dB_graph.o  obj/test/graph/test_file_reader.o obj/test/graph/run_dB_graph_tests.o obj/test/graph/test_cyclic_count.o obj/graph/cyclic_count.o
SV_TRIO_TESTS_OBJ = obj/test/pop_graph/sv_trio/test_file_reader.o obj/test/pop_graph/sv_trio/test_pop_load_and_print.o obj/test/pop_graph/sv_trio/run_sv_trio_tests.o obj/test/pop_graph/sv_trio/supernode_cmp.o obj/test/pop_graph/sv_trio/test_pop_supernode_consensus.o obj/test/pop_graph/sv_trio/test_pop_element.o  obj/test/pop_graph/sv_trio/test_dB_graph_population.o obj/pop_graph/sv_trio/binary_kmer.o obj/pop_graph/sv_trio/element.o obj/pop_graph/sv_trio/seq.o obj/pop_graph/sv_trio/hash_value.o obj/pop_graph/sv_trio/hash_table.o obj/pop_graph/sv_trio/dB_graph.o obj/pop_graph/sv_trio/file_reader.o  obj/pop_graph/sv_trio/dB_graph_population.o



graph_dump : $(GRAPH_OBJ) obj/graph/graph_dump.o
	mkdir -p $(BIN); $(CC) $(OPT)  -o $(BIN)/graph_dump $(GRAPH_OBJ) obj/graph/graph_dump.o

graph : $(GRAPH_OBJ) obj/graph/graph.o
	mkdir -p $(BIN); $(CC) $(OPT)  -o $(BIN)/graph $(GRAPH_OBJ) obj/graph/graph.o

sv_trio : $(SV_TRIO_OBJ)
	mkdir -p $(BIN); $(CC) $(CFLAGS_SV_TRIO) $(OPT)  -o $(BIN)/sv_trio $(SV_TRIO_OBJ)

run_basic_tests : $(BASIC_TESTS_OBJ)
	mkdir -p $(BIN); $(CC) $(OPT) $(CFLAGS_CUNIT)  -o $(BIN)/run_basic_tests $(BASIC_TESTS_OBJ) -lcunit

run_graph_tests : $(GRAPH_TESTS_OBJ)
	mkdir -p $(BIN); $(CC) $(OPT) $(CFLAGS_CUNIT)  -o $(BIN)/run_graph_tests $(GRAPH_TESTS_OBJ) -lcunit 

run_sv_trio_tests : $(SV_TRIO_TESTS_OBJ)
	mkdir -p $(BIN); $(CC) $(CFLAGS_SV_TRIO_TESTS) $(OPT)  -o $(BIN)/run_sv_trio_tests $(SV_TRIO_TESTS_OBJ) -lcunit


#pattern rules

obj/graph/hash_table/open_hash/%.o :  hash_table/open_hash/%.c 
	mkdir -p obj/graph/hash_table/open_hash; $(CC) $(CFLAGS_GRAPH) $(OPT)  -c $? -o $@

obj/graph/hash_table/hash_key/bob_jenkins/%.o :  hash_table/hash_key/bob_jenkins/%.c
	mkdir -p obj/graph/hash_table/hash_key/bob_jenkins; $(CC) $(CFLAGS_GRAPH) $(OPT)  -c $? -o $@

obj/graph/hash_table/hash_key/crc32/%.o :  hash_table/hash_key/crc32/%.c
	mkdir -p obj/graph/hash_table/hash_key/crc32; $(CC) $(CFLAGS_GRAPH) $(OPT)  -c $? -o $@

obj/pop_graph/sv_trio/hash_table/open_hash/%.o :  hash_table/open_hash/%.c 
	mkdir -p obj/pop_graph/sv_trio/hash_table/open_hash; $(CC) $(CFLAGS_GRAPH) $(OPT)  -c $? -o $@

obj/pop_graph/sv_trio/hash_table/hash_key/bob_jenkins/%.o :  hash_table/hash_key/bob_jenkins/%.c
	mkdir -p obj/pop_graph/sv_trio/hash_table/hash_key/bob_jenkins; $(CC) $(CFLAGS_GRAPH) $(OPT)  -c $? -o $@

obj/pop_graph/sv_trio/hash_table/hash_key/crc32/%.o :  hash_table/hash_key/crc32/%.c
	mkdir -p obj/pop_graph/sv_trio/hash_table/hash_key/crc32; $(CC) $(CFLAGS_GRAPH) $(OPT)  -c $? -o $@

obj/basic/%.o : basic/%.c 
	mkdir -p obj/basic/; $(CC) $(CFLAGS_GRAPH) $(OPT)  -c $? -o $@

obj/test/basic/%.o :  test/basic/%.c 
	mkdir -p obj/test/basic; $(CC) $(CFLAGS_GRAPH) $(CFLAGS_BASIC_TESTS) $(OPT)  -c $? -o $@

obj/graph/%.o :  graph/%.c
	mkdir -p obj/graph/; $(CC) $(CFLAGS_GRAPH) $(OPT)  -c $? -o $@

obj/test/graph/%.o :  test/graph/%.c 
	mkdir -p obj/test/graph; $(CC) $(CFLAGS_GRAPH) $(CFLAGS_GRAPH_TESTS) $(OPT)  -c $? -o $@

obj/pop_graph/sv_trio/%.o :  basic/%.c 
	mkdir -p obj/pop_graph/sv_trio; $(CC) $(CFLAGS_SV_TRIO) $(OPT)  -c $? -o $@

obj/pop_graph/sv_trio/%.o :  hash_table/hash_key/bob_jenkins/%.c 
	mkdir -p obj/pop_graph/sv_trio; $(CC) $(CFLAGS_SV_TRIO) $(OPT)  -c $? -o $@

obj/pop_graph/sv_trio/%.o :  hash_table/open_hash/%.c 
	mkdir -p obj/pop_graph/sv_trio; $(CC) $(CFLAGS_SV_TRIO) $(OPT)  -c $? -o $@

obj/pop_graph/sv_trio/%.o :  pop_graph/core/%.c 
	mkdir -p obj/pop_graph/sv_trio; $(CC) $(CFLAGS_SV_TRIO) $(OPT)  -c $? -o $@

obj/pop_graph/sv_trio/%.o :  pop_graph/sv_trio/%.c 
	mkdir -p obj/pop_graph/sv_trio; $(CC) $(CFLAGS_SV_TRIO) $(OPT)  -c $? -o $@

obj/test/pop_graph/sv_trio/%.o :  test/pop_graph/sv_trio/%.c 
	mkdir -p obj/test/pop_graph/sv_trio; $(CC)  $(CFLAGS_SV_TRIO_TESTS) $(OPT)  -c $? -o $@
