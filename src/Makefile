CC = gcc

IDIR_TEST  =../include/test
IDIR_BASIC_TESTS =../include/test/basic
IDIR_GRAPH_COUNT_TESTS =../include/test/graph_count
IDIR_DB_GRAPH_TESTS=../include/test/graph
IDIR_SV_TRIO_TESTS=../include/test/pop_graph/sv_trio
IDIR_MANY_POP_TESTS=../include/test/pop_graph/many_pops

IDIR_BASIC =../include/basic
IDIR_HASH  =../include/basic/hash_table 
IDIR_GRAPH =../include/graph
IDIR_GRAPH_COUNT = ../include/graph_count
IDIR_POP_GRAPH_CORE = ../include/pop_graph/core
IDIR_SV_TRIO = ../include/pop_graph/sv_trio
IDIR_MANY_POP_GRAPH = ../include/pop_graph/many_pops

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

CFLAGS_GRAPH_COUNT = -m64 -Wall  -I$(IDIR_BASIC)  -I$(IDIR_HASH) -I$(IDIR_GRAPH_COUNT) $(MACFLAG) -O0 -g
CFLAGS_POP_GRAPH_CORE = -m64 -Wall  -I$(IDIR_BASIC)  -I$(IDIR_HASH) -I$(IDIR_POP_GRAPH_CORE)  $(MACFLAG) -O0 -g
CFLAGS_SV_TRIO = -m64 -Wall  -I$(IDIR_BASIC)  -I$(IDIR_HASH) -I$(IDIR_POP_GRAPH_CORE) -I$(IDIR_SV_TRIO)  $(MACFLAG) -O0 -g
CFLAGS_MANY_POP_GRAPH = -m64 -Wall  -I$(IDIR_BASIC)  -I$(IDIR_HASH) -I$(IDIR_POP_GRAPH_CORE) -I$(IDIR_MANY_POP_GRAPH)  $(MACFLAG) -O0 -g

CFLAGS_BASIC_TESTS = -m64 -Wall -I$(IDIR_BASIC_TESTS)  -I$(IDIR_BASIC)  -I$(IDIR_HASH) -I$(IDIR_CUNIT) -L/homes/zi/bin/lib  $(MACFLAG) -O0 -g
CFLAGS_COUNT_TESTS = -m64 -Wall -I$(IDIR_COUNT_TESTS)  -I$(IDIR_COUNT) -I$(IDIR_BASIC)  -I$(IDIR_HASH) -I$(IDIR_CUNIT) -L/homes/zi/bin/lib  $(MACFLAG) -O0 -g
CFLAGS_GRAPH_COUNT_TESTS = -m64 -Wall -I$(IDIR_GRAPH_COUNT_TESTS)  -I$(IDIR_GRAPH_COUNT) -I$(IDIR_BASIC) -I$(IDIR_HASH) -I$(IDIR_CUNIT) -L/homes/zi/bin/lib  $(MACFLAG) -O0 -g

CFLAGS_DB_GRAPH_TESTS = -m64 -Wall -I$(IDIR_DB_GRAPH_TESTS)  -I$(IDIR_BASIC) -I$(IDIR_HASH) -I$(IDIR_GRAPH) -I$(IDIR_CUNIT) -L/homes/zi/bin/lib  $(MACFLAG) -O0 -g
CFLAGS_SV_TRIO_TESTS = -m64 -Wall -I$(IDIR_SV_TRIO_TESTS)  -I$(IDIR_BASIC) -I$(IDIR_HASH) -I$(IDIR_POP_GRAPH_CORE) -I$(IDIR_SV_TRIO) -I$(IDIR_CUNIT) -L/homes/zi/bin/lib  $(MACFLAG) -O0 -g
CFLAGS_MANY_POP_TESTS = -m64 -Wall -I$(IDIR_MANY_POP_TESTS)  -I$(IDIR_BASIC) -I$(IDIR_HASH) -I$(IDIR_POP_GRAPH_CORE) -I$(IDIR_MANY_POP_GRAPH) -I$(IDIR_CUNIT) -L/homes/zi/bin/lib  $(MACFLAG) -O0 -g




GRAPH_OBJ_DUMP = obj/basic/binary_kmer.o obj/graph/element.o obj/basic/seq.o obj/basic/hash_table/hash_key/bob_jenkins/hash_value.o obj/basic/hash_table/open_hash/hash_table.o obj/graph/dB_graph.o obj/graph/file_reader.o obj/graph/graph_dump.o 

GRAPH_OBJ      = obj/basic/binary_kmer.o obj/graph/element.o obj/basic/seq.o obj/basic/hash_table/hash_key/bob_jenkins/hash_value.o obj/basic/hash_table/open_hash/hash_table.o obj/graph/dB_graph.o obj/graph/file_reader.o obj/graph/graph.o obj/graph/pq_graph.o obj/graph/cyclic_count.o

BASIC_TESTS_OBJ = obj/basic/binary_kmer.o obj/basic/seq.o obj/test/basic/test_binary_kmer.o obj/test/basic/test_seq.o obj/test/basic/run_basic_tests.o

DB_GRAPH_TESTS_OBJ = obj/basic/binary_kmer.o obj/graph/element.o obj/basic/seq.o obj/basic/hash_table/hash_key/bob_jenkins/hash_value.o obj/basic/hash_table/open_hash/hash_table.o obj/graph/dB_graph.o obj/graph/file_reader.o obj/test/graph/test_dB_graph.o  obj/test/graph/test_file_reader.o obj/test/graph/run_dB_graph_tests.o obj/graph/pq_graph.o obj/test/graph/test_cyclic_count.o obj/graph/cyclic_count.o

SV_TRIO_OBJ = obj/pop_graph/sv_trio/graph.o obj/basic/binary_kmer.o obj/pop_graph/sv_trio/element.o obj/basic/seq.o obj/basic/hash_table/hash_key/bob_jenkins/hash_value.o obj/basic/hash_table/open_hash/hash_table.o obj/pop_graph/core/dB_graph.o obj/pop_graph/core/file_reader.o  obj/pop_graph/core/dB_graph_population.o obj/pop_graph/core/pqueue_pop.o obj/pop_graph/core/clipping_trio.o






graph_dump : $(GRAPH_OBJ_DUMP)
	mkdir -p $(BIN); $(CC) $(OPT) $(LDFLAGS) -o $(BIN)/graph_dump $(GRAPH_OBJ_DUMP)

graph : $(GRAPH_OBJ)
	mkdir -p $(BIN); $(CC) $(OPT) $(LDFLAGS) -o $(BIN)/graph $(GRAPH_OBJ)

run_basic_tests : $(BASIC_TESTS_OBJ)
	mkdir -p $(BIN); $(CC) $(OPT) $(CFLAGS_CUNIT) $(LDFLAGS) -o $(BIN)/run_basic_tests $(BASIC_TESTS_OBJ) -lcunit

run_dB_graph_tests : $(DB_GRAPH_TESTS_OBJ)
	mkdir -p $(BIN); $(CC) $(OPT) $(CFLAGS_CUNIT) $(LDFLAGS) -o $(BIN)/run_dB_graph_tests $(DB_GRAPH_TESTS_OBJ) -lcunit 

run_sv_trio_tests : $(SV_TRIO_TESTS_OBJ)
	mkdir -p $(BIN); $(CC) $(CFLAGS_SV_TRIO_TESTS) $(OPT) $(LDFLAGS) -o $(BIN)/run_sv_trio_tests $(SV_TRIO_TESTS_OBJ) -lcunit



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

obj/test/basic/%.o :  test/basic/%.c 
	mkdir -p obj/test/basic; $(CC) $(CFLAGS_GRAPH) $(CFLAGS_BASIC_TESTS) $(OPT) $(LDFLAGS) -c $? -o $@

obj/test/graph/%.o :  test/graph/%.c 
	mkdir -p obj/test/graph; $(CC) $(CFLAGS_GRAPH) $(CFLAGS_DB_GRAPH_TESTS) $(OPT) $(LDFLAGS) -c $? -o $@

obj/pop_graph/core/%.o :  pop_graph/core/%.c 
	mkdir -p obj/test/pop_graph/core; $(CC) $(CFLAGS_GRAPH) $(CFLAGS_POP_GRAPH_CORE) $(OPT) $(LDFLAGS) -c $? -o $@

obj/pop_graph/sv_trio/%.o :  pop_graph/sv_trio/%.c 
	mkdir -p obj/test/pop_graph/sv_trio; $(CC) $(CFLAGS_GRAPH) $(CFLAGS_SV_TRIO) $(OPT) $(LDFLAGS) -c $? -o $@

