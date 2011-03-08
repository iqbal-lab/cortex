ifndef CC
  CC = gcc	
endif

BIN = bin

ifeq ($(MAXK),31)
   BITFIELDS = 1
endif

ifeq ($(MAXK),63)
   BITFIELDS = 2
endif

ifeq ($(MAXK),95)
   BITFIELDS = 3
endif

ifeq ($(MAXK),127)
   BITFIELDS = 4
endif

ifeq ($(MAXK),160)
   BITFIELDS = 5
endif

ifeq ($(MAXK),192)
   BITFIELDS = 6
endif

ifeq ($(MAXK),223)
   BITFIELDS = 7
endif

ifeq ($(MAXK),255)
   BITFIELDS = 8
endif

ifndef BITFIELDS
   BITFIELDS = 1
   MAXK = 31	
endif


ifndef NUM_COLS
  NUM_COLS = 1
endif

#main program includes
IDIR_BASIC =include/basic
IDIR_BASE_ENCODING = ${IDIR_BASIC}/event_encoding/base_encoding
IDIR_COLOUR_ENCODING = ${IDIR_BASIC}/event_encoding/solid_colour_encoding
IDIR_HASH  =include/hash_table 
IDIR_CORTEX_CON = include/cortex_con
IDIR_CORTEX_VAR =include/cortex_var/many_colours 
IDIR_CORTEX_VAR_CORE = include/cortex_var/core
IDIR_CORTEX_VAR_CMD_LINE =include/cortex_var/many_colours

## test code includes
IDIR_BASIC_TESTS =include/test/basic
IDIR_HASH_TABLE_TESTS =include/test/hash_table
IDIR_CORTEX_CON_TESTS=include/test/graph
IDIR_CORTEX_VAR_TESTS=include/test/cortex_var/many_colours

IDIR_CUNIT = /home/zam/dev/hg/CUnit/CUnit-2.1-0/CUnit/Headers


ifdef MAC
 MACFLAG    = -fnested-functions
 IDIR_CUNIT = /Users/zi/dev/hg/laptop_local/repos/CUnit/CUnit-2.1-0/CUnit/Headers
endif 


ARCH =  -m64 


ifdef 32_BITS
 ARCH =  
endif
OPT  		      = $(ARCH) -Wall  -O3  $(MACFLAG) -DNUMBER_OF_BITFIELDS_IN_BINARY_KMER=$(BITFIELDS) -DNUMBER_OF_COLOURS=$(NUM_COLS)

ifdef DEBUG
 OPT                  = $(ARCH) -Wall -O0 -g $(MACFLAG) -DNUMBER_OF_BITFIELDS_IN_BINARY_KMER=$(BITFIELDS) -DNUMBER_OF_COLOURS=$(NUM_COLS)
endif



CFLAGS_CUNIT          = -L/home/zam/bin/lib

CFLAGS_BASIC          = -I$(IDIR_BASIC) -I$(IDIR_BASE_ENCODING)
CFLAGS_GRAPH          = -I$(IDIR_BASIC) -I$(IDIR_HASH) -I$(IDIR_CORTEX_CON) -I$(IDIR_BASE_ENCODING)
CFLAGS_CORTEX_VAR = -I$(IDIR_CORTEX_VAR_CORE) -I$(IDIR_BASIC)  -I$(IDIR_HASH)  -I$(IDIR_CORTEX_VAR) -I$(IDIR_BASE_ENCODING) 


CFLAGS_BASIC_TESTS = -I$(IDIR_BASIC_TESTS)  -I$(IDIR_BASIC)  -I$(IDIR_BASE_ENCODING) -I$(IDIR_CUNIT) -L/home/zam/bin/lib
CFLAGS_HASH_TABLE_TESTS =  -I$(IDIR_HASH) -I$(IDIR_HASH_TABLE_TESTS) -I$(IDIR_CUNIT)  -I$(IDIR_CORTEX_VAR) -L/home/zam/bin/lib
CFLAGS_GRAPH_TESTS = -I$(IDIR_GRAPH_TESTS)  -I$(IDIR_BASIC) -I$(IDIR_BASE_ENCODING) -I$(IDIR_HASH) -I$(IDIR_CORTEX_CON) -I$(IDIR_CUNIT) -L/home/zam/bin/lib 
CFLAGS_CORTEX_VAR_TESTS = -I$(IDIR_CORTEX_VAR_TESTS) -I$(IDIR_BASIC) -I$(IDIR_BASE_ENCODING) -I$(IDIR_HASH) -I$(IDIR_CORTEX_VAR) -I$(IDIR_CORTEX_VAR_CORE) -I$(IDIR_CUNIT) -L/home/zam/bin/lib
CFLAGS_CORTEX_VAR_CMD_LINE_TESTS = -I$(IDIR_CORTEX_VAR_TESTS) -I$(IDIR_CORTEX_VAR_CMD_LINE) -I$(IDIR_BASIC) -I$(IDIR_BASE_ENCODING) -I$(IDIR_CUNIT) -L/home/zam/bin/lib


GRAPH_OBJ      = src/obj/cortex_con/cmd_line.o src/obj/cortex_con/binary_kmer.o  src/obj/cortex_con/event_encoding.o src/obj/cortex_con/seq.o src/obj/cortex_con/element.o src/obj/cortex_con/hash_value.o src/obj/cortex_con/hash_table.o src/obj/cortex_con/dB_graph.o src/obj/cortex_con/file_reader.o src/obj/cortex_con/cortex_con.o

CORTEX_VAR_OBJ = src/obj/cortex_var/many_colours/cortex_var.o src/obj/cortex_var/many_colours/binary_kmer.o src/obj/cortex_var/many_colours/element.o src/obj/cortex_var/many_colours/seq.o src/obj/cortex_var/many_colours/hash_value.o src/obj/cortex_var/many_colours/hash_table.o src/obj/cortex_var/many_colours/dB_graph.o src/obj/cortex_var/many_colours/file_reader.o  src/obj/cortex_var/many_colours/dB_graph_population.o  src/obj/cortex_var/many_colours/db_variants.o src/obj/cortex_var/many_colours/cmd_line.o src/obj/cortex_var/many_colours/event_encoding.o src/obj/cortex_var/many_colours/graph_info.o src/obj/cortex_var/many_colours/db_differentiation.o

BASIC_TESTS_OBJ = src/obj/basic/binary_kmer.o src/obj/basic/seq.o src/obj/test/basic/test_binary_kmer.o src/obj/test/basic/test_seq.o src/obj/test/basic/run_basic_tests.o src/obj/basic/event_encoding.o

HASH_TABLE_TESTS_OBJ = src/obj/test/hash_table/run_hash_table_tests.o src/obj/cortex_var/many_colours/element.o src/obj/cortex_var/many_colours/hash_value.o src/obj/cortex_var/many_colours/hash_table.o src/obj/test/hash_table/test_hash.o src/obj/basic/binary_kmer.o  src/obj/basic/seq.o src/obj/basic/event_encoding.o

CORTEX_VAR_TESTS_OBJ = src/obj/test/cortex_var/many_colours/test_file_reader.o src/obj/test/cortex_var/many_colours/test_pop_load_and_print.o src/obj/test/cortex_var/many_colours/run_sv_trio_tests.o src/obj/test/cortex_var/many_colours/supernode_cmp.o src/obj/test/cortex_var/many_colours/test_pop_supernode_consensus.o src/obj/test/cortex_var/many_colours/test_pop_element.o src/obj/test/cortex_var/many_colours/test_dB_graph_population.o src/obj/cortex_var/many_colours/binary_kmer.o src/obj/cortex_var/many_colours/element.o src/obj/cortex_var/many_colours/seq.o src/obj/cortex_var/many_colours/hash_value.o src/obj/cortex_var/many_colours/hash_table.o src/obj/cortex_var/many_colours/dB_graph.o src/obj/cortex_var/many_colours/file_reader.o src/obj/cortex_var/many_colours/dB_graph_population.o src/obj/cortex_var/many_colours/db_variants.o src/obj/cortex_var/many_colours/event_encoding.o src/obj/cortex_var/many_colours/graph_info.o

CORTEX_VAR_CMD_LINE_TESTS_OBJ = src/obj/cortex_var/many_colours/cmd_line.o src/obj/test/cortex_var/many_colours/test_cmd_line.o src/obj/test/cortex_var/many_colours/run_cmd_line_tests.o src/obj/cortex_var/many_colours/binary_kmer.o src/obj/cortex_var/many_colours/element.o src/obj/cortex_var/many_colours/seq.o src/obj/cortex_var/many_colours/hash_value.o src/obj/cortex_var/many_colours/hash_table.o src/obj/cortex_var/many_colours/dB_graph.o src/obj/cortex_var/many_colours/file_reader.o src/obj/cortex_var/many_colours/dB_graph_population.o src/obj/cortex_var/many_colours/db_variants.o src/obj/cortex_var/many_colours/event_encoding.o

MAXK_AND_TEXT = $(join "", $(MAXK))
NUMCOLS_AND_TEST = $(join "_c", $(NUM_COLS))


cortex_con : remove_objects $(GRAPH_OBJ) 
	mkdir -p $(BIN); $(CC) -lm $(OPT) $(OPT_COLS) -o $(BIN)/cortex_con_$(MAXK) $(GRAPH_OBJ)

cortex_var : remove_objects $(CORTEX_VAR_OBJ)
	mkdir -p $(BIN); $(CC) $(CFLAGS_CORTEX_VAR) -lm $(OPT) $(OPT_COLS) -o $(BIN)/cortex_var_$(join $(MAXK_AND_TEXT),$(NUMCOLS_AND_TEST)) $(CORTEX_VAR_OBJ)

run_basic_tests : remove_objects $(BASIC_TESTS_OBJ)
	mkdir -p $(BIN); $(CC) $(OPT) $(CFLAGS_BASIC) $(CFLAGS_CUNIT) -o $(BIN)/run_basic_tests_$(MAXK) $(BASIC_TESTS_OBJ) -lcunit

run_hash_table_tests : remove_objects $(HASH_TABLE_TESTS_OBJ)
	mkdir -p $(BIN); $(CC) $(OPT) $(CFLAGS_CUNIT) $(CFLAGS_HASH_TABLE_TESTS) -o $(BIN)/run_hash_table_tests_$(MAXK) $(HASH_TABLE_TESTS_OBJ) -lcunit

run_cortex_var_cmdline_tests : remove_objects $(CORTEX_VAR_CMD_LINE_TESTS_OBJ)
	mkdir -p $(BIN); $(CC) $(CFLAGS_CORTEX_VAR_CMD_LINE_TESTS) $(OPT)  -o $(BIN)/run_cortex_var_cmdline_tests $(CORTEX_VAR_CMD_LINE_TESTS_OBJ) -lcunit

run_cortex_var_tests : remove_objects $(CORTEX_VAR_TESTS_OBJ)
	mkdir -p $(BIN); $(CC) $(CFLAGS_CORTEX_VAR_TESTS) $(OPT)  -o $(BIN)/run_cortex_var_tests_$(MAXK) $(CORTEX_VAR_TESTS_OBJ) -lcunit

.PHONY : clean
clean :
	rm -rf $(BIN)/*
	rm -rf src/obj

remove_objects:
	rm -rf src/obj/*







#pattern rules


src/obj/cortex_con/%.o : src/cortex_con/%.c include/cortex_con/%.h
	mkdir -p src/obj/cortex_con;  $(CC) $(CFLAGS_GRAPH) $(OPT) -c $< -o $@

src/obj/cortex_con/%.o : src/basic/%.c include/basic/%.h
	mkdir -p src/obj/cortex_con;  $(CC) $(CFLAGS_GRAPH) $(OPT) -c $< -o $@

src/obj/cortex_con/%.o : src/basic/event_encoding/base_encoding/%.c include/basic/event_encoding/base_encoding/%.h
	mkdir -p src/obj/cortex_con;  $(CC) $(CFLAGS_GRAPH) $(OPT) -c $< -o $@

src/obj/cortex_con/%.o : src/hash_table/hash_key/bob_jenkins/%.c include/hash_table/hash_value.h
	mkdir -p src/obj/cortex_con;  $(CC) $(CFLAGS_GRAPH) $(OPT) -c $< -o $@

src/obj/cortex_con/%.o : src/hash_table/open_hash/%.c include/hash_table/open_hash/%.h
	mkdir -p src/obj/cortex_con;  $(CC) $(CFLAGS_GRAPH) $(OPT) -c $< -o $@

src/obj/cortex_con/%.o : src/cortex_con/%.c 
	mkdir -p src/obj/cortex_con;  $(CC) $(CFLAGS_GRAPH) $(OPT) -c $? -o $@

src/obj/cortex_var/many_colours/%.o : src/cortex_var/many_colours/%.c include/cortex_var/many_colours/%.h
	mkdir -p src/obj/cortex_var/many_colours; $(CC) $(CFLAGS_CORTEX_VAR) $(OPT)  -c $< -o $@

src/obj/basic/%.o : src/basic/%.c  include/basic/%.h
	mkdir -p src/obj/basic/; $(CC) $(CFLAGS_BASIC) $(OPT) -c $< -o $@

src/obj/basic/%.o : src/basic/event_encoding/base_encoding/%.c  include/basic/event_encoding/base_encoding/%.h
	mkdir -p src/obj/basic/; $(CC) $(CFLAGS_BASIC) $(OPT) -c $< -o $@

src/obj/test/basic/%.o :  src/test/basic/%.c include/test/basic/%.h
	mkdir -p src/obj/test/basic; $(CC) $(CFLAGS_BASIC_TESTS) $(OPT) -c $< -o $@

src/obj/test/hash_table/open_hash/%.o : src/hash_table/open_hash/%.c include/hash_table/open_hash/%.h
	mkdir -p src/obj/test/hash_table/open_hash; $(CC) $(CFLAGS_HASH_TABLE_TESTS) $(OPT)  -c $< -o $@

src/obj/test/hash_table/hash_key/bob_jenkins/%.o : src/hash_table/hash_key/bob_jenkins/%.c include/hash_table/hash_key/bob_jenkins/%.h
	mkdir -p src/obj/test/hash_table/hash_key/bob_jenkins; $(CC) $(CFLAGS_HASH_TABLE_TESTS) $(OPT)  -c $< -o $@

src/obj/graph/hash_table/open_hash/%.o : src/hash_table/open_hash/%.c include/hash_table/open_hash/%.h
	mkdir -p src/obj/graph/hash_table/open_hash; $(CC) $(CFLAGS_GRAPH) $(OPT)  -c $< -o $@

src/obj/graph/hash_table/hash_key/bob_jenkins/%.o : src/hash_table/hash_key/bob_jenkins/%.c
	mkdir -p src/obj/graph/hash_table/hash_key/bob_jenkins; $(CC) $(CFLAGS_GRAPH) $(OPT)  -c $< -o $@

src/obj/test/hash_table/%.o : src/test/hash_table/%.c include/test/hash_table/%.h
	mkdir -p src/obj/test/hash_table; $(CC) $(CFLAGS_CORTEX_VAR) $(CFLAGS_HASH_TABLE_TESTS) $(OPT) -c $< -o $@


src/obj/cortex_var/many_colours/%.o : src/basic/event_encoding/base_encoding/%.c include/basic/event_encoding/base_encoding/%.h
	mkdir -p src/obj/cortex_var/many_colours;  $(CC) $(CFLAGS_CORTEX_VAR) $(OPT) -c $< -o $@

src/obj/cortex_var/many_colours/hash_table/open_hash/%.o :  src/hash_table/open_hash/%.c include/hash_table/open_hash/%.h
	mkdir -p src/obj/cortex_var/many_colours/hash_table/open_hash; $(CC) $(CFLAGS_CORTEX_VAR) $(OPT)  -c $< -o $@

src/obj/cortex_var/many_colours/hash_table/hash_key/bob_jenkins/%.o :  src/hash_table/hash_key/bob_jenkins/%.c include/hash_table/hash_key/bob_jenkins/%.h
	mkdir -p src/obj/cortex_var/many_colours/hash_table/hash_key/bob_jenkins; $(CC) $(CFLAGS_CORTEX_VAR) $(OPT)  -c $< -o $@


src/obj/cortex_var/many_colours/%.o :  src/basic/%.c include/basic/%.h
	mkdir -p src/obj/cortex_var/many_colours; $(CC) $(CFLAGS_CORTEX_VAR) $(OPT)  -c $< -o $@

src/obj/cortex_var/many_colours/%.o :  src/hash_table/hash_key/bob_jenkins/%.c include/hash_table/%.h
	mkdir -p src/obj/cortex_var/many_colours; $(CC) $(CFLAGS_CORTEX_VAR) $(OPT)  -c $< -o $@

src/obj/cortex_var/many_colours/%.o :  src/hash_table/open_hash/%.c include/hash_table/open_hash/%.h
	mkdir -p src/obj/cortex_var/many_colours; $(CC) $(CFLAGS_CORTEX_VAR) $(OPT)  -c $< -o $@

src/obj/cortex_var/many_colours/%.o :  src/cortex_var/core/%.c include/cortex_var/core/%.h
	mkdir -p src/obj/cortex_var/many_colours; $(CC) $(CFLAGS_CORTEX_VAR) $(OPT)  -c $< -o $@


src/obj/test/cortex_var/many_colours/%.o :  src/test/cortex_var/many_colours/%.c include/test/cortex_var/many_colours/%.h
	mkdir -p src/obj/test/cortex_var/many_colours; $(CC)  $(CFLAGS_CORTEX_VAR_TESTS) $(OPT)  -c $< -o $@

