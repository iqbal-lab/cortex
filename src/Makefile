ifndef CC
  CC = gcc	
endif

BIN = ../bin

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

#main program includes
IDIR_BASIC =../include/basic
IDIR_BASE_ENCODING = ${IDIR_BASIC}/event_encoding/base_encoding
IDIR_COLOUR_ENCODING = ${IDIR_BASIC}/event_encoding/solid_colour_encoding
IDIR_HASH  =../include/hash_table 
IDIR_GRAPH = include/cortex_con
IDIR_CORTEX_VAR =../include/pop_graph/sv_trio 
IDIR_CORTEX_VAR_CORE = ../include/pop_graph/core
IDIR_CORTEX_VAR_CMD_LINE =../include/pop_graph/sv_trio

## test code includes
IDIR_BASIC_TESTS =../include/test/basic
IDIR_HASH_TABLE_TESTS =../include/test/hash_table
IDIR_GRAPH_TESTS=../include/test/graph
IDIR_CORTEX_VAR_TESTS=../include/test/pop_graph/sv_trio

IDIR_CUNIT = /home/zam/dev/hg/CUnit/CUnit-2.1-0/CUnit/Headers


ifdef MAC
 MACFLAG    = -fnested-functions
 IDIR_CUNIT = /Users/zi/dev/hg/laptop_local/repos/CUnit/CUnit-2.1-0/CUnit/Headers
endif 


ARCH =  -m64 

ifdef 32_BITS
 ARCH =  
endif
OPT  		      = $(ARCH) -Wall  -O3  $(MACFLAG) -DNUMBER_OF_BITFIELDS_IN_BINARY_KMER=$(BITFIELDS)

ifdef DEBUG
 OPT                  = $(ARCH) -Wall -O0 -g $(MACFLAG) -DNUMBER_OF_BITFIELDS_IN_BINARY_KMER=$(BITFIELDS)
endif


CFLAGS_CUNIT          = -L/home/zam/bin/lib

CFLAGS_BASIC          = -I$(IDIR_BASIC) 
CFLAGS_GRAPH          = -I$(IDIR_BASIC) -I$(IDIR_HASH) -I$(IDIR_GRAPH) 
CFLAGS_SV_TRIO = -I$(IDIR_CORTEX_VAR_CORE) -I$(IDIR_BASIC)  -I$(IDIR_HASH)  -I$(IDIR_CORTEX_VAR) 


CFLAGS_BASIC_TESTS = -I$(IDIR_BASIC_TESTS)  -I$(IDIR_BASIC)  -I$(IDIR_CUNIT) -L/home/zam/bin/lib
CFLAGS_HASH_TABLE_TESTS =  -I$(IDIR_HASH) -I$(IDIR_HASH_TABLE_TESTS) -I$(IDIR_CUNIT)  -I$(IDIR_CORTEX_VAR) -L/home/zam/bin/lib
CFLAGS_GRAPH_TESTS = -I$(IDIR_GRAPH_TESTS)  -I$(IDIR_BASIC) -I$(IDIR_HASH) -I$(IDIR_GRAPH) -I$(IDIR_CUNIT) -L/home/zam/bin/lib 
CFLAGS_SV_TRIO_TESTS = -I$(IDIR_CORTEX_VAR_TESTS) -I$(IDIR_BASIC) -I$(IDIR_HASH) -I$(IDIR_CORTEX_VAR) -I$(IDIR_CORTEX_VAR_CORE) -I$(IDIR_CUNIT) -L/home/zam/bin/lib
CFLAGS_CORTEX_VAR_CMD_LINE_TESTS = -I$(IDIR_CORTEX_VAR_TESTS) -I$(IDIR_CORTEX_VAR_CMD_LINE) -I$(IDIR_BASIC) -I$(IDIR_CUNIT) -L/home/zam/bin/lib




SV_TRIO_OBJ = obj/pop_graph/sv_trio/graph.o obj/pop_graph/sv_trio/binary_kmer.o obj/pop_graph/sv_trio/element.o obj/pop_graph/sv_trio/seq.o obj/pop_graph/sv_trio/hash_value.o obj/pop_graph/sv_trio/hash_table.o obj/pop_graph/sv_trio/dB_graph.o obj/pop_graph/sv_trio/file_reader.o  obj/pop_graph/sv_trio/dB_graph_population.o obj/pop_graph/sv_trio/internal_oxford.o  obj/pop_graph/sv_trio/db_variants.o obj/pop_graph/sv_trio/cmd_line.o


BASIC_TESTS_OBJ = obj/basic/binary_kmer.o obj/basic/seq.o obj/test/basic/test_binary_kmer.o obj/test/basic/test_seq.o obj/test/basic/run_basic_tests.o

HASH_TABLE_TESTS_OBJ = obj/test/hash_table/run_hash_table_tests.o obj/pop_graph/sv_trio/element.o obj/pop_graph/sv_trio/hash_value.o obj/pop_graph/sv_trio/hash_table.o obj/test/hash_table/test_hash.o obj/basic/binary_kmer.o  obj/basic/seq.o 

SV_TRIO_TESTS_OBJ = obj/test/pop_graph/sv_trio/test_file_reader.o obj/test/pop_graph/sv_trio/test_pop_load_and_print.o obj/test/pop_graph/sv_trio/run_sv_trio_tests.o obj/test/pop_graph/sv_trio/supernode_cmp.o obj/test/pop_graph/sv_trio/test_pop_supernode_consensus.o obj/test/pop_graph/sv_trio/test_pop_element.o obj/test/pop_graph/sv_trio/test_dB_graph_population.o obj/pop_graph/sv_trio/binary_kmer.o obj/pop_graph/sv_trio/element.o obj/pop_graph/sv_trio/seq.o obj/pop_graph/sv_trio/hash_value.o obj/pop_graph/sv_trio/hash_table.o obj/pop_graph/sv_trio/dB_graph.o obj/pop_graph/sv_trio/file_reader.o obj/pop_graph/sv_trio/dB_graph_population.o obj/pop_graph/sv_trio/db_variants.o

CORTEX_VAR_CMD_LINE_TESTS_OBJ = obj/pop_graph/sv_trio/cmd_line.o obj/test/pop_graph/sv_trio/test_cmd_line.o obj/test/pop_graph/sv_trio/run_cmd_line_tests.o obj/pop_graph/sv_trio/binary_kmer.o obj/pop_graph/sv_trio/element.o obj/pop_graph/sv_trio/seq.o obj/pop_graph/sv_trio/hash_value.o obj/pop_graph/sv_trio/hash_table.o obj/pop_graph/sv_trio/dB_graph.o obj/pop_graph/sv_trio/file_reader.o obj/pop_graph/sv_trio/dB_graph_population.o obj/pop_graph/sv_trio/db_variants.o 


cortex_con : remove_objects $(GRAPH_OBJ) 
	mkdir -p $(BIN); $(CC) -lm $(OPT) -o $(BIN)/cortex_con_$(MAXK) $(GRAPH_OBJ)

cortex_var : remove_objects $(SV_TRIO_OBJ)
	mkdir -p $(BIN); $(CC) $(CFLAGS_SV_TRIO) -lm $(OPT) -o $(BIN)/cortex_var_$(MAXK) $(SV_TRIO_OBJ)

run_basic_tests : $(BASIC_TESTS_OBJ)
	mkdir -p $(BIN); $(CC) $(OPT) $(CFLAGS_CUNIT) -o $(BIN)/run_basic_tests_$(MAXK) $(BASIC_TESTS_OBJ) -lcunit

run_hash_table_tests : $(HASH_TABLE_TESTS_OBJ)
	mkdir -p $(BIN); $(CC) $(OPT) $(CFLAGS_CUNIT) $(CFLAGS_HASH_TABLE_TESTS) -o $(BIN)/run_hash_table_tests_$(MAXK) $(HASH_TABLE_TESTS_OBJ) -lcunit

run_cortex_var_cmdline_tests : $(CORTEX_VAR_CMD_LINE_TESTS_OBJ)
	mkdir -p $(BIN); $(CC) $(CFLAGS_CORTEX_VAR_CMD_LINE_TESTS) $(OPT)  -o $(BIN)/run_cortex_var_cmdline_tests $(CORTEX_VAR_CMD_LINE_TESTS_OBJ) -lcunit

.PHONY : clean
clean :
	rm -rf $(BIN)/*
	rm -rf obj

remove_objects:
	rm -rf obj/*


#run_sv_trio_tests : $(SV_TRIO_TESTS_OBJ)
#	mkdir -p $(BIN); $(CC) $(CFLAGS_SV_TRIO_TESTS) $(OPT)  -o $(BIN)/run_sv_trio_tests $(SV_TRIO_TESTS_OBJ) -lcunit





#pattern rules

obj/pop_graph/sv_trio/%.o : pop_graph/sv_trio/%.c ../include/pop_graph/sv_trio/%.h
	mkdir -p obj/pop_graph/sv_trio; $(CC) $(CFLAGS_SV_TRIO) $(OPT)  -c $< -o $@

obj/basic/%.o : basic/%.c  ../include/basic/%.h
	mkdir -p obj/basic/; $(CC) $(CFLAGS_BASIC) $(OPT) -c $< -o $@

obj/test/basic/%.o :  test/basic/%.c ../include/test/basic/%.h
	mkdir -p obj/test/basic; $(CC) $(CFLAGS_BASIC_TESTS) $(OPT) -c $< -o $@

obj/%.o : src/basic/event_encoding/base_encoding/%.c ../include/basic/event_encoding/base_encoding/%.h
	mkdir -p obj;  $(CC) $(CFLAGS_GRAPH) $(OPT) -c $< -o $@

obj/test/hash_table/open_hash/%.o : hash_table/open_hash/%.c ../include/hash_table/open_hash/%.h
	mkdir -p obj/test/hash_table/open_hash; $(CC) $(CFLAGS_HASH_TABLE_TESTS) $(OPT)  -c $< -o $@

obj/test/hash_table/hash_key/bob_jenkins/%.o : hash_table/hash_key/bob_jenkins/%.c ../include/hash_table/hash_key/bob_jenkins/%.h
	mkdir -p obj/test/hash_table/hash_key/bob_jenkins; $(CC) $(CFLAGS_HASH_TABLE_TESTS) $(OPT)  -c $< -o $@

#obj/test/hash_table/hash_key/crc32/%.o : hash_table/hash_key/crc32/%.c
#	mkdir -p obj/test/hash_table/hash_key/crc32; $(CC) $(CFLAGS_HASH_TABLE_TESTS) $(OPT) -c $< -o $@

obj/graph/hash_table/open_hash/%.o : hash_table/open_hash/%.c ../include/hash_table/open_hash/%.h
	mkdir -p obj/graph/hash_table/open_hash; $(CC) $(CFLAGS_GRAPH) $(OPT)  -c $< -o $@

obj/graph/hash_table/hash_key/bob_jenkins/%.o : hash_table/hash_key/bob_jenkins/%.c
	mkdir -p obj/graph/hash_table/hash_key/bob_jenkins; $(CC) $(CFLAGS_GRAPH) $(OPT)  -c $< -o $@

obj/test/hash_table/%.o : test/hash_table/%.c ../include/test/hash_table/%.h
	mkdir -p obj/test/hash_table; $(CC) $(CFLAGS_GRAPH) $(CFLAGS_HASH_TABLE_TESTS) $(OPT) -c $< -o $@

obj/pop_graph/sv_trio/hash_table/open_hash/%.o :  hash_table/open_hash/%.c ../include/hash_table/open_hash/%.h
	mkdir -p obj/pop_graph/sv_trio/hash_table/open_hash; $(CC) $(CFLAGS_GRAPH) $(OPT)  -c $< -o $@

obj/pop_graph/sv_trio/hash_table/hash_key/bob_jenkins/%.o :  hash_table/hash_key/bob_jenkins/%.c ../include/hash_table/hash_key/bob_jenkins/%.h
	mkdir -p obj/pop_graph/sv_trio/hash_table/hash_key/bob_jenkins; $(CC) $(CFLAGS_GRAPH) $(OPT)  -c $< -o $@

#obj/pop_graph/sv_trio/hash_table/hash_key/crc32/%.o :  hash_table/hash_key/crc32/%.c
#	mkdir -p obj/pop_graph/sv_trio/hash_table/hash_key/crc32; $(CC) $(CFLAGS_GRAPH) $(OPT)  -c $< -o $@

obj/pop_graph/sv_trio/%.o :  basic/%.c ../include/basic/%.h
	mkdir -p obj/pop_graph/sv_trio; $(CC) $(CFLAGS_SV_TRIO) $(OPT)  -c $< -o $@

obj/pop_graph/sv_trio/%.o :  hash_table/hash_key/bob_jenkins/%.c ../include/hash_table/%.h
	mkdir -p obj/pop_graph/sv_trio; $(CC) $(CFLAGS_SV_TRIO) $(OPT)  -c $< -o $@

obj/pop_graph/sv_trio/%.o :  hash_table/open_hash/%.c ../include/hash_table/open_hash/%.h
	mkdir -p obj/pop_graph/sv_trio; $(CC) $(CFLAGS_SV_TRIO) $(OPT)  -c $< -o $@

obj/pop_graph/sv_trio/%.o :  pop_graph/core/%.c ../include/pop_graph/core/%.h
	mkdir -p obj/pop_graph/sv_trio; $(CC) $(CFLAGS_SV_TRIO) $(OPT)  -c $< -o $@


obj/test/pop_graph/sv_trio/%.o :  test/pop_graph/sv_trio/%.c ../include/test/pop_graph/sv_trio/%.h
	mkdir -p obj/test/pop_graph/sv_trio; $(CC)  $(CFLAGS_SV_TRIO_TESTS) $(OPT)  -c $< -o $@

