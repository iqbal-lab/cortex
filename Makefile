#If on laptop mv for_use_in_laptop_Makefile.inc to Makefile.inc
#Similarly for in ome dir. UGLY. Sorry - open to suggestions.
include Makefile.inc


SRC_DIRS = src
TEST_DIRS = test
DIRS = src test

MARZAM_OBJ = src/main.o
SRC_OBJ = src/binary_kmer.o src/seq.o 
TEST_OBJ = test/test_binary_kmer.o
SRC_DIR = src
TEST_DIR = test

marzam:
	-for d in $(SRC_DIRS); do (cd $$d; $(MAKE) marzam); done

test: 
	-for d in $(SRC_DIRS); do (cd $$d; $(MAKE) marzam); done
	-for d in $(TEST_DIRS); do (cd $$d; $(MAKE) test); done

clean: 
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done





