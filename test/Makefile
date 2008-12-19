include ../Makefile.inc


## You only need the -L when running on EBI network. On the laptop you can make install CUNit so it's in /usr/local/lib or whatever
## Note you also need to buils CUNit so installs to say /homes/zi/bin, and then export LD_LIBRARY_PATH=/homes/zi/bin/lib:$LD_LIBRARY_PATH

run_tests: run_tests.o test_binary_kmer.o test_dB_graph.o
	$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -L/homes/zi/bin/lib -o run_tests ../src/file_reader.o ../src/seq.o ../src/binary_kmer.o ../src/dB_graph_node.o  ../src/dB_graph.o  ../src/priority_queue.o test_binary_kmer.o test_dB_graph.o run_tests.o -lcunit 

test_binary_kmer.o: test_binary_kmer.c
	$(CC) $(CFLAGS) -c $? -o $@

test_dB_graph.o: test_dB_graph.c
	$(CC) $(CFLAGS) -c $? -o $@


run_tests.o: run_tests.c 
	$(CC) $(CFLAGS) -c $? -o $@


clean:
	-rm *.o 
	rm run_tests


