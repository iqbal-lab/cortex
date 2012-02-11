ifdef DEBUG
	FLAGS="-DDEBUG=1"
endif

all:
	gcc -o needleman_wunsch $(FLAGS) -Wall -lz nw_cmdline.c needleman_wunsch.c utility_lib.c string_buffer.c

clean:
	rm needleman_wunsch
