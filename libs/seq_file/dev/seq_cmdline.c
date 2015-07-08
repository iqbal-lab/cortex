#include "seq_file.h"

#define print_prompt() do { printf("#"); fflush(stdout); } while(0)

int main(int argc, char **argv)
{
  (void)argv;

  if(argc != 1) {
    fprintf(stderr, "usage: seq_cmdline\n");
    exit(EXIT_FAILURE);
  }

  print_prompt();
  seq_file_t *file = seq_dopen(fileno(stdin), 0, 0, 0);

  if(file == NULL)
  {
    fprintf(stderr, "Error: couldn't open stdin\n");
    exit(EXIT_FAILURE);
  }

  read_t *read = seq_read_new();
  
  while(seq_read(file, read) > 0)
  {
    printf("%s\n", read->seq.b);
    print_prompt();
  }

  seq_close(file);
  seq_read_free(read);

  return EXIT_SUCCESS;
}
