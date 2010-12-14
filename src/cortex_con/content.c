/*
 * 
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * **********************************************************************
 *
 * This file is part of CORTEX.
 *
 * CORTEX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CORTEX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CORTEX.  If not, see <http://www.gnu.org/licenses/>.
 *
 * **********************************************************************
 */
#include <element.h>
#include <stdio.h>
#include <stdlib.h>
#include <binary_kmer.h>
#include <cmd_line.h>
#include <seq.h>
#include <limits.h>

int main(int argc, char **argv){
  CmdLine cmd_line = parse_cmdline(argc,argv,sizeof(Element));
  long long As=0;
  long long Cs=0;
  long long Gs=0;
  long long Ts=0;
  long long Us=0;

  int max_read_length = 1000;

  FILE* fp = fopen(cmd_line.input_filename, "r"); 
  if (fp == NULL){
    fprintf(stderr,"cannot open file:%s\n",cmd_line.input_filename);
    exit(1); //TODO - prefer to print warning and skip file and return an error code?
  }


  //----------------------------------
  // preallocate the memory used to read the sequences
  //----------------------------------
  Sequence * seq = malloc(sizeof(Sequence));
  if (seq == NULL){
    fputs("Out of memory trying to allocate Sequence\n",stderr);
    exit(1);
  }
  alloc_sequence(seq,max_read_length,LINE_MAX);
  
  int entry_length=0;
  boolean new_entry = true;
  boolean full_entry;
  do {

    switch (cmd_line.input_file_format) {
    case FASTQ:
      entry_length = read_sequence_from_fastq(fp,seq,max_read_length);
      break;
      
    case FASTA:
   
    entry_length = read_sequence_from_fasta(fp,seq,max_read_length,new_entry,&full_entry,0);
    new_entry = full_entry;
    break;
    }

    int i;
    for(i=0;i<entry_length;i++){
      //printf("index %i char %c\n",i,seq->seq[i]);

      switch (seq->seq[i])
	{
	  
	case 'A':
	  As++;
	  break;
	case 'C':
	  Cs++;
	  break;
	case 'G':
	  Gs++;
	  break;
	case 'T':
	  Ts++;
	  break;
	case 'a':
	  As++;
	  break;
	case 'c':
	  Cs++;
	  break;
	case 'g':
	  Gs++;
	  break;
	case 't':
	  Ts++;
	  break;
	default:
	  Us++;
	}
      
    }
    
  } while (entry_length>0);
  
  long long total = As + Cs + Gs + Ts + Us;

  printf("%qd As counted - %5.2f%%\n",As, (As*100.0)/total);
  printf("%qd Cs counted - %5.2f%%\n",Cs,(Cs*100.0)/total);
  printf("%qd Gs counted - %5.2f%%\n",Gs,(Gs*100.0)/total);
  printf("%qd Ts counted - %5.2f%%\n",Ts,(Ts*100.0)/total);
  printf("%qd Us counted - %5.2f%%\n",Us,(Us*100.0)/total);
  
  return 0;
}
 
