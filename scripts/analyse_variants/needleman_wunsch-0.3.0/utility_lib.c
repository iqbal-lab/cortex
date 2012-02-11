/*
 utility_lib.c
 project: utility library
 author: Isaac Turner <isaac.turner@dtc.ox.ac.uk>
 Copyright (C) 7-Nov-2011
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <sys/select.h>

#include "utility_lib.h"

/* Utility functions */
long parse_int(char* c, char* err)
{
  // atoi less dependable that newer strtol (it has no error response!)
  char* strtol_last_char_ptr = c;
  long value = strtol(c, &strtol_last_char_ptr, 10);
  
  // If pointer to end of number string hasn't moved -> error
  if(strtol_last_char_ptr == c)
  {
    fprintf(stderr, err, c);
    exit(-1);
  }
  
  return value;
}

char parse_entire_int(char *str, int *result)
{
  int len = strlen(str);

  char *strtol_last_char_ptr = str;
  *result = strtol(str, &strtol_last_char_ptr, 10);

  return (strtol_last_char_ptr == str+len);
}

int int_cmp(const void *aa, const void *bb)
{
  const int *a = aa, *b = bb;
  return (*a < *b) ? -1 : (*a > *b);
}


// Checks if anything is piping in
int stdin_is_ready()
{
  fd_set fdset;
  struct timeval timeout;

  FD_ZERO(&fdset);
  FD_SET(fileno(stdin), &fdset);
  timeout.tv_sec = 0;
  timeout.tv_usec = 1; // was 0
    
  return select(1, &fdset, NULL, NULL, &timeout) == 1 ? 1 : 0;
}
