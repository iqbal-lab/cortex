/*
 utility_lib.c
 project: utility library
 author: Isaac Turner <turner.isaac@gmail.com>
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
#include <ctype.h> // isspace
#include <limits.h> // UINT_MAX etc.

#include <sys/select.h> // used in stdin_is_ready()

#include "utility_lib.h"

/* Utility functions */
long parse_int(char* c, char* err)
{
  // atoi less dependable than newer strtol (it has no error response!)
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

/* parse entire integer */

char parse_entire_int(const char *str, int *result)
{
  size_t len = strlen(str);

  char *strtol_last_char_ptr = (char*)str;
  long tmp = strtol(str, &strtol_last_char_ptr, 10);

  if(tmp > INT_MAX || tmp < INT_MIN || strtol_last_char_ptr != str+len)
  {
    return 0;
  }
  else
  {
    *result = (int)tmp;
    return 1;
  }
}

char parse_entire_uint(const char *str, unsigned int *result)
{
  size_t len = strlen(str);

  char *strtol_last_char_ptr = (char*)str;
  unsigned long tmp = strtoul(str, &strtol_last_char_ptr, 10);

  if(tmp > UINT_MAX || strtol_last_char_ptr != str+len)
  {
    return 0;
  }
  else
  {
    *result = (unsigned int)tmp;
    return 1;
  }
}

char parse_entire_long(const char *str, long *result)
{
  size_t len = strlen(str);

  char *strtol_last_char_ptr = (char*)str;
  long tmp = strtol(str, &strtol_last_char_ptr, 10);

  if(strtol_last_char_ptr == str+len)
  {
    *result = tmp;
    return 1;
  }
  else
  {
    return 0;
  }
}

char parse_entire_ulong(const char *str, unsigned long *result)
{
  size_t len = strlen(str);

  char *strtol_last_char_ptr = (char*)str;
  unsigned long tmp = strtoul(str, &strtol_last_char_ptr, 10);

  if(strtol_last_char_ptr == str+len)
  {
    *result = tmp;
    return 1;
  }
  else
  {
    return 0;
  }
}

char parse_entire_longlong(const char *str, long long *result)
{
  size_t len = strlen(str);

  char *strtol_last_char_ptr = (char*)str;
  long long tmp = strtoll(str, &strtol_last_char_ptr, 10);

  if(strtol_last_char_ptr == str+len)
  {
    *result = tmp;
    return 1;
  }
  else
  {
    return 0;
  }
}

char parse_entire_ulonglong(const char *str, unsigned long long *result)
{
  size_t len = strlen(str);

  char *strtol_last_char_ptr = (char*)str;
  unsigned long long tmp = strtoull(str, &strtol_last_char_ptr, 10);

  if(strtol_last_char_ptr == str+len)
  {
    *result = tmp;
    return 1;
  }
  else
  {
    return 0;
  }
}

// Compare integers

int int_cmp(const void *aa, const void *bb)
{
  const int *a = aa, *b = bb;
  return (*a < *b) ? -1 : (*a > *b);
}

// Convert an int to a string of 0 or 1 characters
char* int_to_binary(const int x)
{
  char *b = (char*) malloc(sizeof(char)*(32+1));
  char *p = b;

  unsigned int z;

  for(z = (1 << 31); z > 0; z >>= 1)
  {
    *p++ = (x & z) ? '1' : '0';
  }

  *p = '\0';

  return b;
}

// Convert an long to a string of 0 or 1 characters
char* long_to_binary(const long x)
{
  char *b = (char*) malloc(sizeof(char)*(64+1));
  char *p = b;

  unsigned long z;

  for(z = (1l << 63); z > 0; z >>= 1)
  {
    *p++ = (x & z) ? '1' : '0';
  }

  *p = '\0';

  return b;
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

char string_is_all_whitespace(const char* s)
{
  int i;

  for(i = 0; s[i] != '\0'; i++)
  {
    if(!isspace(s[i]))
    {
      return 0;
    }
  }

  return 1;
}

char* next_nonwhitespace(const char* s)
{
  while(*s != '\0')
  {
    if(!isspace(*s))
    {
      return (char*)s;
    }

    s++;
  }

  return NULL;
}

char* trim(char* str)
{
  while(isspace(*str))
  {
    str++;
  }

  size_t len = strlen(str);

  while(isspace(*(str+len-1)))
  {
    len--;
  }

  *(str+len) = '\0';

  return str;
}

long count_strchr(const char* str, const int c)
{
  int count = 0;
  const char *tmp = str;

  while((tmp = strchr(tmp,c)) != NULL)
  {
    tmp++;
    count++;
  }

  return count;
}
