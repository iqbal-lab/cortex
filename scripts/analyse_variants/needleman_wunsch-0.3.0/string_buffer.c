/*
  string_buffer.c
  project: string_buffer
  author: Isaac Turner <isaac.turner@dtc.ox.ac.uk>
  Copyright (C) 09-May-2011
 
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

#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) > (y) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h> // toupper() and tolower()

#include <zlib.h>

#include "string_buffer.h"


STRING_BUFFER* string_buff_init(const t_buf_pos size)
{
  STRING_BUFFER* sbuf = (STRING_BUFFER*) malloc(sizeof(STRING_BUFFER));

  sbuf->buff = (char*) malloc(size);
  sbuf->len = 0;
  sbuf->size = size;

  if(sbuf->buff == NULL)
  {
    fprintf(stderr, "Error: STRING_BUFFER couldn't be created with %lui bytes.",
            size);
    exit(-1);
  }

  return sbuf;
}

STRING_BUFFER* string_buff_create(const char* str)
{
  t_buf_pos str_len = strlen(str);

  STRING_BUFFER* sbuf = string_buff_init(str_len+1);

  strcpy(sbuf->buff, str);
  sbuf->buff[str_len] = '\0';

  sbuf->len = str_len;

  return sbuf;
}

void string_buff_reset(STRING_BUFFER* sbuf)
{
  sbuf->len = 0;
  
  if(sbuf->size > 0)
  {
    sbuf->buff[0] = '\0';
  }
}

void string_buff_free(STRING_BUFFER* sbuf)
{
  free(sbuf->buff);
  free(sbuf);
}

void string_buff_grow(STRING_BUFFER *sbuf, const t_buf_pos new_size)
{
  sbuf->size = new_size;
  char *new_buff = realloc(sbuf->buff, sbuf->size);
  
  if(new_buff == NULL)
  {
    fprintf(stderr, "Error: STRING_BUFFER couldn't be given more memory.  "
                    "Requested %lui bytes.  STRING_BUFFER begins '%s...'",
            sbuf->size, string_buff_substr(sbuf, 0, 5));
    
    free(sbuf->buff);
    
    exit(-1);
  }
  else {
    sbuf->buff = new_buff;
  }
}

void string_buff_shrink(STRING_BUFFER *sbuf, const t_buf_pos new_len)
{
  sbuf->len = new_len;
  sbuf->buff[new_len] = '\0';
}

STRING_BUFFER* string_buff_copy(STRING_BUFFER* sbuf)
{
  // One byte for the string end / null char \0
  STRING_BUFFER* sbuf_cpy = string_buff_init(sbuf->len+1);
  
  strcpy(sbuf_cpy->buff, sbuf->buff);
  sbuf_cpy->buff[sbuf->len] = '\0';

  sbuf_cpy->len = sbuf->len;
  
  return sbuf_cpy;
}


/********************/
/* String functions */
/********************/

void string_buff_add(STRING_BUFFER* sbuf, const char* txt)
{
  int str_len = strlen(txt);
  string_buff_addn(sbuf, txt, str_len);
}

void string_buff_addn(STRING_BUFFER* sbuf, const char* txt, const t_buf_pos len)
{
  // plus 1 for '\0'
  if(sbuf->len + len + 1 > sbuf->size)
  {
    string_buff_grow(sbuf, sbuf->len + len + 1);
  }

  strncpy(sbuf->buff+sbuf->len, txt, len);
  sbuf->len += len;

  sbuf->buff[sbuf->len] = '\0';
}

void string_buff_chomp(STRING_BUFFER *sbuf)
{
  while(sbuf->len >= 1)
  {
    char last_char = sbuf->buff[sbuf->len-1];

    if(last_char == '\n' || last_char == '\r')
    {
      sbuf->buff[--(sbuf->len)] = '\0';
    }
    else
    {
      break;
    }
  }
}

char* string_buff_substr(STRING_BUFFER *sbuf, const t_buf_pos start,
                         const t_buf_pos len)
{
  char* mem_start = sbuf->buff + start;
  t_buf_pos mem_length = MIN(len, sbuf->buff + sbuf->len - mem_start);

  if (mem_length < 0)
  {
    fprintf(stderr, "string_buff_substr(%lui, %lui) end < start on "
                    "STRING_BUFFER with length %lui\n",
            start, len, sbuf->len);
    exit(-1);
  }

  char* new_string = (char*) malloc(mem_length+1);
  strncpy(new_string, mem_start, mem_length);
  new_string[mem_length] = '\0';

  return new_string;
}

void string_buff_to_uppercase(STRING_BUFFER *sbuf)
{
  char* pos;
  char* end = sbuf->buff + sbuf->len;

  for(pos = sbuf->buff; pos < end; pos++)
  {
    *pos = toupper(*pos);
  }
}

void string_buff_to_lowercase(STRING_BUFFER *sbuf)
{
  char* pos;
  char* end = sbuf->buff + sbuf->len;

  for(pos = sbuf->buff; pos < end; pos++)
  {
    *pos = tolower(*pos);
  }
}


/*****************/
/* File handling */
/*****************/

// returns number of characters read
// or -1 if EOF
t_buf_pos string_buff_reset_readline(STRING_BUFFER *sbuf, gzFile *gz_file)
{
  string_buff_reset(sbuf);
  return string_buff_readline(sbuf, gz_file);
}

// returns number of characters read
// or -1 if EOF
t_buf_pos string_buff_readline(STRING_BUFFER *sbuf, gzFile *gz_file)
{
  t_buf_pos init_str_len = sbuf->len;

  // Enlarge *str allocated mem if needed
  // Need AT LEAST 2 free spaces - one for a character and one for \0
  if(sbuf->len+2 >= sbuf->size) {
    // Double buffer size
    string_buff_grow(sbuf, 2*sbuf->size);
  }

  // max characters to read = sbuf.size - sbuf.len
  char* read;
  while((read = gzgets(gz_file,
                       (char*)(sbuf->buff + sbuf->len),
                       sbuf->size - sbuf->len)) != NULL)
  {
    // Check if we hit the end of the line
    t_buf_pos read_length = strlen(sbuf->buff + sbuf->len);
    char* last_char = (char*)(sbuf->buff + sbuf->len + read_length - 1);
    
    // Get the new length of the string buffer
    // count should include the return chars
    sbuf->len += read_length;
    
    if(*last_char == '\n')
    {
      // Return characters read
      return sbuf->len - init_str_len;
    }

    // Double buffer size
    string_buff_grow(sbuf, 2*sbuf->size);
  }

  int bases_read = sbuf->len - init_str_len;

  return (bases_read > 0 ? sbuf->len - init_str_len : -1);
}

/**************************/
/* Other String Functions */
/**************************/
long split_str(const char* split, const char* txt, char*** result)
{
  long split_len = strlen(split);
  long txt_len = strlen(txt);

  // result is temporarily held here
  char** arr;

  if(split_len == 0)
  {
    // Special case
    if(txt_len == 0)
    {
      *result = NULL;
      return 0;
    }
    else
    {
      arr = (char**) malloc(txt_len * sizeof(char*));
    
      int i;
      for(i = 0; i < txt_len; i++)
      {
        arr[i] = (char*) malloc(2 * sizeof(char));
        arr[i][0] = txt[i];
        arr[i][1] = '\0';
      }

      *result = arr;
      return txt_len;
    }
  }
  
  const char* find = txt;
  long count = 1; // must have at least one item

  while((find = strstr(find, split)) != NULL)
  {
    //printf("Found1: '%s'\n", find);
    count++;
    find += split_len;
  }

  // Create return array
  arr = (char**) malloc(count * sizeof(char*));
  
  count = 0;
  const char* last_position = txt;

  while((find = strstr(last_position, split)) != NULL)
  {
    long str_len = find - last_position;

    arr[count] = (char*) malloc((str_len+1) * sizeof(char));
    strncpy(arr[count], last_position, str_len);
    arr[count][str_len] = '\0';
    
    count++;
    last_position = find + split_len;
  }

  // Copy last item
  long str_len = txt + txt_len - last_position;
  arr[count] = (char*) malloc((str_len+1) * sizeof(char));

  if(count == 0)
  {
    strcpy(arr[count], txt);
  }
  else
  {
    strncpy(arr[count], last_position, str_len);
  }

  arr[count][str_len] = '\0';
  count++;
  
  *result = arr;
  
  return count;
}
