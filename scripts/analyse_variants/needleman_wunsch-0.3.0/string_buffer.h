/*
  string_buffer.h
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

#ifndef STRING_BUFFER_FILE_SEEN
#define STRING_BUFFER_FILE_SEEN

#define t_buf_pos unsigned long

// needed for FILE*
//#include <stdio.h>
// needed for gzFile
#include <zlib.h>

typedef struct STRING_BUFFER STRING_BUFFER;

struct STRING_BUFFER {
  char *buff;
  t_buf_pos len; // length of the string
  t_buf_pos size; // buffer size
};

// Creation, reset, free and memory expansion
STRING_BUFFER* string_buff_init(const t_buf_pos size);
STRING_BUFFER* string_buff_create(const char* str);
void string_buff_reset(STRING_BUFFER* sbuf);
void string_buff_free(STRING_BUFFER* sbuf);
void string_buff_grow(STRING_BUFFER *sbuf, const t_buf_pos new_size);
void string_buff_shrink(STRING_BUFFER *sbuf, const t_buf_pos new_len);
STRING_BUFFER* string_buff_copy(STRING_BUFFER* sbuf);

// Useful String functions
void string_buff_add(STRING_BUFFER* sbuf, const char* txt);
void string_buff_addn(STRING_BUFFER* sbuf, const char* txt, const t_buf_pos len);
void string_buff_chomp(STRING_BUFFER *sbuf);
char* string_buff_substr(STRING_BUFFER *sbuf, const t_buf_pos start, const t_buf_pos len);
void string_buff_to_uppercase(STRING_BUFFER *sbuf);
void string_buff_to_lowercase(STRING_BUFFER *sbuf);

// Reading a file
t_buf_pos string_buff_reset_readline(STRING_BUFFER *sbuf, gzFile *gz_file);
t_buf_pos string_buff_readline(STRING_BUFFER *sbuf, gzFile *gz_file);

// Other String functions
long split_str(const char* split, const char* txt, char*** result);

#endif
