/*
 stream_buffer.h
 project: string_buffer
 url: https://github.com/noporpoise/StringBuffer
 author: Isaac Turner <turner.isaac@gmail.com>
 license: Public Domain
 Jan 2015
*/

#ifndef STRING_BUFFER_FILE_SEEN
#define STRING_BUFFER_FILE_SEEN

#include <stdio.h> // needed for FILE
#include <zlib.h> // needed for gzFile
#include <stdarg.h> // needed for va_list

#include "stream_buffer.h"


typedef struct
{
  char *b;
  size_t end; // end is index of \0
  size_t size; // size should be >= end+1 to allow for \0
} StrBuf;


//
// Creation, reset, free and memory expansion
//

// Constructors / Destructors
StrBuf* strbuf_new(size_t len);
static inline void strbuf_free(StrBuf *sb) { free(sb->b); free(sb); }

// Place a string buffer into existing memory. Example:
//   StrBuf buf;
//   strbuf_alloc(&buf, 100);
//   ...
//   strbuf_dealloc(&buf);
StrBuf* strbuf_alloc(StrBuf *sb, size_t len);

static inline void strbuf_dealloc(StrBuf *sb) {
  free(sb->b);
  memset(sb, 0, sizeof(*sb));
}

// Copy a string or existing string buffer
StrBuf* strbuf_create(const char *str);
StrBuf* strbuf_clone(const StrBuf *sb);

// Clear the content of an existing StrBuf (sets size to 0)
static inline void strbuf_reset(StrBuf *sb) {
  if(sb->b) { sb->b[sb->end = 0] = '\0'; }
}

//
// Resizing
//

// Ensure capacity for len characters plus '\0' character - exits on FAILURE
static inline void strbuf_ensure_capacity(StrBuf *sb, size_t len) {
  cbuf_capacity(&sb->b, &sb->size, len);
}

// Same as above, but update pointer if it pointed to resized array
void strbuf_ensure_capacity_update_ptr(StrBuf *sbuf, size_t size, const char **ptr);

// Resize the buffer to have capacity to hold a string of length new_len
// (+ a null terminating character).  Can also be used to downsize the buffer's
// memory usage.  Returns 1 on success, 0 on failure.
char strbuf_resize(StrBuf *sb, size_t new_len);

//
// Useful String functions
//

#define strbuf_len(sb)  ((sb)->end)

#define strbuf_char(sb,idx) (sb)->b[idx]

// Note: in MACROs we use local variables to avoid multiple evaluation of param

// Set string buffer to contain a given string
#define strbuf_set(__sb,__str) do         \
{                                         \
  StrBuf     *_sb  = (__sb);              \
  const char *_str = (__str);             \
  size_t _s = strlen(_str);               \
  strbuf_ensure_capacity(_sb,_s);         \
  memcpy(_sb->b, _str, _s);               \
  _sb->b[_sb->end = _s] = '\0';           \
} while(0)

// Set string buffer to match existing string buffer
#define strbuf_set_buff(__dst,__src) do      \
{                                            \
  StrBuf       *_dst = (__dst);              \
  const StrBuf *_src = (__src);              \
  strbuf_ensure_capacity(_dst, _src->end);   \
  memmove(_dst->b, _src->b, _src->end);      \
  _dst->b[_dst->end = _src->end] = '\0';     \
} while(0)

// Add a character to the end of this StrBuf
#define strbuf_append_char(__sb,__c) do     \
{                                           \
  StrBuf *_sb = (__sb);                     \
  char    _c  = (__c);                      \
  strbuf_ensure_capacity(_sb, _sb->end+1);  \
  _sb->b[_sb->end] = _c;                    \
  _sb->b[++_sb->end] = '\0';                \
} while(0)

// Copy N characters from a character array to the end of this StrBuf
// strlen(__str) must be >= __n
#define strbuf_append_strn(__sb,__str,__n) do                   \
{                                                               \
  StrBuf     *_sb  = (__sb);                                    \
  const char *_str = (__str);                                   \
  size_t      _n   = (__n);                                     \
  strbuf_ensure_capacity_update_ptr(_sb, _sb->end+_n, &_str);   \
  memcpy(_sb->b+_sb->end, _str, _n);                            \
  _sb->b[_sb->end = _sb->end+_n] = '\0';                        \
} while(0)

// Copy a character array to the end of this StrBuf
// name char* _str2 since strbuf_append_strn uses _str
#define strbuf_append_str(__sb,__str) do {        \
  const char *_str2 = (__str);                    \
  strbuf_append_strn(__sb, _str2, strlen(_str2)); \
} while(0)

#define strbuf_append_buff(__sb1,__sb2) do {     \
  const StrBuf *_sb2 = (__sb2);                  \
  strbuf_append_strn(__sb1, _sb2->b, _sb2->end); \
} while(0)

// Convert integers to string to append
void strbuf_append_int(StrBuf *buf, int value);
void strbuf_append_long(StrBuf *buf, long value);
void strbuf_append_ulong(StrBuf *buf, unsigned long value);

// Append a given string in lower or uppercase
void strbuf_append_strn_lc(StrBuf *buf, const char *str, size_t len);
void strbuf_append_strn_uc(StrBuf *buf, const char *str, size_t len);

// Append char `c` `n` times
void strbuf_append_charn(StrBuf *buf, char c, size_t n);

#define strbuf_shrink(__sb,__len) do {                    \
  StrBuf *_sb = (__sb); _sb->b[_sb->end = (__len)] = 0;   \
} while(0)

#define strbuf_dup_str(sb) strdup((sb)->b)

// Remove \r and \n characters from the end of this StrBuf
// Returns the number of characters removed
size_t strbuf_chomp(StrBuf *sb);

// Reverse a string
void strbuf_reverse(StrBuf *sb);

// Get a substring as a new null terminated char array
// (remember to free the returned char* after you're done with it!)
char* strbuf_substr(const StrBuf *sb, size_t start, size_t len);

// Change to upper or lower case
void strbuf_to_uppercase(StrBuf *sb);
void strbuf_to_lowercase(StrBuf *sb);

// Copy a string to this StrBuf, overwriting any existing characters
// Note: dst_pos + len can be longer the the current dst StrBuf
void strbuf_copy(StrBuf *dst, size_t dst_pos,
                 const char *src, size_t len);

// Insert: copy to a StrBuf, shifting any existing characters along
void strbuf_insert(StrBuf *dst, size_t dst_pos,
                   const char *src, size_t len);

// Overwrite dst_pos..(dst_pos+dst_len-1) with src_len chars from src
// if dst_len != src_len, content to the right of dst_len is shifted
// Example:
//   strbuf_set(sb, "aaabbccc");
//   char *data = "xxx";
//   strbuf_overwrite(sb,3,2,data,strlen(data));
//   // sb is now "aaaxxxccc"
//   strbuf_overwrite(sb,3,2,"_",1);
//   // sb is now "aaa_ccc"
void strbuf_overwrite(StrBuf *dst, size_t dst_pos, size_t dst_len,
                      const char *src, size_t src_len);

// Remove characters from the buffer
//   strbuf_set(sb, "aaaBBccc");
//   strbuf_delete(sb, 3, 2);
//   // sb is now "aaaccc"
void strbuf_delete(StrBuf *sb, size_t pos, size_t len);

//
// sprintf
//

// sprintf to the end of a StrBuf (adds string terminator after sprint)
int strbuf_sprintf(StrBuf *sb, const char *fmt, ...)
  __attribute__ ((format(printf, 2, 3)));

// Print at a given position (overwrite chars at positions >= pos)
int strbuf_sprintf_at(StrBuf *sb, size_t pos, const char *fmt, ...)
  __attribute__ ((format(printf, 3, 4)));

int strbuf_vsprintf(StrBuf *sb, size_t pos, const char *fmt, va_list argptr)
  __attribute__ ((format(printf, 3, 0)));

// sprintf without terminating character
// Does not prematurely end the string if you sprintf within the string
// (terminates string if sprintf to the end)
int strbuf_sprintf_noterm(StrBuf *sb, size_t pos, const char *fmt, ...)
  __attribute__ ((format(printf, 3, 4)));

//
// Reading files
//

// Reading a FILE
size_t strbuf_reset_readline(StrBuf *sb, FILE *file);
size_t strbuf_readline(StrBuf *sb, FILE *file);
size_t strbuf_skipline(FILE *file);
size_t strbuf_readline_buf(StrBuf *sb, FILE *file, StreamBuffer *in);
size_t strbuf_skipline_buf(FILE* file, StreamBuffer *in);
size_t strbuf_read(StrBuf *sb, FILE *file, size_t len);

// Reading a gzFile
size_t strbuf_reset_gzreadline(StrBuf *sb, gzFile gz_file);
size_t strbuf_gzreadline(StrBuf *sb, gzFile gz_file);
size_t strbuf_gzskipline(gzFile gz_file);
size_t strbuf_gzreadline_buf(StrBuf *sb, gzFile gz_file, StreamBuffer *in);
size_t strbuf_gzskipline_buf(gzFile file, StreamBuffer *in);
size_t strbuf_gzread(StrBuf *sb, gzFile gz_file, size_t len);

// Read a line that has at least one character that is not \r or \n
// these functions do not call reset before reading
// Returns the number of characters read
size_t strbuf_readline_nonempty(StrBuf *line, FILE *fh);
size_t strbuf_gzreadline_nonempty(StrBuf *line, gzFile gz);

//
// String functions
//

// Trim whitespace characters from the start and end of a string
void strbuf_trim(StrBuf *sb);

// Trim the characters listed in `list` from the left of `sb`
// `list` is a null-terminated string of characters
void strbuf_ltrim(StrBuf *sb, const char *list);

// Trim the characters listed in `list` from the right of `sb`
// `list` is a null-terminated string of characters
void strbuf_rtrim(StrBuf *sb, const char *list);

/**************************/
/* Other String functions */
/**************************/

// `n` is the maximum number of bytes to copy including the NULL byte
// copies at most n bytes from `src` to `dst`
// Always appends a NULL terminating byte, unless n is zero.
// Returns a pointer to dst
char* string_safe_ncpy(char *dst, const char *src, size_t n);

// Replaces `sep` with \0 in str
// Returns number of occurances of `sep` character in `str`
// Stores `nptrs` pointers in `ptrs`
size_t string_split_str(char *str, char sep, char **ptrs, size_t nptrs);

// Replace one char with another in a string. Return number of replacements made
size_t string_char_replace(char *str, char from, char to);

void string_reverse_region(char *str, size_t length);
char string_is_all_whitespace(const char *s);
char* string_next_nonwhitespace(char *s);
char* string_trim(char *str);
// Chomp a string, returns new length
size_t string_chomp(char *str, size_t len);
size_t string_count_char(const char *str, char c);
size_t string_split(const char *split, const char *txt, char ***result);

#endif /* STRING_BUFFER_FILE_SEEN */
