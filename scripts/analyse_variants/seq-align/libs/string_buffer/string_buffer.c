/*
 string_buffer.c
 project: string_buffer
 url: https://github.com/noporpoise/StringBuffer
 author: Isaac Turner <turner.isaac@gmail.com>
 license: Public Domain
 Jan 2014
*/

// POSIX required for kill signal to work
#define _XOPEN_SOURCE 700

#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <signal.h> // kill on error
#include <ctype.h> // toupper() and tolower()

#include "string_buffer.h"

#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) > (y) ? (x) : (y))

#ifndef ROUNDUP2POW
  #define ROUNDUP2POW(x) (0x1UL << (64 - __builtin_clzl(x)))
#endif

#define exit_on_error() do { abort(); exit(EXIT_FAILURE); } while(0)

/*********************/
/*  Bounds checking  */
/*********************/

// Bounds check when inserting (pos <= len are valid)
#define _bounds_check_insert(sbuf,pos) \
        _call_bounds_check_insert(sbuf,pos,__FILE__,__LINE__,__func__)
#define _bounds_check_read_range(sbuf,start,len) \
        _call_bounds_check_read_range(sbuf,start,len,__FILE__,__LINE__,__func__)

static inline
void _call_bounds_check_insert(const StrBuf *sbuf, size_t pos,
                               const char *file, int line, const char *func)
{
  if(pos > sbuf->end)
  {
    fprintf(stderr, "%s:%i:%s() - out of bounds error "
                    "[index: %zu, num_of_bits: %zu]\n",
            file, line, func, pos, sbuf->end);
    errno = EDOM;
    exit_on_error();
  }
}

// Bounds check when reading a range (start+len < strlen is valid)
static inline
void _call_bounds_check_read_range(const StrBuf *sbuf, size_t start, size_t len,
                                   const char *file, int line, const char *func)
{
  if(start + len > sbuf->end)
  {
    fprintf(stderr, "%s:%i:%s() - out of bounds error "
                    "[start: %zu; length: %zu; strlen: %zu; buf:%.*s%s]\n",
            file, line, func, start, len, sbuf->end,
            (int)MIN(5, sbuf->end), sbuf->b, sbuf->end > 5 ? "..." : "");
    errno = EDOM;
    exit_on_error();
  }
}

/******************************/
/*  Constructors/Destructors  */
/******************************/

StrBuf* strbuf_new(size_t len)
{
  StrBuf *sbuf = calloc(1, sizeof(StrBuf));
  if(!sbuf) return NULL;
  if(!strbuf_alloc(sbuf, len)) { free(sbuf); return NULL; }
  return sbuf;
}

StrBuf* strbuf_alloc(StrBuf *sbuf, size_t len)
{
  sbuf->end  = 0;
  sbuf->size = ROUNDUP2POW(len+1);
  sbuf->b    = malloc(sbuf->size);
  if(!sbuf->b) return NULL;
  sbuf->b[0] = '\0';
  return sbuf;
}

StrBuf* strbuf_create(const char *str)
{
  size_t str_len = strlen(str);
  StrBuf *sbuf = strbuf_new(str_len+1);
  if(!sbuf) return NULL;
  memcpy(sbuf->b, str, str_len);
  sbuf->b[sbuf->end = str_len] = '\0';
  return sbuf;
}

StrBuf* strbuf_clone(const StrBuf *sbuf)
{
  // One byte for the string end / null char \0
  StrBuf *cpy = strbuf_new(sbuf->end+1);
  if(!cpy) return NULL;
  memcpy(cpy->b, sbuf->b, sbuf->end);
  cpy->b[cpy->end = sbuf->end] = '\0';
  return cpy;
}

/******************************/
/*  Resize Buffer Functions   */
/******************************/

// Resize the buffer to have capacity to hold a string of length new_len
// (+ a null terminating character).  Can also be used to downsize the buffer's
// memory usage.  Returns 1 on success, 0 on failure.
char strbuf_resize(StrBuf *sbuf, size_t new_len)
{
  size_t capacity = ROUNDUP2POW(new_len+1);
  char *new_buff = realloc(sbuf->b, capacity * sizeof(char));
  if(new_buff == NULL) return 0;

  sbuf->b    = new_buff;
  sbuf->size = capacity;

  if(sbuf->end > new_len)
  {
    // Buffer was shrunk - re-add null byte
    sbuf->end = new_len;
    sbuf->b[sbuf->end] = '\0';
  }

  return 1;
}

void strbuf_ensure_capacity_update_ptr(StrBuf *sbuf, size_t size, const char **ptr)
{
  if(sbuf->size <= size+1)
  {
    size_t oldcap = sbuf->size;
    char *oldbuf  = sbuf->b;

    if(!strbuf_resize(sbuf, size))
    {
      fprintf(stderr, "%s:%i:Error: _ensure_capacity_update_ptr couldn't resize "
                      "buffer. [requested %zu bytes; capacity: %zu bytes]\n",
              __FILE__, __LINE__, size, sbuf->size);
      exit_on_error();
    }

    // ptr may have pointed to sbuf, which has now moved
    if(*ptr >= oldbuf && *ptr < oldbuf + oldcap) {
      *ptr = sbuf->b + (*ptr - oldbuf);
    }
  }
}

/*******************/
/* Append Integers */
/*******************/

/*
 * Integer to string functions adapted from:
 *   https://www.facebook.com/notes/facebook-engineering/three-optimization-tips-for-c/10151361643253920
 */

#define P01 10
#define P02 100
#define P03 1000
#define P04 10000
#define P05 100000
#define P06 1000000
#define P07 10000000
#define P08 100000000
#define P09 1000000000
#define P10 10000000000
#define P11 100000000000
#define P12 1000000000000

/**
 * Return number of digits required to represent `num` in base 10.
 * Uses binary search to find number.
 * Examples:
 *   num_of_digits(0)   = 1
 *   num_of_digits(1)   = 1
 *   num_of_digits(10)  = 2
 *   num_of_digits(123) = 3
 */
static inline size_t num_of_digits(unsigned long v)
{
  if(v < P01) return 1;
  if(v < P02) return 2;
  if(v < P03) return 3;
  if(v < P12) {
    if(v < P08) {
      if(v < P06) {
        if(v < P04) return 4;
        return 5 + (v >= P05);
      }
      return 7 + (v >= P07);
    }
    if(v < P10) {
      return 9 + (v >= P09);
    }
    return 11 + (v >= P11);
  }
  return 12 + num_of_digits(v / P12);
}

void strbuf_append_ulong(StrBuf *buf, unsigned long value)
{
  // Append two digits at a time
  static const char digits[201] =
    "0001020304050607080910111213141516171819"
    "2021222324252627282930313233343536373839"
    "4041424344454647484950515253545556575859"
    "6061626364656667686970717273747576777879"
    "8081828384858687888990919293949596979899";

  size_t num_digits = num_of_digits(value);
  size_t pos = num_digits - 1;

  strbuf_ensure_capacity(buf, buf->end+num_digits);
  char *dst = buf->b + buf->end;

  while(value >= 100)
  {
    size_t v = value % 100;
    value /= 100;
    dst[pos] = digits[v * 2 + 1];
    dst[pos - 1] = digits[v * 2];
    pos -= 2;
  }

  // Handle last 1-2 digits
  if(value < 10) {
    dst[pos] = '0' + value;
  } else {
    dst[pos] = digits[value * 2 + 1];
    dst[pos - 1] = digits[value * 2];
  }

  buf->end += num_digits;
  buf->b[buf->end] = '\0';
}

void strbuf_append_int(StrBuf *buf, int value)
{
  // strbuf_sprintf(buf, "%i", value);
  strbuf_append_long(buf, value);
}

void strbuf_append_long(StrBuf *buf, long value)
{
  // strbuf_sprintf(buf, "%li", value);
  if(value < 0) { strbuf_append_char(buf, '-'); value = -value; }
  strbuf_append_ulong(buf, value);
}


/*
size_t strbuf_num_of_digits(unsigned long num)
{
  size_t digits = 1;
  while(1) {
    if(num < 10) return digits;
    if(num < 100) return digits+1;
    if(num < 1000) return digits+2;
    if(num < 10000) return digits+3;
    num /= 10000;
    digits += 4;
  }
  return digits;
}

void strbuf_append_ulong(StrBuf *buf, unsigned long value)
{
  // strbuf_sprintf(buf, "%lu", value);

  size_t i, num_digits = strbuf_num_of_digits(value);
  strbuf_ensure_capacity(buf, buf->end + num_digits);
  buf->end += num_digits;
  buf->b[buf->end] = '\0';

  for(i = 1; i <= num_digits; i++) {
    buf->b[buf->end - i] = '0' + (value % 10);
    value /= 10;
  }
}
*/

/********************/
/* Append functions */
/********************/

// Append string converted to lowercase
void strbuf_append_strn_lc(StrBuf *buf, const char *str, size_t len)
{
  strbuf_ensure_capacity(buf, buf->end + len);
  char *to = buf->b + buf->end;
  const char *end = str + len;
  for(; str < end; str++, to++) *to = tolower(*str);
  buf->end += len;
  buf->b[buf->end] = '\0';
}

// Append string converted to uppercase
void strbuf_append_strn_uc(StrBuf *buf, const char *str, size_t len)
{
  strbuf_ensure_capacity(buf, buf->end + len);
  char *to = buf->b + buf->end;
  const char *end = str + len;
  for(; str < end; str++, to++) *to = toupper(*str);
  buf->end += len;
  buf->b[buf->end] = '\0';
}

// Append char `c` `n` times
void strbuf_append_charn(StrBuf *buf, char c, size_t n)
{
  strbuf_ensure_capacity(buf, buf->end + n);
  memset(buf->b+buf->end, c, n);
  buf->end += n;
  buf->b[buf->end] = '\0';
}

// Remove \r and \n characters from the end of this StrBuf
// Returns the number of characters removed
size_t strbuf_chomp(StrBuf *sbuf)
{
  size_t old_len = sbuf->end;
  sbuf->end = string_chomp(sbuf->b, sbuf->end);
  return old_len - sbuf->end;
}

// Reverse a string
void strbuf_reverse(StrBuf *sbuf)
{
  string_reverse_region(sbuf->b, sbuf->end);
}

char* strbuf_substr(const StrBuf *sbuf, size_t start, size_t len)
{
  _bounds_check_read_range(sbuf, start, len);

  char *new_string = malloc((len+1) * sizeof(char));
  strncpy(new_string, sbuf->b + start, len);
  new_string[len] = '\0';

  return new_string;
}

void strbuf_to_uppercase(StrBuf *sbuf)
{
  char *pos, *end = sbuf->b + sbuf->end;
  for(pos = sbuf->b; pos < end; pos++) *pos = (char)toupper(*pos);
}

void strbuf_to_lowercase(StrBuf *sbuf)
{
  char *pos, *end = sbuf->b + sbuf->end;
  for(pos = sbuf->b; pos < end; pos++) *pos = (char)tolower(*pos);
}

// Copy a string to this StrBuf, overwriting any existing characters
// Note: dst_pos + len can be longer the the current dst StrBuf
void strbuf_copy(StrBuf *dst, size_t dst_pos, const char *src, size_t len)
{
  if(src == NULL || len == 0) return;

  _bounds_check_insert(dst, dst_pos);

  // Check if dst buffer can handle string
  // src may have pointed to dst, which has now moved
  size_t newlen = MAX(dst_pos + len, dst->end);
  strbuf_ensure_capacity_update_ptr(dst, newlen, &src);

  // memmove instead of strncpy, as it can handle overlapping regions
  memmove(dst->b+dst_pos, src, len * sizeof(char));

  if(dst_pos + len > dst->end)
  {
    // Extended string - add '\0' char
    dst->end = dst_pos + len;
    dst->b[dst->end] = '\0';
  }
}

// Insert: copy to a StrBuf, shifting any existing characters along
void strbuf_insert(StrBuf *dst, size_t dst_pos, const char *src, size_t len)
{
  if(src == NULL || len == 0) return;

  _bounds_check_insert(dst, dst_pos);

  // Check if dst buffer has capacity for inserted string plus \0
  // src may have pointed to dst, which will be moved in realloc when
  // calling ensure capacity
  strbuf_ensure_capacity_update_ptr(dst, dst->end + len, &src);

  char *insert = dst->b+dst_pos;

  // dst_pos could be at the end (== dst->end)
  if(dst_pos < dst->end)
  {
    // Shift some characters up
    memmove(insert + len, insert, (dst->end - dst_pos) * sizeof(char));

    if(src >= dst->b && src < dst->b + dst->size)
    {
      // src/dst strings point to the same string in memory
      if(src < insert) memmove(insert, src, len * sizeof(char));
      else if(src > insert) memmove(insert, src+len, len * sizeof(char));
    }
    else memmove(insert, src, len * sizeof(char));
  }
  else memmove(insert, src, len * sizeof(char));

  // Update size
  dst->end += len;
  dst->b[dst->end] = '\0';
}

// Overwrite dst_pos..(dst_pos+dst_len-1) with src_len chars from src
// if dst_len != src_len, content to the right of dst_len is shifted
// Example:
//   strbuf_set(sbuf, "aaabbccc");
//   char *data = "xxx";
//   strbuf_overwrite(sbuf,3,2,data,strlen(data));
//   // sbuf is now "aaaxxxccc"
//   strbuf_overwrite(sbuf,3,2,"_",1);
//   // sbuf is now "aaa_ccc"
void strbuf_overwrite(StrBuf *dst, size_t dst_pos, size_t dst_len,
                      const char *src, size_t src_len)
{
  _bounds_check_read_range(dst, dst_pos, dst_len);

  if(src == NULL) return;
  if(dst_len == src_len) strbuf_copy(dst, dst_pos, src, src_len);
  size_t newlen = dst->end + src_len - dst_len;

  strbuf_ensure_capacity_update_ptr(dst, newlen, &src);

  if(src >= dst->b && src < dst->b + dst->size)
  {
    if(src_len < dst_len) {
      // copy
      memmove(dst->b+dst_pos, src, src_len * sizeof(char));
      // resize (shrink)
      memmove(dst->b+dst_pos+src_len, dst->b+dst_pos+dst_len,
              (dst->end-dst_pos-dst_len) * sizeof(char));
    }
    else
    {
      // Buffer is going to grow and src points to this buffer

      // resize (grow)
      memmove(dst->b+dst_pos+src_len, dst->b+dst_pos+dst_len,
              (dst->end-dst_pos-dst_len) * sizeof(char));

      char *tgt = dst->b + dst_pos;
      char *end = dst->b + dst_pos + src_len;

      if(src < tgt + dst_len)
      {
        size_t len = MIN((size_t)(end - src), src_len);
        memmove(tgt, src, len);
        tgt += len;
        src += len;
        src_len -= len;
      }

      if(src >= tgt + dst_len)
      {
        // shift to account for resizing
        src += src_len - dst_len;
        memmove(tgt, src, src_len);
      }
    }
  }
  else
  {
    // resize
    memmove(dst->b+dst_pos+src_len, dst->b+dst_pos+dst_len,
            (dst->end-dst_pos-dst_len) * sizeof(char));
    // copy
    memcpy(dst->b+dst_pos, src, src_len * sizeof(char));
  }

  dst->end = newlen;
  dst->b[dst->end] = '\0';
}

void strbuf_delete(StrBuf *sbuf, size_t pos, size_t len)
{
  _bounds_check_read_range(sbuf, pos, len);
  memmove(sbuf->b+pos, sbuf->b+pos+len, sbuf->end-pos-len);
  sbuf->end -= len;
  sbuf->b[sbuf->end] = '\0';
}

/**************************/
/*         sprintf        */
/**************************/

int strbuf_vsprintf(StrBuf *sbuf, size_t pos, const char *fmt, va_list argptr)
{
  _bounds_check_insert(sbuf, pos);

  // Length of remaining buffer
  size_t buf_len = sbuf->size - pos;
  if(buf_len == 0 && !strbuf_resize(sbuf, sbuf->size << 1)) {
    fprintf(stderr, "%s:%i:Error: Out of memory\n", __FILE__, __LINE__);
    exit_on_error();
  }

  // Make a copy of the list of args incase we need to resize buff and try again
  va_list argptr_cpy;
  va_copy(argptr_cpy, argptr);

  int num_chars = vsnprintf(sbuf->b+pos, buf_len, fmt, argptr);
  va_end(argptr);

  // num_chars is the number of chars that would be written (not including '\0')
  // num_chars < 0 => failure
  if(num_chars < 0) {
    fprintf(stderr, "Warning: strbuf_sprintf something went wrong..\n");
    exit_on_error();
  }

  // num_chars does not include the null terminating byte
  if((size_t)num_chars+1 > buf_len)
  {
    strbuf_ensure_capacity(sbuf, pos+(size_t)num_chars);

    // now use the argptr copy we made earlier
    // Don't need to use vsnprintf now, vsprintf will do since we know it'll fit
    num_chars = vsprintf(sbuf->b+pos, fmt, argptr_cpy);
    if(num_chars < 0) {
      fprintf(stderr, "Warning: strbuf_sprintf something went wrong..\n");
      exit_on_error();
    }
  }
  va_end(argptr_cpy);

  // Don't need to NUL terminate, vsprintf/vnsprintf does that for us

  // Update length
  sbuf->end = pos + (size_t)num_chars;

  return num_chars;
}

// Appends sprintf
int strbuf_sprintf(StrBuf *sbuf, const char *fmt, ...)
{
  va_list argptr;
  va_start(argptr, fmt);
  int num_chars = strbuf_vsprintf(sbuf, sbuf->end, fmt, argptr);
  va_end(argptr);

  return num_chars;
}

int strbuf_sprintf_at(StrBuf *sbuf, size_t pos, const char *fmt, ...)
{
  _bounds_check_insert(sbuf, pos);

  va_list argptr;
  va_start(argptr, fmt);
  int num_chars = strbuf_vsprintf(sbuf, pos, fmt, argptr);
  va_end(argptr);

  return num_chars;
}

// Does not prematurely end the string if you sprintf within the string
// (vs at the end)
int strbuf_sprintf_noterm(StrBuf *sbuf, size_t pos, const char *fmt, ...)
{
  _bounds_check_insert(sbuf, pos);

  char last_char;
  size_t len = sbuf->end;

  // Call vsnprintf with NULL, 0 to get resulting string length without writing
  va_list argptr;
  va_start(argptr, fmt);
  int nchars = vsnprintf(NULL, 0, fmt, argptr);
  va_end(argptr);

  if(nchars < 0) {
    fprintf(stderr, "Warning: strbuf_sprintf something went wrong..\n");
    exit_on_error();
  }

  // Save overwritten char
  last_char = (pos+(size_t)nchars < sbuf->end) ? sbuf->b[pos+(size_t)nchars] : 0;

  va_start(argptr, fmt);
  nchars = strbuf_vsprintf(sbuf, pos, fmt, argptr);
  va_end(argptr);

  if(nchars < 0) {
    fprintf(stderr, "Warning: strbuf_sprintf something went wrong..\n");
    exit_on_error();
  }

  // Restore length if shrunk, null terminate if extended
  if(sbuf->end < len) sbuf->end = len;
  else sbuf->b[sbuf->end] = '\0';

  // Re-instate overwritten character
  sbuf->b[pos+(size_t)nchars] = last_char;

  return nchars;
}


/*****************/
/* File handling */
/*****************/

// Reading a FILE
size_t strbuf_readline(StrBuf *sbuf, FILE *file)
{
  return freadline(file, &(sbuf->b), &(sbuf->end), &(sbuf->size));
}

size_t strbuf_gzreadline(StrBuf *sbuf, gzFile file)
{
  return gzreadline(file, &sbuf->b, &sbuf->end, &sbuf->size);
}

// Reading a FILE
size_t strbuf_readline_buf(StrBuf *sbuf, FILE *file, StreamBuffer *in)
{
  return (size_t)freadline_buf(file, in, &sbuf->b, &sbuf->end, &sbuf->size);
}

size_t strbuf_gzreadline_buf(StrBuf *sbuf, gzFile file, StreamBuffer *in)
{
  return (size_t)gzreadline_buf(file, in, &sbuf->b, &sbuf->end, &sbuf->size);
}

size_t strbuf_skipline(FILE* file)
{
  return fskipline(file);
}

size_t strbuf_gzskipline(gzFile file)
{
  return gzskipline(file);
}

size_t strbuf_skipline_buf(FILE* file, StreamBuffer *in)
{
  return (size_t)fskipline_buf(file, in);
}

size_t strbuf_gzskipline_buf(gzFile file, StreamBuffer *in)
{
  return (size_t)gzskipline_buf(file, in);
}

#define _func_read_nonempty(name,type_t,__readline)                            \
  size_t name(StrBuf *line, type_t fh)                                         \
  {                                                                            \
    size_t i, origlen = line->end;                                             \
    while(__readline(line, fh) > 0) {                                          \
      i = origlen;                                                             \
      while(i < line->end && (line->b[i] == '\r' || line->b[i] == '\n'))       \
        i++;                                                                   \
      if(i < line->end) return line->end - origlen;                            \
      line->end = origlen;                                                     \
      line->b[line->end] = '\0';                                               \
    }                                                                          \
    return 0;                                                                  \
  }

_func_read_nonempty(strbuf_readline_nonempty,FILE*,strbuf_readline)
_func_read_nonempty(strbuf_gzreadline_nonempty,gzFile,strbuf_gzreadline)


#define _func_read(name,type_t,__read) \
  size_t name(StrBuf *sbuf, type_t file, size_t len)                           \
  {                                                                            \
    if(len == 0) return 0;                                                     \
    strbuf_ensure_capacity(sbuf, sbuf->end + len);                             \
    long nread;                                                                \
    if((nread = (long)__read(file,sbuf->b+sbuf->end,len)) <= 0) return 0;      \
    sbuf->end += (size_t)nread;                                                \
    return (size_t)nread;                                                      \
  }

_func_read(strbuf_gzread, gzFile, gzread2)
_func_read(strbuf_fread, FILE*, fread2)

// read FILE
// returns number of characters read
// or 0 if EOF
size_t strbuf_reset_readline(StrBuf *sbuf, FILE *file)
{
  strbuf_reset(sbuf);
  return strbuf_readline(sbuf, file);
}

// read gzFile
// returns number of characters read
// or 0 if EOF
size_t strbuf_reset_gzreadline(StrBuf *sbuf, gzFile file)
{
  strbuf_reset(sbuf);
  return strbuf_gzreadline(sbuf, file);
}

/**********/
/*  trim  */
/**********/

// Trim whitespace characters from the start and end of a string
void strbuf_trim(StrBuf *sbuf)
{
  if(sbuf->end == 0)
    return;

  // Trim end first
  while(sbuf->end > 0 && isspace(sbuf->b[sbuf->end-1]))
    sbuf->end--;

  sbuf->b[sbuf->end] = '\0';

  if(sbuf->end == 0)
    return;

  size_t start = 0;

  while(start < sbuf->end && isspace(sbuf->b[start]))
    start++;

  if(start != 0)
  {
    sbuf->end -= start;
    memmove(sbuf->b, sbuf->b+start, sbuf->end * sizeof(char));
    sbuf->b[sbuf->end] = '\0';
  }
}

// Trim the characters listed in `list` from the left of `sbuf`
// `list` is a null-terminated string of characters
void strbuf_ltrim(StrBuf *sbuf, const char *list)
{
  size_t start = 0;

  while(start < sbuf->end && (strchr(list, sbuf->b[start]) != NULL))
    start++;

  if(start != 0)
  {
    sbuf->end -= start;
    memmove(sbuf->b, sbuf->b+start, sbuf->end * sizeof(char));
    sbuf->b[sbuf->end] = '\0';
  }
}

// Trim the characters listed in `list` from the right of `sbuf`
// `list` is a null-terminated string of characters
void strbuf_rtrim(StrBuf *sbuf, const char *list)
{
  if(sbuf->end == 0)
    return;

  while(sbuf->end > 0 && strchr(list, sbuf->b[sbuf->end-1]) != NULL)
    sbuf->end--;

  sbuf->b[sbuf->end] = '\0';
}

/**************************/
/* Other String Functions */
/**************************/

// `n` is the maximum number of bytes to copy including the NULL byte
// copies at most n bytes from `src` to `dst`
// Always appends a NULL terminating byte, unless n is zero.
// Returns a pointer to dst
char* string_safe_ncpy(char *dst, const char *src, size_t n)
{
  if(n == 0) return dst;

  // From The Open Group:
  //   The memccpy() function copies bytes from memory area s2 into s1, stopping
  //   after the first occurrence of byte c is copied, or after n bytes are copied,
  //   whichever comes first. If copying takes place between objects that overlap,
  //   the behaviour is undefined.
  // Returns NULL if character c was not found in the copied memory
  if(memccpy(dst, src, '\0', n-1) == NULL)
    dst[n-1] = '\0';

  return dst;
}

// Replaces `sep` with \0 in str
// Returns number of occurances of `sep` character in `str`
// Stores `nptrs` pointers in `ptrs`
size_t string_split_str(char *str, char sep, char **ptrs, size_t nptrs)
{
  size_t n = 1;

  if(*str == '\0') return 0;
  if(nptrs > 0) ptrs[0] = str;

  while((str = strchr(str, sep)) != NULL) {
    *str = '\0';
    str++;
    if(n < nptrs) ptrs[n] = str;
    n++;
  }
  return n;
}

// Replace one char with another in a string. Return number of replacements made
size_t string_char_replace(char *str, char from, char to)
{
  size_t n = 0;
  for(; *str; str++) {
    if(*str == from) { n++; *str = to; }
  }
  return n;
}

// Reverse a string region
void string_reverse_region(char *str, size_t length)
{
  char *a = str, *b = str + length - 1;
  char tmp;
  while(a < b) {
    tmp = *a; *a = *b; *b = tmp;
    a++; b--;
  }
}

char string_is_all_whitespace(const char *s)
{
  int i;
  for(i = 0; s[i] != '\0' && isspace(s[i]); i++);
  return (s[i] == '\0');
}

char* string_next_nonwhitespace(char *s)
{
  while(*s != '\0' && isspace(*s)) s++;
  return (*s == '\0' ? NULL : s);
}

// Strip whitespace the the start and end of a string.  
// Strips whitepace from the end of the string with \0, and returns pointer to
// first non-whitespace character
char* string_trim(char *str)
{
  // Work backwards
  char *end = str+strlen(str);
  while(end > str && isspace(*(end-1))) end--;
  *end = '\0';

  // Work forwards: don't need start < len because will hit \0
  while(isspace(*str)) str++;

  return str;
}

// Removes \r and \n from the ends of a string and returns the new length
size_t string_chomp(char *str, size_t len)
{
  while(len > 0 && (str[len-1] == '\r' || str[len-1] == '\n')) len--;
  str[len] = '\0';
  return len;
}

// Returns count
size_t string_count_char(const char *str, char c)
{
  size_t count = 0;

  while((str = strchr(str, c)) != NULL)
  {
    str++;
    count++;
  }

  return count;
}

// Returns the number of strings resulting from the split
size_t string_split(const char *split, const char *txt, char ***result)
{
  size_t split_len = strlen(split);
  size_t txt_len = strlen(txt);

  // result is temporarily held here
  char **arr;

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
      arr = malloc(txt_len * sizeof(char*));
    
      size_t i;

      for(i = 0; i < txt_len; i++)
      {
        arr[i] = malloc(2 * sizeof(char));
        arr[i][0] = txt[i];
        arr[i][1] = '\0';
      }

      *result = arr;
      return txt_len;
    }
  }
  
  const char *find = txt;
  size_t count = 1; // must have at least one item

  for(; (find = strstr(find, split)) != NULL; count++, find += split_len) {}

  // Create return array
  arr = malloc(count * sizeof(char*));
  
  count = 0;
  const char *last_position = txt;

  size_t str_len;

  while((find = strstr(last_position, split)) != NULL)
  {
    str_len = (size_t)(find - last_position);

    arr[count] = malloc((str_len+1) * sizeof(char));
    strncpy(arr[count], last_position, str_len);
    arr[count][str_len] = '\0';
    
    count++;
    last_position = find + split_len;
  }

  // Copy last item
  str_len = (size_t)(txt + txt_len - last_position);
  arr[count] = malloc((str_len+1) * sizeof(char));

  if(count == 0) strcpy(arr[count], txt);
  else           strncpy(arr[count], last_position, str_len);

  arr[count][str_len] = '\0';
  count++;
  
  *result = arr;
  
  return count;
}
