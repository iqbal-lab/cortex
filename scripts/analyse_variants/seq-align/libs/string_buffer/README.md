C String Buffer
===============
Library code for handling strings and reading from files  
url: https://github.com/noporpoise/StringBuffer  
author: Isaac Turner <turner.isaac@gmail.com>  
license: Public Domain  
Jan 2015

[![Build Status](https://travis-ci.org/noporpoise/string_buffer.png?branch=master)](https://travis-ci.org/noporpoise/string_buffer)

About
=====

A string buffer library for C. Only has zlib as a dependency. Compiles with gcc
and clang.

Features:
- copying, inserting, appending, substring, chomp, trim
- reverse region, convert to upper/lower case
- sprintf into string buffer
- read a line at a time from a file
- buffered reading (10X faster for older versions of zlib)
- gzip file support

To build the test code:

    $ make
    $ ./strbuf_test

Calling
=======

To use in your code, include the following arguments in your gcc command:

    gcc ... -I$(STRING_BUF_PATH) -L$(STRING_BUF_PATH) ... -lstrbuf -lz

and include in your source code:

    include "string_buffer.h"

Example Code
============

    StrBuf* myBuff = strbuf_new()

    // Read from a file:

    gzFile fgz = gzopen("path/here.txt.gz")

    while(strbuf_gzreadline(myBuff, fgz))
    {
      // Do something with the line

      // e.g. chomp (remove newline)
      strbuf_chomp(myBuff)

      printf("%s\n", myBuff->b)

      // Reset StrBuf so you're not just concatenating all the lines in memory
      strbuf_reset(myBuff)
    }

    // Close up
    gzclose(fgz)

    strbuf_free(myBuff)


String buffers can still be used as input to standard str functions by accessing
the char* in the StrBuf struct. e.g.:

Get the position of the first 'a' in a StrBuf

    char *ptr = strchr(strbuf->b, 'a')
    int pos = (ptr == NULL ? -1 : ptr - strbuf->b)

Test if the StrBuf contains 'hello' from index 12

    if(strncasecmp(strbuf->b+12, "hello", 5) == 0)
      puts("world!\n")


API
===

Struct
------

    typedef struct
    {
      char *b;
      size_t end; // length of the string
      size_t size; // buffer size - includes '\0' (size is always >= len+1)
    } StrBuf;

Creators, destructors etc.
--------------------------

Constructors.  Note capacity increases as needed.

    StrBuf* strbuf_new()
    StrBuf* strbuf_init(const size_t capacity)
    StrBuf* strbuf_create(const char* str)

Place a string buffer into existing memory

    StrBuf* strbuf_alloc(StrBuf *sbuf, size_t capacity)

    // Example:
    StrBuf buf;
    strbuf_alloc(&buf, 100);

    // free malloc'd memory
    free(buf->seq.b);

Destructors

    void strbuf_free(StrBuf* sbuf)

Clone a string buffer (including content)

    StrBuf* strbuf_clone(const StrBuf* sbuf)

Clear the content of an existing StrBuf (sets size to 0)

    void strbuf_reset(StrBuf* sbuf)

Resizing
--------

Ensure capacity for len characters plus '\0' character - exits on FAILURE

    void strbuf_ensure_capacity(StrBuf *sbuf, const size_t len)

Resize the buffer to have capacity to hold a string of length new_len
(+ a null terminating character).  Can also be used to downsize the buffer's
memory usage.  Returns 1 on success, 0 on failure.

    char strbuf_resize(StrBuf *sbuf, const size_t new_size)

Useful functions
----------------

get/set chars

    char strbuf_get_char(const StrBuf *sbuf, const size_t index)
    void strbuf_set_char(StrBuf *sbuf, const size_t index, const char c)

Set string buffer to contain a given string
The string can be a string within the given string buffer

    void strbuf_set(StrBuf *sbuf, const char *str)

Get a copy of this StrBuf as a char array.
Returns NULL if not enough memory.
`strbuf_dup` is also provided as a shorthand.

    char* strbuf_as_str(const StrBuf* sbuf)

Add a character to the end of this StrBuf

    void strbuf_append_char(StrBuf* sbuf, const char txt)

Copy a StrBuf to the end of this StrBuf.
`strbuf_append` is also provided as a shorthand.

    void strbuf_append_buff(StrBuf* dst, StrBuf* src)

Copy a character array to the end of this StrBuf

    void strbuf_append_str(StrBuf* sbuf, const char* txt)

Copy N characters from a character array to the end of this StrBuf

    void strbuf_append_strn(StrBuf* sbuf, const char* txt, const size_t len)

Convert integers to strings and append the the end

    void strbuf_append_int(StrBuf *buf, int value)
    void strbuf_append_long(StrBuf *buf, long value)
    void strbuf_append_ulong(StrBuf *buf, unsigned long value)

Remove \r and \n characters from the end of this StrBuf.
Returns the number of characters removed

    size_t strbuf_chomp(StrBuf *sbuf)

Reverse a string

    void strbuf_reverse(StrBuf *sbuf)

Get a substring as a new null terminated char array
(remember to free the returned char* after you're done with it!)

    char* strbuf_substr(StrBuf *sbuf, const size_t start, const size_t len)

Change to upper or lower case

    void strbuf_to_uppercase(StrBuf *sbuf)
    void strbuf_to_lowercase(StrBuf *sbuf)

Copy a string to this StrBuf, overwriting any existing characters
Note: dst_pos + len can be longer the the current dst StrBuf

    void strbuf_copy(StrBuf* dst, size_t dst_pos,
                     const char* src, size_t len)

Insert: copy to a StrBuf, shifting any existing characters along

    void strbuf_insert(StrBuf* dst, size_t dst_pos,
                       const char* src, size_t len)


Overwrite `dst_pos..(dst_pos+dst_len-1)` with `src_len` chars from `src`.
If `dst_len != src_len`, content to the right of `dst_len` is shifted

    void strbuf_overwrite(StrBuf *dst, size_t dst_pos, size_t dst_len,
                          const char *src, size_t src_len)

    // Example:
    strbuf_set(sbuf, "aaabbccc");
    char *data = "xxx";
    strbuf_overwrite(sbuf, 3, 2, data, strlen(data));
    // sbuf is now "aaaxxxccc"
    strbuf_overwrite(sbuf, 3, 3, "_", 1);
    // sbuf is now "aaa_ccc"

Remove characters from the buffer

    void strbuf_delete(StrBuf *sbuf, size_t pos, size_t len)

    // Example:
    strbuf_set(sbuf, "aaaBBccc");
    strbuf_delete(sbuf, 3, 2);
    // sbuf is now "aaaccc"


Formatted strings (sprintf)
---------------------------

sprintf to a StrBuf (adds string terminator after sprint)

    int strbuf_sprintf(StrBuf *sbuf, const char* fmt, ...)
    int strbuf_sprintf_at(StrBuf *sbuf, const size_t pos, const char* fmt, ...)
    int strbuf_vsprintf(StrBuf *sbuf, const size_t pos,
                        const char* fmt, va_list argptr)

sprintf without terminating character.
Does not prematurely end the string if you sprintf within the string
(terminates string if sprintf to the end)

    int strbuf_sprintf_noterm(StrBuf *sbuf, const size_t pos,
                              const char* fmt, ...)

Reading files
-------------

Reading a FILE

    size_t strbuf_reset_readline(StrBuf *sbuf, FILE *file)
    size_t strbuf_readline(StrBuf *sbuf, FILE *gz_file)

Reading a gzFile

    size_t strbuf_reset_gzreadline(StrBuf *sbuf, gzFile gz_file)
    size_t strbuf_gzreadline(StrBuf *sbuf, gzFile gz_file)

Skip a line and return how many characters were skipped

    size_t strbuf_skipline(FILE *file)
    size_t strbuf_gzskipline(gzFile gz_file)

Read a line but no more than len bytes

    size_t strbuf_read(StrBuf *sbuf, FILE *file, size_t len)
    size_t strbuf_gzread(StrBuf *sbuf, gzFile file, size_t len)

Buffered reading

    size_t strbuf_gzreadline_buf(StrBuf *sbuf, gzFile gz_file, buffer_t *in);
    size_t strbuf_gzskipline_buf(gzFile file, buffer_t *in);

    size_t strbuf_readline_buf(StrBuf *sbuf, FILE *file, buffer_t *in);
    size_t strbuf_skipline_buf(FILE* file, buffer_t *in);

Example of buffered reading:

    gzFile gzf = gzopen("input.txt.gz", "r");
    buffer_t *in = buffer_new(1024); // pass buffer size in bytes
    StrBuf *line = strbuf_new();

    while(strbuf_gzreadline_buf(line, gzf, in) > 0)
    {
      strbuf_chomp(line);
      printf("read: %s\n", line->b);
    }

    strbuf_free(line);
    buffer_free(in);
    gzclose(gzf);

Read a line that has at least one character that is not \r or \n.
These functions do not call reset before reading.
Returns the number of characters read.

    size_t strbuf_readline_nonempty(StrBuf *line, FILE *fh)
    size_t strbuf_gzreadline_nonempty(StrBuf *line, gzFile gz)

Trim characters
---------------

Trim whitespace characters from the start and end of a string

    void strbuf_trim(StrBuf *sbuf)

Trim the characters listed in `list` from the left of `sbuf`.
`list` is a null-terminated string of characters

    void strbuf_ltrim(StrBuf *sbuf, const char* list)

Trim the characters listed in `list` from the right of `sbuf`.
`list` is a null-terminated string of characters

    void strbuf_rtrim(StrBuf *sbuf, const char* list)

Use with sscanf
---------------

To read strings into a string buffer using `sscanf`, first you must ensure the
buffer is big enough, and afterwards you must ensure the length is stored
correctly.

    StrBuf *sbuf = strbuf_new();
    char *input = "I'm sorry Dave I can't do that";
    
    strbuf_ensure_capacity(sbuf, strlen(input));
    sscanf(input, "I'm sorry %s I can't do that", sbuf->b);
    sbuf->end = strlen(sbuf->b);

    printf("Name: '%s'\n", sbuf->b);

    strbuf_free(sbuf);

Buffered input
--------------

`buffered_input.h` also provides generic buffered input functions

    buffer_t* buffer_new(size_t s)
    char buffer_init(buffer_t *b, size_t s)
    void buffer_free(buffer_t *b)
    void buffer_ensure_capacity(buffer_t *buf, size_t s)
    void buffer_append_str(buffer_t *buf, const char *str)
    void buffer_append_char(buffer_t *buf, char c)
    void buffer_terminate(buffer_t *b)
    void buffer_chomp(buffer_t *b)

Standardized gzFile and FILE versions of stream functions

    // Returns non-zero if an error has occurred
    int ferror2(FILE *fh);
    int gzerror2(gzFile gz);

    //
    // Unbuffered reading
    //

    // Return number of bytes read
    // Check ferror/gzerror on return for error
    size_t fread2(FILE *fh, buffer_t *buf, size_t len)
    size_t gzread2(gzFile gz, buffer_t *buf, size_t len)

    // Read up to n bytes or up to a newline
    // Return pointer to buffer read into or NULL if EOF
    // Check ferror/gzerror on return for error
    char* fgets2(FILE *fh, buffer_t *buf, size_t len)
    char* gzgets2(gzFile gz, buffer_t *buf, size_t len)

    // Read a line
    // Returns number of bytes read
    // Check ferror/gzerror on return for error
    size_t freadline(FILE* fh, char **buf, size_t *len, size_t *size)
    size_t gzreadline(gzFile gz, char **buf, size_t *len, size_t *size)

    // Skip a line
    // Returns number of bytes skipped
    // Check ferror/gzerror on return for error
    size_t fskipline(FILE* fh)
    size_t gzskipline(gzFile gz)

    //
    // Buffered reading
    //

    // Get a character
    // Check ferror/gzerror on return for error
    int fgetc_buf(FILE *fh, buffer_t *in)
    int gzgetc_buf(gzFile gz, buffer_t *in)

    // Read up to n bytes or up to a newline
    // Return pointer to buffer read into or NULL if EOF
    // Check ferror/gzerror on return for error
    char* fgets_buf(FILE* fh, buffer_t *in, char* str, size_t len)
    char* gzgets_buf(gzFile gz, buffer_t *in, char* str, size_t len)

    // Read a line
    // Returns number of bytes read
    // Check ferror/gzerror on return for error
    size_t freadline_buf(FILE* fh, buffer_t *in, char **buf, size_t *len, size_t *size)
    size_t gzreadline_buf(gzFile gz, buffer_t *in, char **buf, size_t *len, size_t *size)

    // Skip a line
    // Returns number of bytes skipped
    // Check ferror/gzerror on return for error
    int fskipline_buf(FILE* fh, buffer_t *in)
    int gzskipline_buf(gzFile gz, buffer_t *in)

    //
    // Unbuffered writing
    //

    // Write a single character to a stream
    // Returns a non-negative number, or –1 in case of error.
    int fputc2(FILE *fh, char c)
    int gzputc2(gzFile gz, char c)

    // Writes the given null-terminated string to a stream, excluding the
    // terminating null character
    // Returns a non-negative number, or –1 in case of error.
    int fputs2(FILE *fh, const char *str)
    int gzputs2(gzFile gz, const char *str)

    // Write to a stream
    // Check ferror/gzerror on return for error
    size_t fwrite2(FILE *fh, void *ptr, size_t len);
    int gzwrite2(gzFile gz, void *ptr, size_t len);


Other string functions
----------------------

These work on `char*` not `StrBuf`, but they're here because they're useful. 

Safely copy a string.
`n` is the maximum number of bytes to copy including the NULL byte.
Copies at most n bytes from `src` to `dst`.
Always appends a NULL terminating byte, unless n is zero.
Returns a pointer to `dst`.

    char* string_safe_ncpy(char *dst, const char *src, size_t n)

Other functions:

    size_t string_split_str(char *str, char sep, char **ptrs, size_t nptrs)
    size_t string_char_replace(char *str, char from, char to)
    void string_reverse_region(char *str, size_t length)
    char string_is_all_whitespace(const char* s)
    char* string_next_nonwhitespace(char* s)
    char* string_trim(char* str)
    size_t string_chomp(char* str)
    size_t string_count_char(const char* str, const int c)
    long string_split(const char* split, const char* txt, char*** result)


License
=======

This software is in the *Public Domain*. That means you can do whatever you like
with it. That includes being used in proprietary products without attribution or
restrictions. There are no warranties and there may be bugs. 

Formally we are using CC0 - a Creative Commons license to place this work in the
public domain. A copy of CC0 is in the LICENSE file. 

    "CC0 is a public domain dedication from Creative Commons. A work released
    under CC0 is dedicated to the public domain to the fullest extent permitted
    by law. If that is not possible for any reason, CC0 also provides a lax,
    permissive license as a fallback. Both public domain works and the lax
    license provided by CC0 are compatible with the GNU GPL."
      - http://www.gnu.org/licenses/license-list.html#CC0


Development
===========

short term goals: none - please suggest some!

I like to hear about how you're using it, what bugs you've found and what
features you'd like to see!  Contact me: Isaac Turner <turner.isaac@gmail>
