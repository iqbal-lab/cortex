/*
 strbuf_test.c
 project: string_buffer
 url: https://github.com/noporpoise/StringBuffer
 author: Isaac Turner <turner.isaac@gmail.com>
 license: Public Domain
 Jan 2015
*/

// request decent POSIX version
#define _XOPEN_SOURCE 700
#define _BSD_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>

#include "string_buffer.h"

#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define MIN(x,y) ((x) <= (y) ? (x) : (y))

const char tmp_file1[] = "tmp.strbuf.001.txt";
const char tmp_file2[] = "tmp.strbuf.002.txt";
const char tmp_gzfile1[] = "tmp.strbuf.001.txt.gz";
const char tmp_gzfile2[] = "tmp.strbuf.002.txt.gz";

/*********************/
/*  Testing toolkit  */
/*********************/

const char *suite_name;
char suite_pass;
size_t suites_run = 0, suites_failed = 0, suites_empty = 0;
size_t tests_passed = 0, tests_failed = 0;
size_t total_tests_passed = 0, total_tests_failed = 0;

#define QUOTE(str) #str

#define ASSERT(x)  do {                                                        \
    if(x) tests_passed++;                                                      \
    else {                                                                     \
      warn("failed assert [%s:%i] %s", __FILE__, __LINE__, QUOTE(x));          \
      suite_pass = 0; tests_failed++;                                          \
    }                                                                          \
  } while(0)

#define SUITE_START(x) {suite_pass = 1; suite_name = x; \
                        suites_run++; tests_passed = tests_failed = 0;}

#define SUITE_END()  do { \
    printf("Testing %s ", suite_name);                                         \
    size_t suite_i, tests_in_suite = tests_passed+tests_failed;                \
    for(suite_i = strlen(suite_name); suite_i<80-9-5; suite_i++) printf(".");  \
    printf("%s", suite_pass ? " pass" : " fail");                              \
    if(tests_in_suite == 0) { printf(" (empty)\n"); suites_empty++;}           \
    else if(tests_passed == 0 || tests_failed == 0) {                          \
      printf(" (%zu)\n", tests_in_suite);                                      \
    } else printf(" [%zu/%zu]\n", tests_passed, tests_in_suite);               \
    if(!suite_pass) suites_failed++;                                           \
    total_tests_failed += tests_failed;                                        \
    total_tests_passed += tests_passed;                                        \
  } while(0)

#define TEST_STATS()  do { \
    printf(" %zu / %zu suites failed\n", suites_failed, suites_run);           \
    printf(" %zu / %zu suites empty\n", suites_empty, suites_run);             \
    printf(" %zu / %zu tests failed\n", total_tests_failed,                    \
           total_tests_passed+total_tests_failed);                             \
  } while(0)

/* Test MACROs specifically for strbuf */

#define ASSERT_VALID(x) do {              \
    ASSERT((x)->end == strlen((x)->b));   \
    ASSERT((x)->end < (x)->size);         \
    ASSERT((x)->b[(x)->end] == '\0');     \
  } while(0)

/**********************/
/*  Output functions  */
/**********************/

void die(const char *fmt, ...)
__attribute__((format(printf, 1, 2)))
__attribute__((noreturn));

void warn(const char *fmt, ...)
__attribute__((format(printf, 1, 2)));

void die(const char *fmt, ...)
{
  fflush(stdout);

  // Print error
  fprintf(stderr, "Error: ");

  va_list argptr;
  va_start(argptr, fmt);
  vfprintf(stderr, fmt, argptr);
  va_end(argptr);

  // Check if we need to print a newline
  if(*(fmt + strlen(fmt) - 1) != '\n') fprintf(stderr, "\n");

  exit(EXIT_FAILURE);
}

void warn(const char *fmt, ...)
{
  fflush(stdout);

  // Print warning
  fprintf(stderr, "Warning: ");

  va_list argptr;
  va_start(argptr, fmt);
  vfprintf(stderr, fmt, argptr);
  va_end(argptr);

  // Check if we need to print a newline
  if(*(fmt + strlen(fmt) - 1) != '\n') fprintf(stderr, "\n");

  fflush(stderr);
}

// Random string of length len plus a null terminator
static void random_str(char *str, size_t len)
{
  size_t i, r = 0;
  for(i = 0; i < len; i++) {
    if(!(i & 7)) { r = (size_t)rand(); }
    str[i] = ' ' + ((r&0xff)%('~'-' '));
  }
  str[len] = '\0';
}

/* Tests! */

void test_roundup2pow()
{
  SUITE_START("roundup2pow");
  ASSERT(ROUNDUP2POW(0) == 0);
  ASSERT(ROUNDUP2POW(1) == 1);
  ASSERT(ROUNDUP2POW(2) == 2);
  ASSERT(ROUNDUP2POW(3) == 4);
  ASSERT(ROUNDUP2POW(4) == 4);
  ASSERT(ROUNDUP2POW(5) == 8);
  ASSERT(ROUNDUP2POW(6) == 8);
  ASSERT(ROUNDUP2POW(7) == 8);
  ASSERT(ROUNDUP2POW(8) == 8);
  ASSERT(ROUNDUP2POW(255) == 256);
  ASSERT(ROUNDUP2POW(256) == 256);
  SUITE_END();
}

/************************/
/* Buffered input tests */
/************************/

void test_buffers()
{
  SUITE_START("buffers");

  // Test buffer_init, cbuf_append_str, cbuf_append_char etc

  char *buf = malloc(4);
  size_t len = 0, size = 4;

  cbuf_append_char(&buf, &len, &size, 'a');
  cbuf_append_char(&buf, &len, &size, 'b');
  cbuf_append_char(&buf, &len, &size, 'c');
  cbuf_append_char(&buf, &len, &size, 'd');
  cbuf_append_char(&buf, &len, &size, 'e');

  ASSERT(len == 5);
  ASSERT(size >= 6); // strlen + '\0'
  ASSERT(strcmp(buf, "abcde") == 0);

  // Causes expansion -- tests ensure_capacity
  const char addstr[] = "fghijklmnopqrstuvwxyz";
  cbuf_append_str(&buf, &len, &size, addstr, strlen(addstr));

  ASSERT(len == 26);
  ASSERT(size >= 27); // strlen + '\0'
  ASSERT(strcmp(buf, "abcdefghijklmnopqrstuvwxyz") == 0);

  cbuf_append_char(&buf, &len, &size, '\r');
  cbuf_append_char(&buf, &len, &size, '\n');
  cbuf_chomp(buf, &len);

  ASSERT(len == 26);
  ASSERT(size >= 27); // strlen + '\0'
  ASSERT(strcmp(buf, "abcdefghijklmnopqrstuvwxyz") == 0);

  free(buf);

  SUITE_END();
}

// Compare buffered vs unbuffered + gzfile vs FILE
void test_buffered_reading()
{
  SUITE_START("buffered reading (getc/ungetc/gets/read/readline/skipline)");

  // generate file
  char *tmp = malloc(10000);
  char *tmpptr = tmp;
  strcpy(tmpptr,"hi\nThis is\nOur file\r\n");
  tmpptr += strlen("hi\nThis is\nOur file\r\n");
  int i;
  for(i = 0; i < 1000; i++) *(tmpptr++) = 'a';
  *(tmpptr++) = '\n';
  strcpy(tmpptr, "That's all folks!");

  // Save the string to four files
  FILE *file1 = fopen(tmp_file1, "w");
  fputs2(file1, tmp);
  fclose(file1);

  FILE *file2 = fopen(tmp_file2, "w");
  fputs2(file2, tmp);
  fclose(file2);

  gzFile gzfile1 = gzopen(tmp_gzfile1, "w");
  gzputs(gzfile1, tmp);
  gzclose(gzfile1);

  gzFile gzfile2 = gzopen(tmp_gzfile2, "w");
  gzputs(gzfile2, tmp);
  gzclose(gzfile2);

  free(tmp);

  // Open the files for reading
  file1 = fopen(tmp_file1, "r");
  file2 = fopen(tmp_file2, "r");
  gzfile1 = gzopen(tmp_gzfile1, "r");
  gzfile2 = gzopen(tmp_gzfile2, "r");

  StreamBuffer *fbuf = strm_buf_new(12);
  StreamBuffer *gzbuf = strm_buf_new(12);

  ASSERT(fbuf != NULL);
  ASSERT(gzbuf != NULL);

  StrBuf *st1 = strbuf_new(10);
  StrBuf *st2 = strbuf_new(10);
  StrBuf *st3 = strbuf_new(10);
  StrBuf *st4 = strbuf_new(10);

  // getc
  int c1 = fgetc(file1);
  int c2 = fgetc_buf(file2, fbuf);
  int c3 = gzgetc(gzfile1);
  int c4 = gzgetc_buf(gzfile2, gzbuf);

  ASSERT(c1 == 'h');
  ASSERT(c2 == 'h');
  ASSERT(c3 == 'h');
  ASSERT(c4 == 'h');

  // ungetc
  ASSERT(ungetc(c1, file1) == c1);
  ASSERT(ungetc_buf(c2, fbuf) == c2);
  ASSERT(gzungetc(c3, gzfile1) == c3);
  ASSERT(ungetc_buf(c4, gzbuf) == c4);

  // getc again
  c1 = fgetc(file1);
  c2 = fgetc_buf(file2, fbuf);
  c3 = gzgetc(gzfile1);
  c4 = gzgetc_buf(gzfile2, gzbuf);

  ASSERT(c1 == 'h');
  ASSERT(c2 == 'h');
  ASSERT(c3 == 'h');
  ASSERT(c4 == 'h');

  // readline
  ASSERT(freadline(file1, &(st1->b), &(st1->end), &(st1->size)) > 0);
  ASSERT(freadline_buf(file2, fbuf, &(st2->b), &(st2->end), &(st2->size)) > 0);
  ASSERT(gzreadline(gzfile1, &(st3->b), &(st3->end), &(st3->size)) > 0);
  ASSERT(gzreadline_buf(gzfile2, gzbuf, &(st4->b), &(st4->end), &(st4->size)) > 0);

  ASSERT(strcmp(st1->b, "i\n") == 0);
  ASSERT(strcmp(st2->b, "i\n") == 0);
  ASSERT(strcmp(st3->b, "i\n") == 0);
  ASSERT(strcmp(st4->b, "i\n") == 0);
  ASSERT(st1->end == 2);
  ASSERT(st2->end == 2);
  ASSERT(st3->end == 2);
  ASSERT(st4->end == 2);

  const char *lines[] = {"This is\n","Our file\r\n"};

  st1->end = st2->end = st3->end = st4->end = 0;

  // readline
  freadline(file1, &st1->b, &st1->end, &st1->size);
  freadline_buf(file2, fbuf, &st2->b, &st2->end, &st2->size);
  gzreadline(gzfile1, &st3->b, &st3->end, &st3->size);
  gzreadline_buf(gzfile2, gzbuf, &st4->b, &st4->end, &st4->size);

  ASSERT(strcmp(st1->b, lines[0]) == 0);
  ASSERT(strcmp(st2->b, lines[0]) == 0);
  ASSERT(strcmp(st3->b, lines[0]) == 0);
  ASSERT(strcmp(st4->b, lines[0]) == 0);
  ASSERT(st1->end == strlen(lines[0]));
  ASSERT(st2->end == strlen(lines[0]));
  ASSERT(st3->end == strlen(lines[0]));
  ASSERT(st4->end == strlen(lines[0]));

  // skipline
  fskipline(file1);
  fskipline_buf(file2, fbuf);
  gzskipline(gzfile1);
  gzskipline_buf(gzfile2, gzbuf);

  // gets
  ASSERT(fgets2(file1, st1->b, 10) != NULL);
  ASSERT(fgets_buf(file2, fbuf, st2->b, 10) != NULL);
  ASSERT(gzgets2(gzfile1, st3->b, 10) != NULL);
  ASSERT(gzgets_buf(gzfile2, gzbuf, st4->b, 10) != NULL);

  const char expected[] = "aaaaaaaaa";
  ASSERT(strcmp(st1->b, expected) == 0);
  ASSERT(strcmp(st2->b, expected) == 0);
  ASSERT(strcmp(st3->b, expected) == 0);
  ASSERT(strcmp(st4->b, expected) == 0);

  st1->end = strlen(st1->b);
  st2->end = strlen(st1->b);
  st3->end = strlen(st1->b);
  st4->end = strlen(st1->b);

  // readline
  freadline(file1, &st1->b, &st1->end, &st1->size);
  freadline_buf(file2, fbuf, &st2->b, &st2->end, &st2->size);
  gzreadline(gzfile1, &st3->b, &st3->end, &st3->size);
  gzreadline_buf(gzfile2, gzbuf, &st4->b, &st4->end, &st4->size);

  for(i = 0; i < 1000; i++)
  {
    ASSERT(st1->b[i] == 'a');
    ASSERT(st2->b[i] == 'a');
    ASSERT(st3->b[i] == 'a');
    ASSERT(st4->b[i] == 'a');
  }
  ASSERT(st1->b[i] == '\n');
  ASSERT(st2->b[i] == '\n');
  ASSERT(st3->b[i] == '\n');
  ASSERT(st4->b[i] == '\n');

  st1->end = st2->end = st3->end = st4->end = 0;

  // gets
  ASSERT(fgets2(file1, st1->b, (unsigned int)st1->size) != NULL);
  ASSERT(fgets_buf(file2, fbuf, st2->b, (unsigned int)st2->size) != NULL);
  ASSERT(gzgets2(gzfile1, st3->b, (unsigned int)st3->size) != NULL);
  ASSERT(gzgets_buf(gzfile2, gzbuf, st4->b, (unsigned int)st4->size) != NULL);

  ASSERT(strcmp(st1->b, "That's all folks!") == 0);
  ASSERT(strcmp(st2->b, "That's all folks!") == 0);
  ASSERT(strcmp(st3->b, "That's all folks!") == 0);
  ASSERT(strcmp(st4->b, "That's all folks!") == 0);

  // Check file/buffers empty
  ASSERT(fgetc(file1) == -1);
  ASSERT(fgetc_buf(file2, fbuf) == -1);
  ASSERT(gzgetc(gzfile1) == -1);
  ASSERT(gzgetc_buf(gzfile2, gzbuf) == -1);

  // close files
  fclose(file1);
  fclose(file2);
  gzclose(gzfile1);
  gzclose(gzfile2);

  // Open the files for reading again
  file1 = fopen(tmp_file1, "r");
  file2 = fopen(tmp_file2, "r");
  gzfile1 = gzopen(tmp_gzfile1, "r");
  gzfile2 = gzopen(tmp_gzfile2, "r");

  // Reset strings
  st1->end = st2->end = st3->end = st4->end = 0;

  // Rest buffers
  fbuf->begin = fbuf->end = gzbuf->begin = gzbuf->end = 0;

  // Read lines from file1 and compare to fread results on other files
  while(freadline(file1, &st1->b, &st1->end, &st1->size) > 0)
  {
    ASSERT(fread_buf(file2, st2->b, st1->end, fbuf) == st1->end);
    ASSERT(gzread(gzfile1, st3->b, (unsigned int)st1->end) == (int)st1->end);
    ASSERT(gzread_buf(gzfile2, st4->b, st1->end, gzbuf) == st1->end);

    // Null terminate since fread doesn't do that
    st2->b[st1->end] = st3->b[st1->end] = st4->b[st1->end] = '\0';

    ASSERT(strcmp(st1->b,st2->b) == 0);
    ASSERT(strcmp(st1->b,st3->b) == 0);
    ASSERT(strcmp(st1->b,st4->b) == 0);

    st1->end = 0;
  }

  // close files
  fclose(file1);
  fclose(file2);
  gzclose(gzfile1);
  gzclose(gzfile2);

  // free
  strbuf_free(st1);
  strbuf_free(st2);
  strbuf_free(st3);
  strbuf_free(st4);
  strm_buf_free(gzbuf);
  strm_buf_free(fbuf);

  SUITE_END();
}


/***********************/
/* String buffer tests */
/***********************/

void _test_clone(const char *str)
{
  StrBuf *a = strbuf_create(str);
  StrBuf *b = strbuf_clone(a);

  ASSERT(strcmp(a->b,b->b) == 0);
  ASSERT(a->size == b->size);
  ASSERT(a->end == b->end);
  ASSERT_VALID(a);
  ASSERT_VALID(b);

  strbuf_free(a);
  strbuf_free(b);
}

void test_clone()
{
  SUITE_START("clone");

  _test_clone("");
  _test_clone("ASDFASDFASDFASDF");
  _test_clone("                                                               "
              "                                                               "
              "                                                               "
              "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n"
              "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa");
  _test_clone("0");
  _test_clone("\n");
  _test_clone(" ");

  SUITE_END();
}


void _test_reset(const char *str)
{
  StrBuf *a = strbuf_create(str);
  size_t capacity = a->size;
  strbuf_reset(a);

  ASSERT(a->b[0] == '\0');
  ASSERT(a->end == 0);
  ASSERT(a->size == capacity);

  strbuf_free(a);
}

void test_reset()
{
  SUITE_START("reset");

  _test_reset("a");
  _test_reset("ab");
  _test_reset("abc");
  _test_clone("                                                              ");
  _test_reset("");

  SUITE_END();
}

// Test resize and ensure_capacity
// Test that resize can shrink strbuf safely
void _test_resize(const char *str, size_t new_len)
{
  StrBuf *sbuf = strbuf_create(str);
  ASSERT(strcmp(str, sbuf->b) == 0);
  ASSERT_VALID(sbuf);

  strbuf_resize(sbuf, new_len);
  ASSERT_VALID(sbuf);
  ASSERT(sbuf->end == MIN(new_len, strlen(str)));
  ASSERT(strncmp(str, sbuf->b, MIN(new_len, strlen(str))) == 0);

  strbuf_free(sbuf);
}
void test_resize()
{
  SUITE_START("resize");

  _test_resize("", 0);
  _test_resize("", 10000);
  _test_resize("abc", 10000);
  _test_resize("abc", 0);
  _test_resize("abc", 1);
  _test_resize("abc", 3);
  _test_resize("abcdefghijklmnopqrstuvwxyz", 0);
  _test_resize("abcdefghijklmnopqrstuvwxyz", 1);
  _test_resize("abcdefghijklmnopqrstuvwxyz", 10);
  _test_resize("abcdefghijklmnopqrstuvwxyz", 255);
  _test_resize("abcdefghijklmnopqrstuvwxyz", 10000);

  SUITE_END();
}

void test_get_set_char()
{
  SUITE_START("get_char / set_char");

  StrBuf *sbuf = strbuf_create("abcd");

  strbuf_char(sbuf, 0) = 'z';
  strbuf_char(sbuf, 1) = 'y';
  ASSERT(strcmp(sbuf->b, "zycd") == 0);
  strbuf_char(sbuf, 2) = 'x';
  strbuf_char(sbuf, 3) = 'w';
  ASSERT(strcmp(sbuf->b, "zyxw") == 0);
  strbuf_append_char(sbuf, 'v');
  strbuf_append_char(sbuf, 'u');
  ASSERT(strcmp(sbuf->b, "zyxwvu") == 0);

  ASSERT(strbuf_char(sbuf, 0) == 'z');
  ASSERT(strbuf_char(sbuf, 1) == 'y');
  ASSERT(strbuf_char(sbuf, 2) == 'x');
  ASSERT(strbuf_char(sbuf, 3) == 'w');
  ASSERT(strbuf_char(sbuf, 4) == 'v');
  ASSERT(strbuf_char(sbuf, 5) == 'u');

  strbuf_free(sbuf);

  SUITE_END();
}


void _test_set(StrBuf *sbuf, const char *str)
{
  strbuf_set(sbuf, str);
  ASSERT(strcmp(sbuf->b, str) == 0);
  ASSERT_VALID(sbuf);
}
void test_set()
{
  SUITE_START("set");
  StrBuf *sbuf = strbuf_create("abcd");

  _test_set(sbuf, "abcd");
  _test_set(sbuf, "");
  _test_set(sbuf, "a");
  _test_set(sbuf, "abcdefghijklmnopqrstuvwxyzABCDEFGHIJLKMNOPQRSTUVWXY");
  _test_set(sbuf, "ab");
  _test_set(sbuf, "abc");
  _test_set(sbuf, "");

  strbuf_free(sbuf);
  SUITE_END();
}

void _test_as_str(StrBuf *sbuf, const char *str)
{
  strbuf_set(sbuf, str);
  char *tmp = strbuf_dup_str(sbuf);
  ASSERT(strcmp(str, sbuf->b) == 0);
  ASSERT(strcmp(tmp, str) == 0);
  ASSERT_VALID(sbuf);
  free(tmp);
}
void test_as_str()
{
  SUITE_START("as_str");
  StrBuf *sbuf = strbuf_new(10);

  _test_as_str(sbuf, "");
  _test_as_str(sbuf, "a");
  _test_as_str(sbuf, "ab");
  _test_as_str(sbuf, "abc");
  _test_as_str(sbuf, "abcdefghijklmnopqrstuvwxyzABCDEFGHIJLKMNOPQRSTUVWXY");
  _test_as_str(sbuf, "abc");

  strbuf_free(sbuf);
  SUITE_END();
}

void _test_append(StrBuf* sbuf, char c, const char *str, const char *str2,
                  size_t n, const StrBuf *append)
{
  size_t len = sbuf->end;
  size_t extend = append->end;
  size_t end = len + extend;

  ASSERT(strlen(str2) >= n);

  strbuf_append_buff(sbuf, append);
  strbuf_append_char(sbuf, c);
  strbuf_append_str(sbuf, str);
  strbuf_append_strn(sbuf, str2, n);

  size_t str1len = strlen(str);

  ASSERT(strncmp(sbuf->b+len, append->b, extend) == 0);
  ASSERT(sbuf->b[end] == c);
  ASSERT(strncmp(sbuf->b+end+1, str, str1len) == 0);
  ASSERT(strncmp(sbuf->b+end+1+str1len, str2, n) == 0);

  ASSERT(sbuf->end == len + 1 + str1len + n + extend);
  ASSERT_VALID(sbuf);
}
void test_append()
{
  SUITE_START("append_char / append_buff / append_str / append_strn");
  StrBuf *sbuf = strbuf_new(10);

  _test_append(sbuf, 'a', "", "", 0, sbuf);
  _test_append(sbuf, 'b', "a", "xxy", 1, sbuf);
  _test_append(sbuf, 'c', "a", "xxy", 3, sbuf);
  _test_append(sbuf, 'd', "abcdefghijklmno", "abcdefghijklmno", 0, sbuf);
  _test_append(sbuf, 'd', "abcdefghijklmno", "abcdefghijklmno", 15, sbuf);
  _test_append(sbuf, 'd', "", "", 0, sbuf);

  StrBuf *empty = strbuf_new(10);
  _test_append(sbuf, 'd', "", "", 0, empty);
  _test_append(empty, 'd', "", "", 0, sbuf);

  strbuf_free(empty);
  strbuf_free(sbuf);
  SUITE_END();
}

void _test_append_long(StrBuf *sbuf, long x)
{
  char truth[50];
  sprintf(truth, "%li", x);
  strbuf_reset(sbuf);
  // Initialise the string with some noise
  strbuf_append_charn(sbuf, 'x', rand()&0xffff); // 65536
  size_t start = sbuf->end;
  strbuf_append_long(sbuf, x);
  ASSERT(strcmp(sbuf->b+start, truth) == 0);
}

void test_append_int()
{
  SUITE_START("using append_int()");
  long repeat, x;

  for(repeat = 0; repeat < 20; repeat++)
  {
    StrBuf *sbuf = strbuf_create("");

    for(x = -300; x <= 2000; x++) {
      _test_append_long(sbuf, x);
    }

    _test_append_long(sbuf, 123456789);
    _test_append_long(sbuf, 9999999);

    _test_append_long(sbuf, LONG_MIN);
    _test_append_long(sbuf, LONG_MAX);

    strbuf_free(sbuf);
  }

  SUITE_END();
}

void _test_chomp(const char *str)
{
  size_t len = strlen(str);
  size_t trim = len;
  while(trim > 0 && (str[trim-1] == '\r' || str[trim-1] == '\n')) trim--;

  StrBuf *sbuf = strbuf_create(str);
  ASSERT_VALID(sbuf);
  size_t buf_len = sbuf->end;

  strbuf_chomp(sbuf);
  ASSERT_VALID(sbuf);
  size_t buf_trim = sbuf->end;

  ASSERT(buf_len == len);
  ASSERT(buf_trim == trim);

  strbuf_free(sbuf);
}
void test_chomp()
{
  SUITE_START("chomp");

  _test_chomp("\n");
  _test_chomp("");
  _test_chomp("\r\n");
  _test_chomp("asdfa\nasdf");
  _test_chomp("asdfa\n\r");
  _test_chomp("asdfa\r\n");
  _test_chomp("asdfa\n");
  _test_chomp("asdfa\n ");

  SUITE_END();
}

void _test_reverse(const char *str)
{
  size_t len = strlen(str);
  StrBuf *sbuf = strbuf_create(str);
  strbuf_reverse(sbuf);
  ASSERT(sbuf->end == len);
  ASSERT_VALID(sbuf);

  size_t i;
  for(i = 0; i < len; i++) ASSERT(str[i] == sbuf->b[len-i-1]);

  strbuf_free(sbuf);
}
void test_reverse()
{
  SUITE_START("reverse");

  _test_reverse("");
  _test_reverse("ASDFASDF");
  _test_reverse("   ");
  _test_reverse("\n\n\n");
  _test_reverse("abcdefghijklmnopqrstuvwxyz");

  SUITE_END();
}


void _test_substr(const char *str, size_t start, size_t len)
{
  StrBuf *sbuf = strbuf_create(str);
  ASSERT(strcmp(sbuf->b, str) == 0);
  ASSERT_VALID(sbuf);

  char *tmp = strbuf_substr(sbuf, start, len);
  ASSERT(strncmp(tmp, sbuf->b+start, len) == 0);
  ASSERT(len == strlen(tmp));

  free(tmp);
  strbuf_free(sbuf);
}
void test_substr()
{
  SUITE_START("substr");

  _test_substr("", 0, 0);
  _test_substr("a", 0, 1);
  _test_substr("a", 0, 0);
  _test_substr("a", 1, 0);
  _test_substr("abcdef", 3, 0);
  _test_substr("abcdef", 3, 1);
  _test_substr("abcdef", 3, 3);
  _test_substr("abcdef", 0, 6);
  _test_substr("abcdef", 2, 2);
  _test_substr("abcdefghijklmnopqrstuvwxyz", 0, 1);
  _test_substr("abcdefghijklmnopqrstuvwxyz", 25, 1);
  _test_substr("abcdefghijklmnopqrstuvwxyz", 5, 5);

  SUITE_END();
}

void _test_change_case(const char *str)
{
  StrBuf *sbuf = strbuf_create(str);

  strbuf_to_uppercase(sbuf);
  ASSERT_VALID(sbuf);
  char *upper = strbuf_dup_str(sbuf);

  strbuf_to_lowercase(sbuf);
  ASSERT_VALID(sbuf);
  char *lower = strbuf_dup_str(sbuf);

  // Length checks
  size_t len = strlen(str);
  ASSERT(strlen(upper) == len);
  ASSERT(strlen(lower) == len);
  ASSERT(sbuf->end == len);

  size_t i;
  for(i = 0; i < len; i++) {
    ASSERT(upper[i] == toupper(str[i]));
    ASSERT(lower[i] == tolower(str[i]));
  }

  free(upper);
  free(lower);
  strbuf_free(sbuf);
}
void test_change_case()
{
  SUITE_START("uppercase / lowercase");

  _test_change_case("");
  _test_change_case("asdfasdf");
  _test_change_case("asdf");
  _test_change_case("ASDFASDF:. asdfasdf \nasdfasdf'aougyqvo23=-=12#");

  SUITE_END();
}


void _test_copy(StrBuf *sbuf, size_t pos, const char *from, size_t len)
{
  char *frmcpy = strdup(from);
  size_t orig_len = sbuf->end;

  char *orig = strbuf_dup_str(sbuf);
  ASSERT_VALID(sbuf);
  ASSERT(strcmp(sbuf->b, orig) == 0);

  strbuf_copy(sbuf, pos, from, len);

  ASSERT(sbuf->end == MAX(orig_len, pos+len));
  ASSERT_VALID(sbuf);

  ASSERT(strncmp(sbuf->b, orig, pos) == 0);
  ASSERT(strncmp(sbuf->b+pos, frmcpy, len) == 0);
  ASSERT(strncmp(sbuf->b+pos+len, orig+pos+len, sbuf->end-pos-len) == 0);

  free(frmcpy);
  free(orig);
}
void test_copy()
{
  SUITE_START("copy");
  StrBuf *sbuf = strbuf_create("");

  _test_copy(sbuf, 0, "", 0);
  _test_copy(sbuf, 0, "asdf", 0);
  _test_copy(sbuf, 0, "asdf", 1);

  strbuf_set(sbuf, "");
  _test_copy(sbuf, 0, "asdf", 4);

  size_t i, j;
  for(i = 0; i <= 4; i++)
  {
    strbuf_set(sbuf, "asdf");
    _test_copy(sbuf, i, "asdf", 2);
    strbuf_set(sbuf, "asdf");
    _test_copy(sbuf, i, "asdf", 4);
  }

  strbuf_set(sbuf, "asdfasdfasdf");
  _test_copy(sbuf, 8, "df", 2);
  strbuf_set(sbuf, "asdfasdfasdf");
  _test_copy(sbuf, 8, "", 0);

  strbuf_set(sbuf, "asdfasdfasdf");
  _test_copy(sbuf, 8, sbuf->b, sbuf->end);

  for(i = 0; i <= 4; i++)
  {
    for(j = 0; j <= 4; j++)
    {
      strbuf_set(sbuf, "asdf");
      _test_copy(sbuf, i, sbuf->b, j);
    }
  }

  strbuf_free(sbuf);
  SUITE_END();
}

void _test_insert(StrBuf *sbuf, size_t pos, size_t len, const char *from)
{
  char *frmcpy = (from == NULL ? calloc(1,1) : strdup(from));
  size_t orig_len = sbuf->end;

  char *orig = strbuf_dup_str(sbuf);
  ASSERT_VALID(sbuf);
  ASSERT(strcmp(sbuf->b, orig) == 0);

  strbuf_insert(sbuf, pos, from, len);
  ASSERT_VALID(sbuf);

  ASSERT(sbuf->end == orig_len + len);
  ASSERT(sbuf->end < sbuf->size);

  ASSERT(strncmp(sbuf->b, orig, pos) == 0);
  ASSERT(strncmp(sbuf->b+pos, frmcpy, len) == 0);
  ASSERT(strncmp(sbuf->b+pos+len, orig+pos, orig_len-pos) == 0);

  free(frmcpy);
  free(orig);
}
void test_insert()
{
  SUITE_START("insert");
  StrBuf *sbuf = strbuf_create("");

  _test_insert(sbuf, 0, 0, NULL);
  _test_insert(sbuf, 0, 0, "");
  _test_insert(sbuf, 0, 0, "asdf");
  _test_insert(sbuf, 0, 1, "asdf");

  strbuf_set(sbuf, "");
  _test_insert(sbuf, 0, 4, "asdf");

  size_t i, j;
  for(i = 0; i <= 4; i++)
  {
    strbuf_set(sbuf, "asdf");
    _test_insert(sbuf, i, 2, "asdf");
    strbuf_set(sbuf, "asdf");
    _test_insert(sbuf, i, 4, "asdf");
  }

  strbuf_set(sbuf, "asdfasdfasdf");
  _test_insert(sbuf, 8, 2, "df");
  strbuf_set(sbuf, "asdfasdfasdf");
  _test_insert(sbuf, 8, 0, "");

  strbuf_set(sbuf, "asdfasdfasdf");
  _test_insert(sbuf, 8, sbuf->end, sbuf->b);

  for(i = 0; i <= 4; i++)
  {
    for(j = 0; j <= 4; j++)
    {
      strbuf_set(sbuf, "asdf");
      _test_insert(sbuf, i, j, sbuf->b);
    }
  }

  strbuf_set(sbuf, "abcdefghij");
  strbuf_insert(sbuf, 3, sbuf->b+1, 5);
  ASSERT(strcmp(sbuf->b, "abcbcdefdefghij") == 0);
  ASSERT_VALID(sbuf);

  char *long_str = malloc(501);
  random_str(long_str, 500);

  const char *short_str = "GGTTCTTCTTGGCTTCTTCTTTTCATTGCC";
  strbuf_set(sbuf, long_str);
  strbuf_insert(sbuf, 0, short_str, strlen(short_str));
  ASSERT(strncmp(sbuf->b, short_str, strlen(short_str)) == 0);
  ASSERT(strcmp(sbuf->b+strlen(short_str), long_str) == 0);
  ASSERT_VALID(sbuf);
  free(long_str);

  strbuf_free(sbuf);
  SUITE_END();
}

void test_overwrite()
{
  SUITE_START("overwrite");
  StrBuf *sbuf = strbuf_new(10);

  strbuf_set(sbuf, "aaabbccc");

  strbuf_overwrite(sbuf, 3, 2, "BBB", 3);
  ASSERT(strcmp(sbuf->b, "aaaBBBccc") == 0);
  ASSERT_VALID(sbuf);

  strbuf_overwrite(sbuf, 3, 3, "_x", 1);
  ASSERT(strcmp(sbuf->b, "aaa_ccc") == 0);
  ASSERT_VALID(sbuf);

  strbuf_set(sbuf, "abcdefghijklmnopqrstuvwxyz");
  // replace de with abcdef
  strbuf_overwrite(sbuf, 3, 2, sbuf->b, 6);
  ASSERT(strcmp(sbuf->b, "abcabcdeffghijklmnopqrstuvwxyz") == 0);
  ASSERT_VALID(sbuf);

  // replace abcdef with de
  strbuf_overwrite(sbuf, 3, 6, sbuf->b+6, 2);
  ASSERT(strcmp(sbuf->b, "abcdefghijklmnopqrstuvwxyz") == 0);
  ASSERT_VALID(sbuf);

  // do nothing
  strbuf_overwrite(sbuf, 3, 0, sbuf->b+6, 0);
  ASSERT(strcmp(sbuf->b, "abcdefghijklmnopqrstuvwxyz") == 0);
  ASSERT_VALID(sbuf);

  // delete b
  strbuf_overwrite(sbuf, 1, 1, sbuf->b, 0);
  ASSERT(strcmp(sbuf->b, "acdefghijklmnopqrstuvwxyz") == 0);
  ASSERT_VALID(sbuf);

  // swap ghij with hi
  strbuf_overwrite(sbuf, 5, 4, sbuf->b+6, 2);
  ASSERT(strcmp(sbuf->b, "acdefhiklmnopqrstuvwxyz") == 0);
  ASSERT_VALID(sbuf);

  // replace o with z
  strbuf_overwrite(sbuf, 11, 1, sbuf->b+22, 1);
  ASSERT(strcmp(sbuf->b, "acdefhiklmnzpqrstuvwxyz") == 0);
  ASSERT_VALID(sbuf);

  // replace pq with stuv
  strbuf_overwrite(sbuf, 12, 2, sbuf->b+15, 4);
  ASSERT(strcmp(sbuf->b, "acdefhiklmnzstuvrstuvwxyz") == 0);
  ASSERT_VALID(sbuf);

  // replace stuv with e
  strbuf_overwrite(sbuf, 12, 4, sbuf->b+3, 1);
  ASSERT(strcmp(sbuf->b, "acdefhiklmnzerstuvwxyz") == 0);
  ASSERT_VALID(sbuf);

  // replace lmn with "A"
  strbuf_overwrite(sbuf, 8, 3, "AB", 1);
  ASSERT(strcmp(sbuf->b, "acdefhikAzerstuvwxyz") == 0);
  ASSERT_VALID(sbuf);

  // replace A with "XYZ"
  strbuf_overwrite(sbuf, 8, 1, "XYZ", 3);
  ASSERT(strcmp(sbuf->b, "acdefhikXYZzerstuvwxyz") == 0);
  ASSERT_VALID(sbuf);

  // replace XYZ with "Zz"
  strbuf_overwrite(sbuf, 8, 3, sbuf->b+10, 2);
  ASSERT(strcmp(sbuf->b, "acdefhikZzzerstuvwxyz") == 0);
  ASSERT_VALID(sbuf);

  // replace zer with "zz"
  strbuf_overwrite(sbuf, 10, 3, sbuf->b+9, 2);
  ASSERT(strcmp(sbuf->b, "acdefhikZzzzstuvwxyz") == 0);
  ASSERT_VALID(sbuf);

  strbuf_free(sbuf);
  SUITE_END();
}

void test_delete()
{
  SUITE_START("delete");
  StrBuf *sbuf = strbuf_new(10);

  strbuf_set(sbuf, "aaaBBccc");

  strbuf_delete(sbuf, 3, 2);
  ASSERT(strcmp(sbuf->b, "aaaccc") == 0);
  ASSERT_VALID(sbuf);

  strbuf_delete(sbuf, 3, 0);
  ASSERT(strcmp(sbuf->b, "aaaccc") == 0);
  ASSERT_VALID(sbuf);

  strbuf_delete(sbuf, 0, 0);
  ASSERT(strcmp(sbuf->b, "aaaccc") == 0);
  ASSERT_VALID(sbuf);

  strbuf_delete(sbuf, 0, 1);
  ASSERT(strcmp(sbuf->b, "aaccc") == 0);
  ASSERT_VALID(sbuf);

  strbuf_delete(sbuf, 4, 1);
  ASSERT(strcmp(sbuf->b, "aacc") == 0);
  ASSERT_VALID(sbuf);

  strbuf_delete(sbuf, 4, 0);
  ASSERT(strcmp(sbuf->b, "aacc") == 0);
  ASSERT_VALID(sbuf);

  strbuf_free(sbuf);
  SUITE_END();
}

void test_sprintf()
{
  SUITE_START("sprintf");
  StrBuf *sbuf = strbuf_new(10);

  // although valid, GCC complains about formatted strings of length 0
  #ifdef __clang__
  strbuf_sprintf(sbuf, "");
  ASSERT(strcmp(sbuf->b, "") == 0);
  ASSERT_VALID(sbuf);
  #endif

  strbuf_sprintf(sbuf, "hi. ");
  ASSERT(strcmp(sbuf->b, "hi. ") == 0);
  ASSERT_VALID(sbuf);

  // Note: strbuf_sprintf appends -> so we still have 'hi. '
  strbuf_sprintf(sbuf, "A dozen is another way of saying %i, except for bakers "
                       "where it means %lu for some reason.  No other "
                       "profession is known to have its own dozen", 12, 13UL);

  const char ans[] = "hi. A dozen is another way of saying 12, except for bakers "
                     "where it means 13 for some reason.  No other "
                     "profession is known to have its own dozen";

  ASSERT(strcmp(sbuf->b, ans) == 0);
  ASSERT_VALID(sbuf);

  strbuf_reset(sbuf);
  strbuf_sprintf(sbuf, "woot %s %i %c", "what excitement", 12, '?');
  ASSERT(strcmp(sbuf->b, "woot what excitement 12 ?") == 0);
  ASSERT_VALID(sbuf);

  strbuf_reset(sbuf);
  strbuf_sprintf(sbuf, "bye");
  ASSERT(strcmp(sbuf->b, "bye") == 0);
  ASSERT_VALID(sbuf);

  strbuf_free(sbuf);
  SUITE_END();
}

void test_sprintf_at()
{
  SUITE_START("sprintf_at");
  StrBuf *sbuf = strbuf_new(10);

  // although valid, GCC complains about formatted strings of length 0
  #ifdef __clang__
  strbuf_sprintf_at(sbuf, 0, "");
  ASSERT(strcmp(sbuf->b, "") == 0);
  ASSERT_VALID(sbuf);
  #endif

  strbuf_sprintf_at(sbuf, 0, "hi. ");
  ASSERT(strcmp(sbuf->b, "hi. ") == 0);
  ASSERT_VALID(sbuf);

  strbuf_sprintf_at(sbuf, 2, " bye. ");
  ASSERT(strcmp(sbuf->b, "hi bye. ") == 0);
  ASSERT_VALID(sbuf);

  strbuf_sprintf_at(sbuf, 0, "woot %s %i %c", "what excitement", 12, '?');
  ASSERT(strcmp(sbuf->b, "woot what excitement 12 ?") == 0);
  ASSERT_VALID(sbuf);

  strbuf_sprintf_at(sbuf, 5, "moo %i", 6);
  ASSERT(strcmp(sbuf->b, "woot moo 6") == 0);
  ASSERT_VALID(sbuf);

  strbuf_free(sbuf);
  SUITE_END();
}

void test_sprintf_noterm()
{
  SUITE_START("sprintf_noterm");
  StrBuf *sbuf = strbuf_new(10);

  // although valid, GCC complains about formatted strings of length 0
  #ifdef __clang__
  strbuf_sprintf_noterm(sbuf, 0, "");
  ASSERT(strcmp(sbuf->b, "") == 0);
  ASSERT_VALID(sbuf);
  #endif

  strbuf_sprintf_noterm(sbuf, 0, "hi. ");
  ASSERT(strcmp(sbuf->b, "hi. ") == 0);
  ASSERT_VALID(sbuf);

  strbuf_sprintf_noterm(sbuf, 2, " bye. ");
  ASSERT(strcmp(sbuf->b, "hi bye. ") == 0);
  ASSERT_VALID(sbuf);

  strbuf_sprintf_noterm(sbuf, 0, "woot %s %i %c", "what excitement", 12, '?');
  ASSERT(strcmp(sbuf->b, "woot what excitement 12 ?") == 0);
  ASSERT_VALID(sbuf);

  strbuf_sprintf_noterm(sbuf, 5, "moo %i", 6);
  ASSERT(strcmp(sbuf->b, "woot moo 6excitement 12 ?") == 0);
  ASSERT_VALID(sbuf);

  // although valid, GCC complains about formatted strings of length 0
  #ifdef __clang__
  StrBuf *sbuf2 = strbuf_clone(sbuf);
  strbuf_sprintf_noterm(sbuf, 5, "");
  ASSERT(strcmp(sbuf->b, sbuf2->b) == 0);
  ASSERT_VALID(sbuf2);
  strbuf_free(sbuf2);
  #endif

  strbuf_free(sbuf);
  SUITE_END();
}

#define ftest(fname,type_t,__open,__close,__puts,__readline,__skipline) \
  void fname(const char *path) \
  { \
    /* Generate file */ \
    int i; \
    type_t out = __open(path, "w"); \
    if(out == NULL) die("Couldn't open: %s", path); \
    __puts(out, "hi\nthis is\nour file\n"); \
    /* 4th line: print 1000x'a' */ \
    for(i = 0; i < 1000; i++) __puts(out, "a"); \
    __close(out); \
    type_t file = __open(path, "r"); \
    if(file == NULL) die("Couldn't open: %s", path); \
    StrBuf *line = strbuf_new(10); \
    __readline(line, file); \
    ASSERT(strcmp(line->b, "hi\n") == 0); \
    ASSERT(line->end == strlen(line->b)); \
    ASSERT(line->end < line->size); \
    strbuf_chomp(line); \
    ASSERT(strcmp(line->b, "hi") == 0); \
    ASSERT(line->end == strlen(line->b)); \
    ASSERT(line->end < line->size); \
    __skipline(file); \
    strbuf_reset(line); \
    __readline(line, file); \
    ASSERT(strcmp(line->b, "our file\n") == 0); \
    ASSERT(line->end == strlen(line->b)); \
    ASSERT(line->end < line->size); \
    strbuf_chomp(line); \
    ASSERT(strcmp(line->b, "our file") == 0); \
    ASSERT(line->end == strlen(line->b)); \
    ASSERT(line->end < line->size); \
    strbuf_reset(line); \
    __readline(line, file); \
    ASSERT(line->end == 1000); \
    ASSERT(line->end < line->size); \
    for(i = 0; i < 1000; i++) { ASSERT(line->b[i] == 'a'); } \
    ASSERT(line->b[1000] == '\0'); \
    strbuf_free(line); \
    __close(file); \
  }

ftest(test_gzfile,gzFile,gzopen,gzclose,gzputs2,strbuf_gzreadline,strbuf_gzskipline)
ftest(test_file,FILE*,fopen,fclose,fputs2,strbuf_readline,strbuf_skipline)

void test_read_gzfile()
{
  SUITE_START("gzreadline / gzskipline");
  test_gzfile(tmp_gzfile1);
  SUITE_END();
}

void test_read_file()
{
  SUITE_START("readline / skipline");
  test_file(tmp_file1);
  SUITE_END();
}

void test_read_nonempty()
{
  SUITE_START("read nonempty");

  FILE *fh = fopen(tmp_file1, "w");
  if(fh == NULL) die("Cannot write tmp output file: %s", tmp_file1);
  fprintf(fh, "hi\n\r\n\r\n"
              "bye\nx\ny\n"
              "\n\n\nz\n\n\n");
  fclose(fh);

  fh = fopen(tmp_file1, "r");
  if(fh == NULL) die("Cannot read tmp output file: %s", tmp_file1);

  StrBuf *sbuf = strbuf_new(10);
  ASSERT(strbuf_readline_nonempty(sbuf, fh));
  strbuf_chomp(sbuf);
  ASSERT(strcmp(sbuf->b,"hi")==0);

  strbuf_reset(sbuf);
  ASSERT(strbuf_readline_nonempty(sbuf, fh));
  strbuf_chomp(sbuf);
  ASSERT(strcmp(sbuf->b,"bye")==0);

  strbuf_reset(sbuf);
  ASSERT(strbuf_readline(sbuf, fh));
  strbuf_chomp(sbuf);
  ASSERT(strcmp(sbuf->b,"x")==0);

  strbuf_reset(sbuf);
  ASSERT(strbuf_readline_nonempty(sbuf, fh));
  strbuf_chomp(sbuf);
  ASSERT(strcmp(sbuf->b,"y")==0);
  ASSERT_VALID(sbuf);

  strbuf_reset(sbuf);
  ASSERT(strbuf_readline(sbuf, fh));
  strbuf_chomp(sbuf);
  ASSERT(strcmp(sbuf->b,"")==0);

  strbuf_reset(sbuf);
  ASSERT(strbuf_readline_nonempty(sbuf, fh));
  strbuf_chomp(sbuf);
  ASSERT(strcmp(sbuf->b,"z")==0);

  strbuf_reset(sbuf);
  ASSERT(strbuf_readline_nonempty(sbuf, fh) == 0);
  ASSERT(sbuf->end == 0);

  strbuf_reset(sbuf);
  ASSERT(strbuf_readline(sbuf,fh) == 0);
  ASSERT(sbuf->end == 0);

  ASSERT_VALID(sbuf);
  strbuf_free(sbuf);
  fclose(fh);

  SUITE_END();
}

// test trim, ltrim, rtrim
// trim removes whitespace (isspace(c)) from both sides of a str
void _test_trim(const char *str, const char *ans)
{
  StrBuf *sbuf = strbuf_create(str);
  strbuf_trim(sbuf);
  ASSERT(strcmp(sbuf->b, ans) == 0);
  ASSERT(sbuf->end == strlen(sbuf->b));
  ASSERT(sbuf->end < sbuf->size);
  strbuf_free(sbuf);
}
void _test_trim2(const char *str, const char *alphabet, const char *ans,
                     void (*trim)(StrBuf *sbuf, const char *list))
{
  StrBuf *sbuf = strbuf_create(str);
  trim(sbuf, alphabet);
  ASSERT(strcmp(sbuf->b, ans) == 0);
  ASSERT(sbuf->end == strlen(sbuf->b));
  ASSERT(sbuf->end < sbuf->size);
  strbuf_free(sbuf);
}
void test_trim()
{
  SUITE_START("trim");

  // Trim whitespace from either side
  _test_trim("","");
  _test_trim("   ","");
  _test_trim("\r\n","");
  _test_trim("\r \n\t","");
  _test_trim(":\r\n\t.",":\r\n\t.");
  _test_trim("\r \n \t.:",".:");
  _test_trim(".:\r \n \t",".:");
  _test_trim(" abcdefghi\r\njklmn opqrst\tu\nvwxyz",
             "abcdefghi\r\njklmn opqrst\tu\nvwxyz");
  _test_trim(" abcdefghi\r\njklmn opqrst\tu\nvwxyz\n",
             "abcdefghi\r\njklmn opqrst\tu\nvwxyz");
  _test_trim("abcdefghi\r\njklmn opqrst\tu\nvwxyz",
             "abcdefghi\r\njklmn opqrst\tu\nvwxyz");

  // Trim a given alphabet from the left hand side
  _test_trim2("","","", strbuf_ltrim);
  _test_trim2("","abc","", strbuf_ltrim);
  _test_trim2("zabc","abc","zabc", strbuf_ltrim);
  _test_trim2("abacbz","abc","z", strbuf_ltrim);
  _test_trim2("ab:c\nadzb:d\n asdf","abc : \n","dzb:d\n asdf", strbuf_ltrim);

  // Trim a given alphabet from the right hand side
  _test_trim2("","","", strbuf_rtrim);
  _test_trim2("","abc","", strbuf_rtrim);
  _test_trim2("abcz","abc","abcz", strbuf_rtrim);
  _test_trim2("zabacb","abc","z", strbuf_rtrim);
  _test_trim2("ab:c:d\n asdfacb:\n  a","abc : \n", "ab:c:d\n asdf", strbuf_rtrim);

  SUITE_END();
}

void test_safe_ncpy()
{
  SUITE_START("using string_safe_ncpy()");

  const char input[] = "I'm sorry Dave I can't do that";
  char out[100];

  memset(out, 1, sizeof(out));
  string_safe_ncpy(out, input, 0);
  ASSERT(out[0] == 1); // haven't changed out
  string_safe_ncpy(out, input, 1);
  ASSERT(out[0] == 0); // added NULL byte to out[0]
  ASSERT(out[1] == 1); // haven't changed out[1]
  string_safe_ncpy(out, input, 2);
  ASSERT(out[0] == 'I'); // copied one byte to out[0]
  ASSERT(out[1] == 0); // added NULL byte to out[1]
  ASSERT(out[2] == 1); // haven't changed out[2]

  // Copy whole string
  memset(out, 1, sizeof(out));
  string_safe_ncpy(out, input, sizeof(out));
  ASSERT(strcmp(out,input) == 0);
  ASSERT(out[strlen(out)+1] == 1);

  // Copy whole string
  memset(out, 1, sizeof(out));
  string_safe_ncpy(out, input, strlen(input)+1);
  ASSERT(strcmp(out,input) == 0); // copied string+null byte
  ASSERT(out[strlen(out)+1] == 1); // Haven't modified output buffer after string

  // Copy whole string except last byte
  memset(out, 1, sizeof(out));
  string_safe_ncpy(out, input, strlen(input));
  ASSERT(strncmp(out,input,strlen(input)-1) == 0);
  ASSERT(out[strlen(input)-1] == 0); // copied all but last byte
  ASSERT(out[strlen(out)+1] == 1); // Haven't modified output buffer after string

  // Test copying empty string
  memset(out, 1, sizeof(out));
  string_safe_ncpy(out, "", 0);
  ASSERT(out[0] == 1); // haven't changed out
  string_safe_ncpy(out, "", 1);
  ASSERT(out[0] == '\0'); // add NULL byte to out[0]
  ASSERT(out[1] == 1); // haven't changed out[1]
  memset(out, 1, sizeof(out));
  string_safe_ncpy(out, "", 2);
  ASSERT(out[0] == '\0'); // add NULL byte to out[0]
  ASSERT(out[1] == 1); // haven't changed out[1]
  ASSERT(out[2] == 1); // haven't changed out[2]
  ASSERT(out[3] == 1); // haven't changed out[3]

  SUITE_END();
}

void test_split_str()
{
  SUITE_START("using string_split_str()");

  const char input[] = "I'm sorry Dave I can't do that...";
  char buf[100];

  char *ptrs[10];
  size_t i, n;

  memset(ptrs, 0, sizeof(ptrs));
  memcpy(buf, input, sizeof(input));
  n = string_split_str(buf, ' ', ptrs, 3);
  ASSERT(n == 7);
  ASSERT(strcmp(ptrs[0],"I'm") == 0);
  ASSERT(strcmp(ptrs[1],"sorry") == 0);
  ASSERT(strcmp(ptrs[2],"Dave") == 0);

  for(i = 3; i < sizeof(ptrs) / sizeof(ptrs[0]); i++)
    ASSERT(ptrs[i] == NULL);

  memset(ptrs, 0, sizeof(ptrs));
  memcpy(buf, input, sizeof(input));
  n = string_split_str(buf, ' ', ptrs, 10);
  ASSERT(n == 7);
  ASSERT(strcmp(ptrs[0],"I'm") == 0);
  ASSERT(strcmp(ptrs[1],"sorry") == 0);
  ASSERT(strcmp(ptrs[2],"Dave") == 0);
  ASSERT(strcmp(ptrs[3],"I") == 0);
  ASSERT(strcmp(ptrs[4],"can't") == 0);
  ASSERT(strcmp(ptrs[5],"do") == 0);
  ASSERT(strcmp(ptrs[6],"that...") == 0);

  for(i = 7; i < sizeof(ptrs) / sizeof(ptrs[0]); i++)
    ASSERT(ptrs[i] == NULL);

  // Test empty strs
  memset(ptrs, 0, sizeof(ptrs));
  strcpy(buf, "");
  n = string_split_str(buf, ' ', ptrs, 10);
  ASSERT(n == 0);

  for(i = 0; i < sizeof(ptrs) / sizeof(ptrs[0]); i++)
    ASSERT(ptrs[i] == NULL);

  // Test separators at beginning and end
  memset(ptrs, 0, sizeof(ptrs));
  strcpy(buf, ":");
  n = string_split_str(buf, ':', ptrs, 10);
  ASSERT(n == 2);
  ASSERT(strcmp(ptrs[0],"") == 0);
  ASSERT(strcmp(ptrs[1],"") == 0);

  for(i = 2; i < sizeof(ptrs) / sizeof(ptrs[0]); i++)
    ASSERT(ptrs[i] == NULL);

  // Test separators at beginning and end
  memset(ptrs, 0, sizeof(ptrs));
  strcpy(buf, "::and::this::");
  n = string_split_str(buf, ':', ptrs, 10);
  ASSERT(n == 7);
  ASSERT(strcmp(ptrs[0],"") == 0);
  ASSERT(strcmp(ptrs[1],"") == 0);
  ASSERT(strcmp(ptrs[2],"and") == 0);
  ASSERT(strcmp(ptrs[3],"") == 0);
  ASSERT(strcmp(ptrs[4],"this") == 0);
  ASSERT(strcmp(ptrs[5],"") == 0);
  ASSERT(strcmp(ptrs[6],"") == 0);

  for(i = 7; i < sizeof(ptrs) / sizeof(ptrs[0]); i++)
    ASSERT(ptrs[i] == NULL);

  SUITE_END();
}

/* Non-function tests */

void test_sscanf()
{
  SUITE_START("using sscanf");
  StrBuf *sbuf = strbuf_new(10);
  const char *input = "I'm sorry Dave I can't do that";
  
  strbuf_ensure_capacity(sbuf, strlen(input));
  sscanf(input, "I'm sorry %s I can't do that", sbuf->b);
  sbuf->end = strlen(sbuf->b);

  ASSERT(strcmp(sbuf->b, "Dave") == 0);

  strbuf_free(sbuf);
  SUITE_END();
}

/* Old tests */

void test_all_whitespace_old()
{
  printf("Test string_is_all_whitespace:\n");
  const char *str = "  \tasdf";
  printf("string_is_all_whitespace('%s'): %i\n", str, string_is_all_whitespace(str));
  str = "  \t ";
  printf("string_is_all_whitespace('%s'): %i\n", str, string_is_all_whitespace(str));
}

void _test_split_old(const char *split, const char *txt)
{
  char **results;
  
  printf("split '%s' by '%s': (", txt, split);
  
  long count = string_split(split, txt, &results);

  if(count > 0)
  {
    printf("'%s'", results[0]);
    free(results[0]);
  
    int i;
    for(i = 1; i < count; i++)
    {
      printf(", '%s'", results[i]);
      free(results[i]);
    }

    free(results);
  }

  printf(")\n");
}

void test_split_old()
{
  _test_split_old("/", "a/b");
  _test_split_old("/", "/");
  _test_split_old("/", "/b");
  _test_split_old("/", "a/");
  _test_split_old("/", "asdf");
  _test_split_old("/", "");
  _test_split_old("", "asdf");
  _test_split_old("", "");
}


int main()
{
  test_roundup2pow();

  test_buffers();
  test_buffered_reading();

  test_clone();
  test_reset();
  test_resize();

  test_get_set_char();
  test_set();
  test_as_str();
  test_append();
  test_append_int();
  test_chomp();
  test_trim();
  test_reverse();
  test_substr();
  test_change_case();

  test_copy();
  test_insert();
  test_overwrite();
  test_delete();

  test_sprintf();
  test_sprintf_at();
  test_sprintf_noterm();

  test_read_gzfile();
  test_read_file();
  test_read_nonempty();

  test_sscanf();

  // string_ functions
  test_safe_ncpy();
  test_split_str();

  printf("\n");
  TEST_STATS();

  printf("\n THE END.\n");

  // Old tests
  // test_all_whitespace_old();
  // test_split_old();
  
  return total_tests_failed ? 1 : 0;
}
