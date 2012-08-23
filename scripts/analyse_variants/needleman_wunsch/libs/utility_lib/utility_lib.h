/*
 utility_lib.h
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

#ifndef UTILITY_HEADER_SEEN
#define UTILITY_HEADER_SEEN

//
// Integers
//

// parse an int value - exit with some error message if invalid value
long parse_int(char* c, char* err);

// Returns 1 if whole string is int, 0 otherwise
char parse_entire_int(const char *str, int *result);
char parse_entire_uint(const char *str, unsigned int *result);
char parse_entire_long(const char *str, long *result);
char parse_entire_ulong(const char *str, unsigned long *result);
char parse_entire_longlong(const char *str, long long *result);
char parse_entire_ulonglong(const char *str, unsigned long long *result);

// Sort comparison function on int
int int_cmp(const void *aa, const void *bb);

//
// Binary
//

// Convert an int/long to a string of zeros and ones
char* int_to_binary(const int x);
char* long_to_binary(const long x);

//
// Strings
//

// Checks if anything is piping in
int stdin_is_ready();

char string_is_all_whitespace(const char* s);

char* next_nonwhitespace(const char* s);

// Move pointer to first nonwhitespace position
// Remove trailing whitespace characters (trim with '\0')
char* trim(char* str);

long count_strchr(const char* str, const int c);

#endif
