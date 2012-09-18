
/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo
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

/*
 * global.h
 *
 */

#ifndef GLOBAL_H_
#define GLOBAL_H_

#include <err.h>

typedef signed char boolean;
#ifndef true
#define true 1
#define false 0
#endif

typedef enum{
  forward = 0,
  reverse = 1
} Orientation;


#define MAX_READ_NAME_LEN 300
#define LENGTH_FILENAME 300
#define VERSION 1
#define SUBVERSION 0
#define SUBSUBVERSION 5
#define SUBSUBSUBVERSION 12
boolean DEBUG;


boolean test_file_existence(char* file);

void die(const char* fmt, ...)
  __attribute__ ((format(printf, 1, 2)))
  __attribute__ ((noreturn));

// 

#endif /* GLOBAL_H_ */
