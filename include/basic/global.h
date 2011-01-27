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


#define LENGTH_FILENAME 300
#define VERSION 1
#define SUBVERSION 0
#define SUBSUBVERSION 3
boolean DEBUG;




#endif /* GLOBAL_H_ */
