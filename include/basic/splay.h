/*
 * splay.h
 *
 *  Created on: Sep 17, 2008
 *      Author: mariocaccamo
 */


#include "element.h"

#ifndef SPLAY_H_
#define SPLAY_H_

struct SplayNode_st {
	Element * element;
	struct SplayNode_st *left;
	struct SplayNode_st *right;
};

typedef struct SplayNode_st SplayTree;


SplayTree *new_splay_tree();

boolean splay_apply_or_insert(Element *, void (*f)(Element*,Element*),SplayTree **);
void splay_inorder_traverse_apply(void(*f)(Element*),SplayTree *);

#endif /* SPLAY_H_ */
