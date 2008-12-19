/*
 * splay.c
 *
 *  Created on: Sep 19, 2008
 *      Author: mariocaccamo
 */

#include <recycleBin.h>
#include <splay.h>


#define CHUNKSIZE 10000



typedef struct SplayNode_st SplayNode;

static RecycleBin *treeMemory = NULL;


static SplayNode *allocate_splay_node()
{
	if (treeMemory == NULL)
		treeMemory = newRecycleBin(sizeof(SplayNode), CHUNKSIZE);

	return (SplayNode *) allocatePointer(treeMemory);
}

static void deallocate_splay_node(SplayNode * node)
{
	deallocatePointer(treeMemory, node);
}


SplayTree *new_splay_tree()
{
	return NULL;
}

void destroy_all_splay_trees()
{
	destroyRecycleBin(treeMemory);
	treeMemory = NULL;
}

void destroy_splay_tree(SplayTree * tree)
{
	if (tree == NULL)
		return;

	destroy_splay_tree(tree->left);
	destroy_splay_tree(tree->right);
	deallocate_splay_node(tree);
}

/* This function can be called only if K2 has a left child */
/* Perform a rotate between a node (K2) and its left child */
/* Update heights, then return new root */

static SplayNode * single_rotate_with_left(SplayNode * k2)
{
	SplayNode *k1;

	if (k1->left != NULL){
		k1 = k2->left;
		k2->left = k1->right;
		k1->right = k2;
	}
	else{
		k1 = k2;
	}

	return k1;		/* New root */
}

/* This function can be called only if K1 has a right child */
/* Perform a rotate between a node (K1) and its right child */
/* Update heights, then return new root */

static SplayNode * single_rotate_with_right(SplayNode * k1)
{
	SplayNode *k2;

	if (k1->right != NULL){
		k2 = k1->right;
		k1->right = k2->left;
		k2->left = k1;
	}
	else{
		k2 = k1;
	}

	return k2;		/* New root */
}


/* Top-down splay procedure, */
/* not requiring element to be in tree */

static SplayTree * splay(Element * element, SplayTree * tree)
{
	SplayNode header;
	SplayNode * left_tree_max, * right_tree_min;

	if (tree == NULL)
		return NULL;

	header.left = header.right = NULL;
	left_tree_max = right_tree_min = &header;

	//printf("T1\n");
	while ( ! element_equal(element,tree->element)) {
		//printf("T2\n");
		if ( element_smaller(element,tree->element)) {
			//printf("T3\n");
			if (tree->left == NULL)
				break;
			if (element_smaller(element,tree->left->element)){
				///printf("T3a\n");
				tree = single_rotate_with_left(tree);
			}
			if (tree->left == NULL)
				break;
			/* Link right */
			right_tree_min->left = tree;
			right_tree_min = tree;
			tree = tree->left;
		} else {
			//printf("T4\n");
			if (tree->right == NULL)
				break;
			if (element_bigger(element,tree->right->element))
				tree = single_rotate_with_right(tree);
			if (tree->right == NULL)
				break;
			/* Link left */
			left_tree_max->right = tree;
			left_tree_max = tree;
			tree = tree->right;
		}
	}			/* while element != tree->element */
	//printf("T5\n");
	/* Reassemble */
	left_tree_max->right = tree->left;
	right_tree_min->left = tree->right;
	tree->left  = header.right;
	tree->right = header.left;

	return tree; //returns the root I guess
}

boolean find_in_tree(Element * element, SplayTree ** tree)
{
	*tree = splay(element, *tree);
	return element_equal(element,(*tree)->element);
}

//apply function if element is in tree, add element otherwise
boolean splay_apply_or_insert(Element *element, void (*f)(Element*,Element*),SplayTree ** tree)
{
	SplayNode *new_node;
	boolean found = false;
	//printf("S1\n");

	if (*tree == NULL) {
		//printf("S2\n");
		new_node = allocate_splay_node();
		new_node->element = element;
		new_node->left = new_node->right = NULL;
		*tree = new_node;
	}
	else{ //tree is not null
		//printf("S3\n");
		*tree = splay(element, *tree);

		if (element_smaller(element, (*tree)->element)) {
	    	new_node = allocate_splay_node();
	    	new_node->element = element;
	    	new_node->left = (*tree)->left;
	    	new_node->right = *tree;
	    	(*tree)->left = NULL;
	    	*tree = new_node;
	    } else if (element_smaller((*tree)->element,element)) {
	    	new_node = allocate_splay_node();
	    	new_node->element = element;
	    	new_node->right = (*tree)->right;
	    	new_node->left = *tree;
	    	(*tree)->right = NULL;
	    	*tree = new_node;
			}
			else{
				found = true;
				//element in tree
				f(element,(*tree)->element);
			}
	}
	//printf("S4\n");
	return found;
}

//this consumes lots of stack doesn't it?
void splay_inorder_traverse_apply(void(*f)(Element*),SplayTree * tree){

	if (tree !=NULL){
		splay_inorder_traverse_apply(f,tree->left);
		f(tree->element);
		splay_inorder_traverse_apply(f,tree->right);
	}
}

