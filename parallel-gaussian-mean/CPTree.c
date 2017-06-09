/********************************************************************
CPTree.c  -  Binary tree structure to store change points found
*************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include "changepoint.h"
#include <math.h>

struct changepoint* rightRotate();
struct changepoint* leftRotate();
struct changepoint* maxValue();
int getBal();
int height();

// Return maximum of a and b
int max(int a, int b) {
	return (a > b)? a : b;
}

// Return height of cp node p
int height(struct changepoint *p) {
	if (p == NULL)
		return 0;
	return p->height;
}

// Get balance of cp node p
int getBal(struct changepoint *p) {
	return height(p->right) - height(p->left);
}

// Right rotate cp node p
struct changepoint* rightRotate(struct changepoint *p) {
	struct changepoint* newP = p->left;
	struct changepoint* newC = newP->right;

	newP->right = p;
	p->left = newC;
	p->height = max(height(p->left), height(p->right)) + 1;
	newP->height = max(height(newP->left), height(newP->right)) + 1;
	return newP;
}

// Left Rotate cp node p
struct changepoint* leftRotate(struct changepoint *p) {
	struct changepoint* newP = p->right;
	struct changepoint* newC = newP->left;

	newP->left = p;
	p->right = newC;
	p->height = max(height(p->left), height(p->right)) + 1;
	newP->height = max(height(newP->left), height(newP->right)) + 1;
	return newP;
}

// Return min value of the subtree with root at p
struct changepoint * minValue(struct changepoint * p) {
	struct changepoint* current = p;

	while(current->left != NULL)
		current = current->left;
	return current;
}

// Add cp k with boundaries kl and kr to tree with root p
struct changepoint * AddCPNode(int k, int kl, int kr, struct changepoint *p, 
	int *Ncp, enum bool *h) {
	int ol, or, bal;

	if (p == NULL) {                         // Not in tree, insert
		p = (struct changepoint *) malloc(sizeof(struct changepoint));
		*h = true;
		p->i = k;        
		p->il = kl;  
		p->ir = kr;
		p->left = NULL;
		p->right = NULL;
		p->height = 1;
		(*Ncp)++;
	}
	else if ( k < p->i )
		p->left = AddCPNode(k,kl,kr,p->left, Ncp, h);
	else if ( k > p->i ) 
		p->right = AddCPNode(k,kl,kr,p->right,Ncp, h);
	else // Node alredy exists
		*h = false;

	// New node inserted, rebalance tree
	if(*h) {
		p->height = 1+max(height(p->left), height(p->right));
		bal = getBal(p);

		if(bal < -1 && k < p->left->i)
			return rightRotate(p);
		else if(bal < -1 && k > p->left->i) {
			p->left = leftRotate(p->left);
			return rightRotate(p);
		}

		if(bal > 1 && k > p->right->i)
			return leftRotate(p);
		else if(bal > 1 && k < p->right->i) {
			p->right = rightRotate(p->right);
			return leftRotate(p);
		}
	}
	return p;
}

// Delete node with cp k in tree with root p
struct changepoint *DelCPNode(int k, struct changepoint*p, enum bool *h) {
	struct changepoint * temp;
	int bal;

	// Node does not exist in tree
	if(p == NULL) {
		fprintf(stderr, "Changepoint at %d does not exist in tree\n", k);
		exit(1);
	}

	if(k < p->i)
		p->left = DelCPNode(k, p->left, h);
	else if(k > p->i)
		p->right = DelCPNode(k, p->right, h);
	else {

		// One child or no child
		if((p->left == NULL) || (p->right == NULL) ) {
			temp = p->left ? p->left : p->right;

			// No child
			if(temp == NULL) {
				temp = p;
				p = NULL;
			}
			// Copy contents of not null child
			else
				*p = *temp;

			free(temp);
			*h = true;
		}
		else {
			// Take leftmost node (smallest node) in right subtree of p
			temp = minValue(p->right);

			p->i = temp->i;
			p->ir = temp->ir;
			p->il = temp->il;
			p->right = DelCPNode(temp->i, p->right, h);
		}
	}

	if(p == NULL)
		return p;

	// Update heights of current node
	p->height = 1+max(height(p->left), height(p->right));
	bal = getBal(p);

	// Deleted from left node
	if(bal < -1 && getBal(p->left) < 0)
		return rightRotate(p);
	else if (bal < -1 && getBal(p->left) > 0) {
		p->left = leftRotate(p->left);
		return rightRotate(p);
	}

	// Deleted from right node
	else if(bal > 1 && getBal(p->right) > 0)
		return leftRotate(p);
	else if(bal > 1 && getBal(p->right) < 0) {
		p->right = rightRotate(p->right);
		return leftRotate(p);
	}
	return p;
}
