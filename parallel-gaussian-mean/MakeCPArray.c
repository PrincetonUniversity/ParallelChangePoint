/*************************************************************************
MakeCPArray.c  -  Make arrays of change points from specified tree
*************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include "changepoint.h"
 
/* Make array of change points and confidence intervals from tree with root at
change point node cp inorder i.e. lowest to highest */
void MakeCPArray(struct changepoint *p, int **cps,int **cpsl, int **cpsr, int *Ncp) {
    if (p != NULL) {
        MakeCPArray(p->left,cps,cpsl,cpsr, Ncp);
        (*Ncp)++;
        (*cps) = realloc((*cps),((*Ncp)+1) * sizeof(int));
        (*cpsl) = realloc((*cpsl),((*Ncp)+1) * sizeof(int));
        (*cpsr) = realloc((*cpsr),((*Ncp)+1) * sizeof(int));
        (*cps)[*Ncp] = p->i;
        (*cpsl)[*Ncp] = p->il;
        (*cpsr)[*Ncp] = p->ir;
        MakeCPArray(p->right,cps,cpsl,cpsr,Ncp);
    }
}

// Return size of tree with root at change point node p
void SizeCPArray(struct changepoint *p, int *N) {
    if(p!=NULL){
        (*N)++;
        SizeCPArray(p->left,N);
        SizeCPArray(p->right,N);
    }
}
