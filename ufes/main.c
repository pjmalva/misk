#include <stdio.h>
#include "algoritmos-c/CommonFiles/protos.h"
#include "gmres.c"


void main(int argc, char* argv[]){
    
    MAT *A = (MAT*) malloc (sizeof(MAT));						
	MATRIX_readCSR(A,argv[1]);
    
    *F = ;
    
    *V = ;

    GMRES(&A, &F, &V);
}


