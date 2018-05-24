#include<stdio.h>
#include <stdlib.h>

double **CriarMatrix(int linhas, int colunas){    
    double **m[linhas][colunas];
    int i, j = 0;
    for(i = 0; i < colunas; i++){
        for(j = 0; j < linhas; j++){
            m[i][j] = (double **) malloc(sizeof(double)); 
        }
    }
    return **m;
}

void PopularDiagonalPrincipal(double **matrix, int linhas, int colunas){
    int i = 0, TOTAL = linhas * colunas, l = 0, c = 0;
    double *ptr = &matrix[0][0];

    for (i = 0; i < TOTAL; i++) {
        if(l == c)
            *ptr = 1;
        else
            *ptr = 0;
        l++;
        c++;
        *ptr++;
    }
}

void ImprimirMatrix(double **matrix, int linhas, int colunas){

    int i = 0, TOTAL = linhas * colunas;

    double *ptr = &matrix[0][0];

    for (i = 0; i < TOTAL; i++) {
        printf("%f ", *(ptr + i));
    }

}

void main(){
    int c = 5;
    int l = 5;

    double **matrix = CriarMatrix(l, c);
    
    printf("%d", (int) sizeof(matrix));

    PopularDiagonalPrincipal(matrix, l, c);
    // ImprimirMatrix(matrix, l, c);
}