#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#include "./CommonFiles/protos.h"

#include "./CommonFiles/graph.c"
#include "./CommonFiles/ilup.c"
#include "./CommonFiles/linked_list.c"
#include "./CommonFiles/matrix.c"
#include "./CommonFiles/rcm.c"

// FUNCAO DE TEMPO
double get_time (){
	struct timeval tv; gettimeofday(&tv, 0);
	return (double)(tv.tv_sec * 100.0 + tv.tv_usec / 10000.0);
}

// ESTRUTURA VECTOR
typedef struct Vector {
    /* how many elements? */
    unsigned int size;
    /* the actual vector */
    double *v;
} Vector;

// MULTIPLICACAO MATRIZ-VETOR CSR
void matrix_vector_multiply_CSR(MAT* A, Vector b, Vector result){
    int i, j, l_begin, l_end;
    int n = b.size;

    double *AA = A->AA;
    int *IA = A->IA;
    int *JA = A->JA;
    double *r = result.v;
    double *bv = b.v;

    for(i = 0; i < n; i++)
    {
        r[i] = 0.0;

        l_begin = IA[i];
        l_end = IA[i+1];

        for (j = l_begin; j < l_end; j++)
        {
            r[i] += AA[j] * bv[JA[j]];
        }
    }
}

// FUNCAO DE CRIACAO DE UM VETOR
Vector BuildVector(unsigned int s) {
    Vector v;
    v.size = s;
    
    if (0 < s) {
        v.v = (double*) calloc(s, sizeof(double));
    }

    return(v);
}

// GERACAO DO VETOR COM VALOR
Vector BuildVectorWithValue(double size, double value){
    Vector vector = BuildVector(size);
    for(int i; i < vector.size; i++){
        vector.v[i] = value;
    }
    return vector;
}

// NORMA EUCLIDIANA
double EuclideanNorm(Vector vector) {
    double norm = 0;

    if (0 < vector.size) {
        for (int i = 0; i < vector.size; i++) {
            norm += pow(vector.v[i], 2);
        }
    } else {
        return 0.0;
    }
    
	return sqrt(norm);
}

// NORMALIZACAO DO VETOR
void ScaleVector(Vector vector, double value) {
    unsigned int v_size = vector.size, i;

    double *v = vector.v;

    for (i = 0; i < v_size; i++) {
        v[i] *= value;
    }
}

// PRODUTO VETORIAL ENTRE DOIS VETORES
Vector CrossProduct(Vector a, Vector b) {

    /* build the new vector */
    Vector newVector = BuildVector(3);

    newVector.v[0] = a.v[1]*b.v[2] - a.v[2]*b.v[1];
    newVector.v[1] = a.v[2]*b.v[0] - a.v[0]*b.v[2];
    newVector.v[2] = a.v[0]*b.v[1] - a.v[1]*b.v[0];

    return newVector;
}

// PRODUTO ESCALAR ENTRE DOIS VETORES
double InnerProduct(Vector a, Vector b) {
    unsigned int i, size = a.size;
    double result = 0.0f;

    for (i = 0; i < size; i++) {
        result += a.v[i] * b.v[i];
    }
    return(result);
}

// FUNCAO QUE CALCULA A MEDIA COM OS VALORES DO VETOR
double media(Vector a) {
	unsigned int i, size = a.size;
	double m, soma = 0.0;
	
	for (i = 0; i < size; i++) {
		soma = soma + a.v[i];
	}
	m = soma/size;
	return m;
}

// FUNCAO QUE DELETA O VETOR
void DeleteVector(Vector vector) {
    free(vector.v);
    return;
}

Vector MultiplicaTermoVec(Vector vector, double value) {
	Vector solucao = BuildVector(vector.size);
		
	for (int j = 0; j < vector.size; j++) {
		solucao.v[j] = vector.v[j] * value;
	}
	
	return solucao;
}

// CALCULO DA SOMA ENTRE VECTOR
Vector calculo_soma(Vector a, Vector b) {
	
	unsigned i = a.size;
	int j;
	
	double *av = a.v;
	double *bv = b.v;
	
	Vector solucao = BuildVector(i);
	
	for (j = 0;j < i;j++) {
		solucao.v[j] = av[j] + bv[j];
	}
		
	return solucao;
}

// CALCULO DA SUBTRACAO ENTRE VECTOR
Vector calculo_sub(Vector a, Vector b) {
	unsigned i = a.size;
	int j;
	
	Vector solucao = BuildVector(i);
	for (j = 0;j < i; j++) {
		solucao.v[j] = a.v[j] - b.v[j];
	}
	return solucao;
}

// FUNCAO DE MULTIPLICACAO MATRIZ x VETOR PARA O METODO LCD
Vector matrix_vector_multiply_CSR_LCD(MAT* A, Vector b)
{
    int i, j, l_begin, l_end;
    int n = b.size;

    double *AA = A->AA;
    int *IA = A->IA;
    int *JA = A->JA;
    double *bv = b.v;

	Vector solucao = BuildVector(n);

    for(i = 0; i < n; i++)
    {
        solucao.v[i] = 0.0;

        l_begin = IA[i];
        l_end = IA[i+1];

        for (j = l_begin; j < l_end; j++){
            solucao.v[i] += AA[j] * bv[JA[j]];
        }
    }

    return solucao;
}

// ESTRUTURA DA SOLUCAO
typedef struct Solution{
    /* the output vector */
    Vector x;

    /* how many iterations */
    unsigned int iterations;

    /* spent time */
    double time;

    /* vector with rho infos */
    Vector rhos;

    /* the rhos history vector size */
    unsigned int rho_last_index;
} Solution;

// FUNCAO PRINT PARA OBSERVACAO DOS VALORES OBTIDOS
void print_vetor_solucao(Vector x) {
	
	int unsigned i;
	int a = x.size;
	
	double *bv = x.v;
	
	for (i = 0;i < a; i++) {
		printf("x(%d): %10.5f \n",i, bv[i]);
	}
	return;	
}

// FUNCAO QUE DELETA A SOLUCAO
void delete_solution(Solution s)
{
    /* remove the solution vector */
    DeleteVector(s.x);

    /* remove the rho history vector */
    DeleteVector(s.rhos);
}

/****************  FUNCAO GMRES ******************/
Solution gmres_solver(MAT *A, Vector b, double tol, unsigned int kmax, unsigned int lmax) {
    /* the common helpers */
    int i, iplus1, j, jplus1, k, kmax1, iter, iter_gmres;
    unsigned int n = b.size;
    double rho, r, tmp;
    double *hv, *prev_v, *next_v;
    double hji, hj1i;

    kmax1 = kmax + 1;

    /* syntatic sugar */
    /* access the vector b */
    double *bv = b.v;

    /* get the initial epsilon value */
    double epsilon = min(tol, tol * EuclideanNorm(b));
    if (0 == epsilon){
        epsilon = tol;
    }

    /* the solution */
    Solution sol;

    /* build the temporary x vector */
    /* this vector will be constantly */
    /* updated until we find the solution */
    sol.x = BuildVector(n);

    /* build the rho history vector */
    sol.rhos = BuildVector(lmax);

    //direct access:
    double* rhov = sol.rhos.v;

    /* reset th rho vector valid elements */
    unsigned int rho_last_index = 0;

    /* the direct access  */
    double *x0 = sol.x.v;

    /* allocate the c and s array */
    double *c = (double*) calloc(kmax, sizeof(double));
    double *s = (double*) calloc(kmax, sizeof(double));

    /* allocate the y array */
    double *y = (double*) calloc(kmax1, sizeof(double));

    /* allocate the error array, starts with zero value */
    double *e = (double*) malloc(kmax1*sizeof(double));

    /* build the u vector array */
    Vector u[kmax1];

    /* build the h vector array */
    Vector h[kmax1];

    /* allocate each u and h vector inside the array'*/
    for (i = 0; i < kmax1; ++i)
    {
        u[i] = BuildVector(n);
        h[i] = BuildVector(kmax1);
    }

    /* get the residual direct access */
    double *rv = u[0].v;

    /* the GMRES main outside loop */
    /* we are going to break this loop */
    /* if we get a tiny little error */
    for (iter = 0; iter < lmax; ++iter)
    {
        /* first let's find the residual/error vector */
        /* start */

        matrix_vector_multiply_CSR(A, sol.x, u[0]);

        /* let's remember: r = b - Ax */
        /* but we got now just the Ax, so ... */
        for (j = 0; j < n; ++j) {
            /* subtract = b - Ax */
            rv[j] = bv[j] - rv[j];
        }

        /* end */

        /* we need the euclidean norm (i.e the scalar error value) */
        rho = EuclideanNorm(u[0]);

        if (rho != 0.0) {

            /* let's normalize the residual vector */
            ScaleVector(u[0], 1.0/rho);

        } else {

            /* we got a trivial solution */
            for (i = 0; i < n; i++)
            {
                x0[i] = 0.0;
            }

            iter += 1;

            break;
        }

        /* update the rho value, the current error */
        rhov[rho_last_index] = rho;

        /* update the last index */
        rho_last_index += 1;

        /* update the rho value */
        e[0] = rho;
        /* reset the error */
        for (j = 1; j < kmax1; ++j)
        {
            e[j] = 0.0;
        }

        /* reset the h matrix */
        for (j = 0; j < kmax1; ++j)
        {
            hv = h[j].v;
            for (k = 0; k < kmax1; ++k)
            {
                hv[k] = 0.0;
            }
        }

        /* reset the i counter */
        i = 0;

        /* the internal loop, for restart purpose */
        while (rho > epsilon && i < kmax)
        {
            /* set the iplus value */
            iplus1 = i + 1;

            /* get the next direction vector */
            matrix_vector_multiply_CSR(A, u[i], u[iplus1]);

            /* Gram-Schmidt process */
            /* GRAM-SCHMIDT begin */
            /* update the h vector */
            hv = h[i].v;

            /* get the next vector direct access */
            next_v = u[iplus1].v;

            for (j = 0; j < iplus1; ++j)
            {
                /* get the prev vector direct access */
                prev_v = u[j].v;

                tmp = InnerProduct(u[j], u[iplus1]);
                hv[j] = tmp;

                /* update the current direction vector */
                for (k = 0; k < n; k++)
                {
                    next_v[k] -= tmp * prev_v[k];
                }
            }

            /* get the euclidean norm of the last direction vector */
            tmp = EuclideanNorm(u[iplus1]);

            /* update the next h vector */
            hv[iplus1] = tmp;

            if (0 == tmp)
            {
                kmax = i;
                i += 1;
                break;
            }

            /* normalize the direction vector */
            ScaleVector(u[iplus1], 1.0/tmp);

            /* GRAM-SCHMIDT end */

            /* QR algorithm */
            if (0 < i)
            {
                for (j = 0, jplus1 = 1; j < i; ++j, ++jplus1)
                {
                    hji = c[j]*hv[j] + s[j]*hv[jplus1];
                    hj1i = -s[j]*hv[j] + c[j]*hv[jplus1];
                    hv[j] = hji;
                    hv[jplus1] = hj1i;
                }

            }

            /* update the residual value */
            r = sqrt(hv[i]*hv[i] + hv[iplus1]*hv[iplus1]);

            /* update the cosine value */
            c[i] = hv[i]/r;

            /* update the sine value */
            s[i] = hv[iplus1]/r;

            /* update the current position inside the h vector */
            hv[i] = r;

            /* update the next position inside the h vector */
            hv[iplus1] = 0.0;

            /* rotate the error  */
            e[iplus1] = -s[i]*e[i];
            e[i] *= c[i];

            /* set the new rho value */
            rho = fabs(e[iplus1]);

            /* update the i counter */
            i += 1;
        }

        /* get the iteration counter */
        iter_gmres = i - 1;

        /* update the y value */
        y[iter_gmres] = e[iter_gmres]/h[iter_gmres].v[iter_gmres];

        for (i = iter_gmres - 1; 0 <= i; --i)
        {
            /* update the h vector direct access */
            hv = h[i].v;

            for (j = iter_gmres; j > i; --j)
            {
                e[i] -= h[j].v[i] * y[j];
            }
            y[i] = e[i]/hv[i];

        }

        iter_gmres += 1;
        for (j = 0; j < n; ++j)
        {
            for (k = 0; k < iter_gmres; ++k)
            {
                /* update the solution vector */
                /* it's so tricky! */
                x0[j] += y[k] * u[k].v[j];
            }

        }

        if (rho < epsilon)
        {
            /* update the rho value, the current error */
            rhov[rho_last_index] = rho;

            /* update the last index */
            rho_last_index += 1;

            break;
        }

    }
	
	if (kmax*lmax <= iter_gmres*lmax)
	{
		iter_gmres = 0;
	}
    /* update the iteration counter*/
    sol.iterations = iter*kmax + iter_gmres;

    /* get the rho last index */
    sol.rho_last_index = rho_last_index;

    free(c);
    free(s);
    free(y);
    free(e);
    //DeleteVector(vRho_tmp);

    /* remove the h and u vectors */
    for (i = 0; i < kmax1; ++i)
    {
        DeleteVector(u[i]);
        DeleteVector(h[i]);
    }

    return sol;
}

/**************** FUNCAO LCD ****************/
Solution lcd_solver(MAT *A, Vector b, double tol, unsigned int kmax, unsigned int lmax){

    int n = b.size;
    
    Vector q[kmax], p[kmax];

    for (int i = 0; i < kmax; i++) {
		q[i] = BuildVector(n);
        p[i] = BuildVector(n);
	}
    
    Vector res = b;
    p[0] = res;

    double rnorm = EuclideanNorm(b);
    double epson = tol * rnorm;     
    
    int itertotal = 0; 
    int l = 0;
    int niter = 0;
    
    Solution sol;
    sol.x = BuildVector(n);

    while((l <= lmax)&&(rnorm > epson)){
        niter = 0;
        q[0] = matrix_vector_multiply_CSR_LCD(A, p[0]);

        for(int i = 0; i < kmax; i++){
            double alpha = InnerProduct(p[i], res) / InnerProduct(p[i], q[i]);

            sol.x = calculo_soma(sol.x, MultiplicaTermoVec(p[i], alpha));
            res = calculo_sub(res, MultiplicaTermoVec(q[i], alpha));
            rnorm = EuclideanNorm(res);
            itertotal++;

            if(rnorm < epson){
                break;
            }

            p[i+1] = res;
            q[i+1] = matrix_vector_multiply_CSR_LCD(A, p[i+1]);

            for(int j = 0; j <= i; j++){
                double beta = -(InnerProduct(q[i+1], p[j]) / InnerProduct(q[j], p[j])); 
                p[i+1] = calculo_soma(p[i+1], MultiplicaTermoVec(p[j], beta));
    			q[i+1] = calculo_soma(q[i+1], MultiplicaTermoVec(q[j], beta));
            }

            niter++;
        }

        p[0] = res;
        l++;
    }

    sol.iterations = l * niter;
   	return sol;
}

int main (int argc, char* argv[]){
	double time;
  	int kmax = 50 , lmax = 10000;
  	double tol = 1e-1;
  
	if (argc != 2)
	{
		printf("\n Erro! Sem arquivo da matriz (.mtx)"); 
		printf("\n Modo de usar: ./program <nome_da_matriz> Saindo... [main]\n\n");
		return 0;
	}
	
	MAT *A = (MAT*) malloc (sizeof(MAT));						
	MATRIX_readCSR (A,argv[1]);
	
	Vector ones = BuildVectorWithValue(A->m, 1.0);
	Vector b = BuildVector(A->m);
	
	matrix_vector_multiply_CSR(A, ones, b);
	
	Solution sol;
	int escolha = 0;
	
	printf("\n k = %d \n lmax = %d", kmax, lmax);

//	printf("\n --> Digite o valor de k e l <-- \n");
//	scanf("%i %i", &kmax, &lmax);
		
//	printf("\n Escolha uma opcao: \n");
//	printf("\n 1. GMRES: \n");
//	printf("\n 2. LCD: \n");
		
//	scanf("%i", &escolha);
		
		// SOLUCAO GMRES
//	if (escolha == 1)
//	{
		time = get_time();
		sol = gmres_solver(A, b, tol, kmax, lmax);
		time = (get_time() - time)/100.0;
		
		printf("\n ###############################");
		printf("\n\n Solucao GMRES:");
		printf("\n Tempo de solucao: %.6f s", time);
		printf("\n Iteracoes: %4d", sol.iterations);
		printf("\n Solucao: %7.6f", media(sol.x));
		printf("\n");
			
		//delete_solution(sol);
//	}
				
		// SOLUCAO LCD
//	if (escolha == 2)
//	{
		time = get_time();
		sol = lcd_solver(A, b, tol, kmax, lmax);
		time = (get_time() - time)/100.0;
			
		printf("\n ###############################");
		printf("\n\n Solucao LCD:");
		printf("\n Tempo de solucao: %.6f s", time);
		printf("\n Iteracoes: %4d", sol.iterations);
		printf("\n Solucao: %7.6f", media(sol.x));
		printf("\n");
			
		//delete_solution(sol);
//	}
	

	MATRIX_clean(A);

	return 0;
}
