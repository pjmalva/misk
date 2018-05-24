#include<stdio.h>
#include <stdlib.h>

/* A number of global data structures and variables relating to
the mesh and its partition are declared here. */

typedef struct
{
	double Pos[2];
	int Local;
	int Type;
	int Global;
}NodeType;

typedef struct
{
	int Vertex[3];
}ElementType;

NodeType *Node;
ElementType *Element;
int Nodes, Elements, IntNodes;

/* The following function calculates the inner product of two
   vectors distributed in a prescribed manner. */

double InnerProduct(double *A1, double *B1)
{
	int I;
	double IP = 0.0;

	for (I=0; I<IntNodes; I++)
		IP += A1[I] * B1[I];
	return(IP);
}

void SolveSystem(double **A, double *x, double *b, int n)
{
	int i,j;
	double soma;
	
	x[n] = b[n]/A[n][n];
	for (i=n-1;i>=0;i--)
	{
		soma = 0;
		for (j=i+1;j<=n;j++)
			soma += A[j][i]*x[j];	
		x[i] = (1./A[i][i])*(b[i]-soma);
	}
	
}

		
int GMRES(double **Ap, double *Fp, double *Vp)
{
	double **up, **h, *Beta, *c, *s, *e, *y, Eps, inv_aux, r, res, n_u, t, Tol = 1.0e-8;
	int i, j, I, J, iter, lmax = 50000, k=10;

	Beta = (double*) malloc(k*sizeof(double));
	c = (double*) malloc(k*sizeof(double));
	s = (double*) malloc(k*sizeof(double));
	e = (double*) malloc((k+1)*sizeof(double));
	y = (double*) malloc(k*sizeof(double));
	up = (double**) malloc((k+1)*sizeof(double));
	for (I=0;I<k+1;I++)
		up[I] = (double*) malloc(IntNodes*sizeof(double));
	h = (double**) malloc(k*sizeof(double*));
	for (I=0;I<k;I++)
		h[I] = (double*) malloc((I+2)*sizeof(double));
	
	Eps = Tol*sqrt(InnerProduct(Fp,Fp));
	printf("Eps = %lf\n",Eps);
	res = Eps + 1.0;
	iter = 0;
	while (iter<lmax && res>Eps)
	{
		for (I=0; I<IntNodes; I++)
		{
			for (J=0, up[0][I] = Fp[I]; J<IntNodes; J++)
				up[0][I] -= Ap[I][J] * Vp[J];
		}
		e[0] = sqrt(InnerProduct(up[0], up[0]));
		
		inv_aux = 1./e[0];

		for (I=0; I<IntNodes; I++)
			up[0][I] = inv_aux*up[0][I];

		i = 0;
		while (i<k && res>Eps )
		{
				
			for (I=0; I<IntNodes; I++)
				for (J=0, up[i+1][I]=0.; J<IntNodes; J++)
					up[i+1][I] += Ap[I][J] * up[i][J];
		
			
			for (j=0;j<i+1;j++)
			{
				Beta[j] = InnerProduct(up[i+1],up[j]);	
				
				for (I=0; I<IntNodes; I++)
					up[i+1][I] -= Beta[j] * up[j][I];
				
				h[i][j] = Beta[j];
			
			}
			
			n_u = sqrt((InnerProduct(up[i+1],up[i+1])));
			h[i][i+1] = n_u;
			inv_aux = 1./n_u;
			
			for (I=0; I<IntNodes; I++)
				up[i+1][I] = inv_aux*up[i+1][I];

			
			for (j=0;j<i;j++)
			{
				t = h[i][j];
				h[i][j] = c[j]*t + s[j]*h[i][j+1];	
				h[i][j+1] = -s[j]*t + c[j]*h[i][j+1];
			}
	
			r = sqrt(h[i][i]*h[i][i] + h[i][i+1]*h[i][i+1]);
			inv_aux = 1./r;
			c[i] = inv_aux*h[i][i];
			s[i] = inv_aux*h[i][i+1];
			h[i][i] = r;
			h[i][i+1]= 0.0;
			e[i+1] = -s[i]*e[i];
			e[i] = c[i]*e[i];
			res = fabs(e[i+1]);		
			i++;
		}
		i--;
		
		SolveSystem(h,y,e,i);	
		
		for (j=0;j<i+1;j++)
			for (I=0;I<IntNodes;I++)
				Vp[I] += y[j]*up[j][I];
		
		iter++;		
		printf("iteracoes = %d, erro [%d]= %.15lf\n", iter, i+1, res); 	

	}
	return(0);
}	