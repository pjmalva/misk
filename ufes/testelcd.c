#include "program.c"

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