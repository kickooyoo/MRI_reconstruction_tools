/* mex file for tridiagonal solver 2013-11-09 */

#include "mex.h"
#include "matrix.h"
//#include "def,type.h"
#include "defs-env.h"
#include "jf,mex,def.h" 
#include "jf,thread1.h"
#include "pthread.h"

//#if !defined(Need_tridiag_inv_mex_gateway)
//#include "tridiag_inv,def.h"
//#endif

#define Usage "usage error. see above"
#define NUM_THREADS 4 // number of cores

//pthread_mutex_t mutexout; // global var for locking

//#if defined(Mmex)

static void tridiag_inv_mex_help(void)
{
	printf("\n\
	Usage for tridiag_inv_mex: \n\
	output = tridiag_inv_mex(subdiag,diagvals,supdiag,argum) \n\
	\n\
	subdiag has length n-1	\n\
	expect tridiag matrix to be real \n\
	\n");
}

// does one col of x at a time
static sof tridiag_inv_mex(
double *a, double *b, double *c, double *d, mwSize N, double *x)
{
	double a_prev, new_c_prev;

	mxArray *new_c; 
	mxArray *new_d;
	//mxArray *x;
	double *new_c_db;
	double *new_d_db;

	int ndims;
	int output_dims[ndims];

    mwSize i;
    for (i=0; i<N; i++){
        x[i] = b[i] + d[i];
    }


	new_c = mxCreateDoubleMatrix(N, 1, mxREAL);// todo: set to zero
	new_c_db = mxGetPr(new_c);
	new_d = mxCreateDoubleMatrix(N, 1, mxREAL);// todo: set to zero
	new_d_db = mxGetPr(new_d);


	// need to do the following for real and imag components of d

	*new_c_db = (*c) / (*b); //new_c[0] = c[0]/b[0];
	*new_d_db = (*d) / (*b); //new_d[0] = d[0]/b[0];
    
	int ii;
	for (ii = 1; ii <= N-2; ii++) {
		a_prev =  *(a + ii -1);
		new_c_prev =  *(new_c_db + ii -1);
		*(new_c_db + ii) = *(c + ii) / (*(b + ii) - new_c_prev * a_prev);
		//new_c[ii] = c[ii]/(b[ii]-new_c[ii-1]*a[ii-1]);	
		*(new_d_db + ii) = (*(d + ii) - *(new_d_db + ii -1) * a_prev) / (*(b + ii) - new_c_prev * a_prev);
		//new_d[ii] = (d[ii]-new_d[ii-1]*a[ii-1])/(b[ii]-new_c[ii-1]*a[ii-1]);
	}
    
    

	*(new_d_db + N - 1) = (*(d + N - 1) - *(new_d_db + N - 2) * (*(a + N -2))) / (*(b + N - 1) - *(new_c_db + N - 2) * (*(a + N - 2)));
	//new_d[N-1] = (d[N-1]-new_d[N-2]*a[N-1])/(b[N-1]-new_c[N-2]*a[N-2]);

	*(x + N - 1) = *(new_d_db + N -1);
    
	//x[N-1] = new_d[N-1];
	for (ii = N-2; ii >= 0; ii--) {
		*(x + ii) = *(new_d_db + ii) - *(new_c_db + ii) * (*(x + ii +1));
		//x[ii] = new_d[ii]-new_c[ii]*x[ii+1];
	}
    
    /*ndims = 2;
	output_dims[0] = N; 
	output_dims[1] = mxGetN(d); // should be 1
	*/ //FIGURE THIS SHIT OUT
	//Call(plhs[0] = mxCreateNumericArray, (ndims, output_dims, mxDOUBLE_CLASS,mxCOMPLEX))

    
	//Call(tridiag_inv, (args)) // why???
	Ok


    
}

static sof check_types_and_sizes(
Const mxArray *subdiag, 
Const mxArray *diagvals, 
Const mxArray *supdiag, 
Const mxArray *rhs, 
Const mxArray *block_size_ptr)
{
	int pass = 1;
	if (!mxIsDouble(subdiag)) {
		printf("subdiag is not double type \n");
		pass = 0;
	}
	if (!mxIsDouble(supdiag)) {
		printf("supdiag is not double type \n");
		pass = 0;
	}
	if (!mxIsDouble(diagvals)) {
		printf("diagvals is not double type \n");
		pass = 0;
	}
	if (!mxIsDouble(rhs)) {
		printf("rhs is not double type \n");
		pass = 0;
	}
	// todo: check lengths using mxGetM
	if (pass) {
		printf(" all are double types \n");
	}
	printf("\n");
    
    Ok
}

// intermediate GateWay routine 
static sof tridiag_inv_mex_gw(
int nlhs, mxArray *plhs[],
int nrhs, Const mxArray *prhs[])
{
    double *sub;              /* input subdiagonal */
    double *diag;               /* 1xN input diagonal */
    double *sup;
    double *rhs;
    size_t N;                   /* size of tridiag matrix */
    size_t M;                   /* numcols of rhs matrix */
    double *x;                  /* output */
//   double *new_c; /* debug output */
//   double *new_d; /* debug output */
    
	if (nrhs != 4 ) { // hard coding :(
		tridiag_inv_mex_help();
		//Call(mxu_arg, (nrhs, prhs))
        printf("Fail in gw \n");
		Fail(Usage)
	}

	// todo: check lengths of vectors
	//
    
    sub = mxGetPr(prhs[0]);
    diag = mxGetPr(prhs[1]);
    sup = mxGetPr(prhs[2]);
    rhs = mxGetPr(prhs[3]);
    N = mxGetM(prhs[1]); //remember col vec!
    M = mxGetN(prhs[3]);
    
    plhs[0] = mxCreateDoubleMatrix((mwSize)N, (mwSize) M, mxREAL);
    x = mxGetPr(plhs[0]);
//    plhs[1] = mxCreateDoubleMatrix((mwSize)(N-1),1,mxREAL);
//    new_c = mxGetPr(plhs[1]);
//    plhs[2] = mxCreateDoubleMatrix((mwSize)N,1,mxREAL);
//    new_d = mxGetPr(plhs[2]);
    //printf("N is %d \n", N);
    
    int jj;
    for (jj=0; jj<M; jj++){
        tridiag_inv_mex(sub, diag, sup, rhs + N*jj, N, x + N*jj);
    }
	
	Ok
}


//#if defined (Need_tridiag_inv_mex_gateway)
// gateway routine 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (!nlhs && !nrhs) {
		tridiag_inv_mex_help();
		return;
	}
	if (!tridiag_inv_mex_gw(nlhs, plhs, nrhs, prhs))
		mexErrMsgTxt("tridiag_inv_mex");
}

//#endif


