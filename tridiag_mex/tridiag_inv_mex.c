/* tridiag_inv_mex.c
 Matlab mex file for tridiagonal solver 2013-11-09
 T * x = y
 assumes real tridiagonal matrix T, rhs y is formulated as NxM matrix, computes M columns of x in parallel
*/

#include "mex.h"
#include "defs-env.h"
#include "jf,mex,def.h"
//#include "def/def,mexarg.h" // in above
#include "pthread.h"

#define Usage "usage error. see above"
#define NUM_THREADS 2 // number of cores // 4 for iv1, 2 for vega
#define VERBOSE false


static void tridiag_inv_mex_help(void)
{
    printf("\n\
	Usage for tridiag_inv_mex: \n\
	output = tridiag_inv_mex(subdiag, diagvals, supdiag, argum) \n\
	\n\
           T * x = y \n\
           \n\
           subdiag: (single, real) [N-1 1] -1st subdiagonal values of T \n\
           diagvals: (single, real) [N 1] diagonal values of T \n\
           supdiag: (single, real) [N-1 1] 1st diagonal values of T \n\
           argum: (single) [N M] rhs of inverse problem, y \n\
	\n");
}

struct thread_data
{
    int thread_id;
	int block_size;
    int num_blocks_for_me;
	float *subdiag_ptr;
	float *diagvals_ptr;
	float *supdiag_ptr;
	float *rhsr_ptr;
    float *rhsi_ptr;
	float *outr_ptr;
    float *outi_ptr;
};

struct thread_data thread_data_array[NUM_THREADS];

// tridiag_inv()
// tridiag solver over one block
static sof tridiag_inv(float *a, float *b, float *c, float *d, int N, float *x, float *new_c, float *new_d)
{

	new_c[0] = c[0] / b[0];
	new_d[0] = d[0] / b[0];

	int ii;
	float a_prev, new_c_prev;
	for (ii = 1; ii <= N-2; ii++) {
		a_prev = a[ii - 1];
		new_c_prev = new_c[ii - 1];
		new_c[ii] = c[ii] / (b[ii] - new_c_prev * a_prev);
		new_d[ii] = (d[ii] - new_d[ii - 1] * a_prev) / (b[ii] - new_c_prev * a_prev);
	}
	new_d[N - 1] = (d[N - 1] - new_d[N - 2] * (a[N - 2])) / (b[N - 1] - new_c[N - 2] * (a[N - 2]));

	x[N - 1] = new_d[N - 1];
	for (ii = N-2; ii >= 0; ii--) {
        x[ii] = new_d[ii] - new_c[ii] * x[ii + 1];
	}
	Ok
}

// tridiag_inv_loop_thr()
// each thread loops over tridiag_inv
static sof tridiag_inv_loop_thr(void *threadarg)
{
    struct thread_data *my_data;
    int taskid;
    int num_runs;
    int N;
    float *a;
    float *b;
    float *c;
    float *dr;
    float *xr;
    float *di;
    float *xi;
    float *new_c;
    float *new_d;
    
    my_data = (struct thread_data *) threadarg;
    taskid = my_data -> thread_id;
    N = my_data -> block_size;
    a = my_data -> subdiag_ptr;
    b = my_data -> diagvals_ptr;
    c = my_data -> supdiag_ptr;
    dr = my_data -> rhsr_ptr;
    xr = my_data -> outr_ptr;
    di = my_data -> rhsi_ptr;
    xi = my_data -> outi_ptr;
    num_runs = my_data -> num_blocks_for_me;
    new_c = (float *) calloc (N - 1, sizeof(float));
    new_d = (float *) calloc (N, sizeof(float));
    
#if VERBOSE
    printf("taskid: %d, N: %d, inside looper func with num runs %d \n", taskid, N, num_runs);
#endif
    for (int ii = 0; ii < num_runs; ii++) {
#if VERBOSE
        printf(" looping over index %d /%d \n", ii, num_runs);
        
        printf("address of rhs: %d, incr by %d,  new addr: %d, float size: %d \n", dr, N, dr+N, sizeof(float));
        printf("rhs vals: \n");
/*        for (int jj = 0; jj < N; jj++) {
            printf("rhs[%d] = %d \n", jj, dr[jj]);
        } */
#endif
        tridiag_inv(a, b, c, dr, N, xr, new_c, new_d);
        dr += N;
        xr += N;
        if (xi != NULL) {
            tridiag_inv(a, b, c, di, N, xi, new_c, new_d);
            di += N;
            xi += N;
        }
    }
    free(new_c);
    free(new_d);
    
    pthread_exit(NULL);
    Ok
}

// check_types_and_sizes()
static sof check_types_and_sizes(
Const mxArray *prhs[]
)
{
    int Nsub;
    int Msub;
    int Ndiag;
    int Mdiag;
    int Nsup;
    int Msup;
    int Nrhs;
    
    Nsub = mxGetM(prhs[0]);
    Msub = mxGetN(prhs[0]);
    Ndiag = mxGetM(prhs[1]);
    Mdiag = mxGetN(prhs[1]);
    Nsup = mxGetM(prhs[2]);
    Msup = mxGetN(prhs[2]);
    Nrhs = mxGetM(prhs[3]);
    
	int pass = 1;
	if (!((Nsub == Nrhs - 1) && (Msub == 1)) && !((Nsub == 1) && (Msub == Nrhs - 1))) {
		printf("subdiag size [%d %d] does not match rhs length of %d \n", Nsub, Msub, Nrhs);
		pass = 0;
	}
    if (!((Ndiag == Nrhs) && (Mdiag == 1)) && !((Ndiag == 1) && (Mdiag == Nrhs))) {
        printf("diag size [%d %d] does not match rhs length of %d \n", Ndiag, Mdiag, Nrhs);
        pass = 0;
    }
    if (!((Nsup == Nrhs - 1) && (Msup == 1)) && !((Nsup == 1) && (Msup == Nrhs - 1))) {
        printf("supdiag size [%d %d] does not match rhs length of %d \n", Nsup, Msup, Nrhs);
        pass = 0;
    }
    if (mxIsComplex(prhs[0])) {
        printf("subdiag cannot be complex \n");
        pass = 0;
    }
    if (!mxIsClass(prhs[0], "single")) {
        printf("subdiag must be single \n");
        pass = 0;
    }
    if (mxIsComplex(prhs[1])) {
        printf("diag cannot be complex \n");
        pass = 0;
    }
    if (!mxIsClass(prhs[1], "single")) {
        printf("diag must be single \n");
        pass = 0;
    }
    if (mxIsComplex(prhs[2])) {
        printf("supdiag cannot be complex \n");
        pass = 0;
    }
    if (!mxIsClass(prhs[2], "single")) {
        printf("supdiag must be single \n");
        pass = 0;
    }
    if (!mxIsClass(prhs[3], "single")) {
        printf("rhs must be single \n");
        pass = 0;
    }
	if (!pass) {
		tridiag_inv_mex_help();
        Fail(Usage)
	}
	printf("\n");
    Ok
}


// currently assume all inputs real
// wrapper function for thread
static sof tridiag_inv_mex_thr(
float *subdiag_ptr, float *diagvals_ptr, float *supdiag_ptr, float *rhs_real_ptr, float *rhs_imag_ptr, mwSize block_size, mwSize nblocks, float *out_real_ptr, float *out_imag_ptr)
{
	int rc;
    int block_ndx;
    int blocks_per_thread[NUM_THREADS];
    int cum_blocks[NUM_THREADS + 1];
    int remainder_blocks;
    
    pthread_attr_t attr;
	pthread_t threads[NUM_THREADS];
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    remainder_blocks = nblocks - (nblocks/NUM_THREADS)*NUM_THREADS;
#if VERBOSE
    printf("NUM_THREADS: %d, nblocks: %d \n", NUM_THREADS, nblocks);
    printf("remainder_blocks: %d \n", remainder_blocks);
#endif
    
    cum_blocks[0] = 0;
    for (int th_ndx = 0; th_ndx < NUM_THREADS; th_ndx++) {
        blocks_per_thread[th_ndx] = nblocks/NUM_THREADS;
        if (th_ndx < remainder_blocks) {
            blocks_per_thread[th_ndx]++;
        }
        cum_blocks[th_ndx + 1] = cum_blocks[th_ndx] + blocks_per_thread[th_ndx];
#if VERBOSE
        printf("cum_blocks_per_thread for %d : %d \n", th_ndx+1, cum_blocks[th_ndx+1]);
#endif
    }
#if VERBOSE
    
    printf("rhs vals: \n");
    for (int jj = 0; jj < block_size; jj++) {
        printf("rhs[%d] = %d \n", jj, *(rhs_real_ptr + jj) );
    }
#endif
    
    for (int th_id = 0; th_id < NUM_THREADS; th_id++) {
        thread_data_array[th_id].thread_id = th_id;
        thread_data_array[th_id].block_size = block_size;
        thread_data_array[th_id].subdiag_ptr = subdiag_ptr;
        thread_data_array[th_id].diagvals_ptr = diagvals_ptr;
        thread_data_array[th_id].supdiag_ptr = supdiag_ptr;
        thread_data_array[th_id].rhsr_ptr = rhs_real_ptr + cum_blocks[th_id] * block_size;
        thread_data_array[th_id].outr_ptr = out_real_ptr + cum_blocks[th_id] * block_size;
        if (rhs_imag_ptr != NULL) {
            thread_data_array[th_id].rhsi_ptr = rhs_imag_ptr + cum_blocks[th_id] * block_size;
            thread_data_array[th_id].outi_ptr = out_imag_ptr + cum_blocks[th_id] * block_size;
        } else {
            thread_data_array[th_id].rhsi_ptr = NULL;
            thread_data_array[th_id].outi_ptr = NULL;
        }
        thread_data_array[th_id].num_blocks_for_me = blocks_per_thread[th_id];
#if VERBOSE
        printf(" th_id: %d, rhs ptr: %d, curr ptr: %d \n \n", th_id, rhs_real_ptr, thread_data_array[th_id].rhsr_ptr);
#endif
        rc = pthread_create(&threads[th_id], &attr, (void *) &tridiag_inv_loop_thr, (void *) &(thread_data_array[th_id]));
    }
    Ok
}


// intermediate GateWay routine 
static sof tridiag_inv_mex_gw(
int nlhs, mxArray *plhs[],
int nrhs, Const mxArray *prhs[])
{
    float *sub;              /* input subdiagonal */
    float *diag;               /* 1xN input diagonal */
    float *sup;
    float *rhs;
    float *rhs_imag;
    size_t N;                   /* size of tridiag matrix */
    size_t M;                   /* numcols of rhs matrix */
    float *x_real;                  /* output */
    float *x_imag;                  /* output */
    
	if (nrhs != 4 ) {
		tridiag_inv_mex_help();
		Fail(Usage)
	}

    check_types_and_sizes(prhs);
    
    sub = (float *) mxGetData(prhs[0]);
    diag = (float *)  mxGetData(prhs[1]);
    sup =  (float *) mxGetData(prhs[2]);
    rhs = (float *) mxGetData(prhs[3]);
    N = mxGetM(prhs[1]) + mxGetN(prhs[1]) - 1; // agnostic to col or row
    M = mxGetN(prhs[3]);

#if VERBOSE
    printf("M: %d, N: %d \n", M, N);
    printf(" address of rhs in gw: %d \n", rhs);
    printf("rhs vals: \n");
    for (int jj = 0; jj < N; jj++) {
        printf("rhs[%d] = %f \n", jj, *(rhs + jj) );
    }
#endif
    
    if (mxIsComplex(prhs[3])) {
        rhs_imag = (float *) mxGetImagData(prhs[3]);
        plhs[0] = mxCreateNumericMatrix((int) N, (int) M, mxSINGLE_CLASS, mxCOMPLEX);
    } else {
        rhs_imag = NULL;
        plhs[0] = mxCreateNumericMatrix((int) N, (int) M, mxSINGLE_CLASS, mxREAL);
    }
    x_real = (float *)mxGetData(plhs[0]);
    x_imag = (float *)mxGetImagData(plhs[0]); // NULL when rhs is real
    
    if ((rhs_imag == NULL) ^ (x_imag == NULL)) {
        printf("problem: only one NULL vec, should be both or neither");
    }
    
	Call(tridiag_inv_mex_thr, (sub, diag, sup, rhs, rhs_imag, N, M, x_real, x_imag));
	Ok
}

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



