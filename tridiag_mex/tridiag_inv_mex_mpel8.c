/* mex file for tridiagonal solver 2013-11-09 
 assuming block diagonal with identical, real tridiagonal blocks
 so rhs is formulated as NxM matrix */


#include "mex.h"
#include "def/defs-env-local.h"
#include "pthread.h"

//#if !defined(Need_tridiag_inv_mex_gateway)
//#include "tridiag_inv,def.h"
//#endif

#define Usage "usage error. see above"
#define NUM_THREADS 20 // number of cores // 4 for iv1, 2 for vega

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

/*typedef struct 
{
	float *yy;
	int Ny; // etc??
} tridiag_inv_worker_args; */


struct thread_data
{
	int thread_id;
	int block_size; 
	double *subdiag_ptr;
	double *diagvals_ptr;
	double *supdiag_ptr;
	double *rhs_ptr;
	double *out_ptr;
};

struct thread_data thread_data_array[NUM_THREADS];

//static sof tridiag_inv_thr(void *threadarg)
void *tridiag_inv_thr(void *threadarg)
{
	//tridiag_inv_worker_args args;
	int taskid;
	int N;
	double *a;
	double *b;
	double *c;
       	double *d;
       	double *x;
	double *new_c;
	double *new_d;
	struct thread_data *my_data;
	my_data = (struct thread_data *) threadarg;
	taskid = my_data -> thread_id;
	N = my_data -> block_size;
	a = my_data -> subdiag_ptr;
	b = my_data -> diagvals_ptr;
	c = my_data -> supdiag_ptr;
	d = my_data -> rhs_ptr;
	x = my_data -> out_ptr;

	// need to malloc for x_real, x_imag, and then free them ??
	new_c = (double *) calloc (N - 1, sizeof(double));
	new_d = (double *) calloc (N, sizeof(double));

	*new_c = *c / *b; //new_c[0] = c[0]/b[0];
	*new_d = *d / *b; //new_d[0] = d[0]/b[0];

	int ii;
	double a_prev, new_c_prev;
	for (ii = 1; ii <= N-2; ii++) {
		/*
		a_prev = *(a + ii - 1);
		new_c_prev = *(new_c + ii - 1);
		*(new_c + ii) = *(c + ii) / (*(b + ii) - new_c_prev * a_prev);
		*(new_d + ii) = (*(d + ii) - *(new_d + ii - 1) * a_prev) / (*(b + ii) - new_c_prev * a_prev);
		*/
		a_prev = a[ii - 1];
		new_c_prev = new_c[ii - 1];
		new_c[ii] = c[ii] / (b[ii] - new_c_prev * a_prev);
		new_d[ii] = (d[ii] - new_d[ii - 1] * a_prev) / (b[ii] - new_c_prev * a_prev);
	}
	//*(new_d + N - 1) = (*(d + N - 1) - *(new_d + N - 2) * (*(a + N - 1))) / (*(b + N - 1) - *(new_c + N - 2) * (*(a + N - 2)));
	new_d[N - 1] = (d[N - 1] - new_d[N - 2] * (a[N - 2])) / (b[N - 1] - new_c[N - 2] * (a[N - 2]));

	//*(x + N - 1) = *(new_d + N - 1);
	x[N - 1] = new_d[N - 1];
	for (ii = N-2; ii >= 0; ii--) {
		//*(x + ii) = *(new_d + ii) - *(new_c + ii) * (*(x + ii + 1));
		x[ii] = new_d[ii] - new_c[ii] * x[ii + 1];
	}

//	printf("task id: %d \n", taskid);

//	pthread_mutex_lock(&mutexout);
	
//	pthread_mutex_unlock(&mutexout);
	free(new_c);
	free(new_d);

	pthread_exit(NULL);
}

/*
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
} */


// currently assume all inputs real
// wrapper function for thread
static sof tridiag_inv_mex_thr(
double *subdiag_ptr, double *diagvals_ptr, double *supdiag_ptr, double *rhs_real_ptr, double *rhs_imag_ptr, mwSize block_size, mwSize nblocks, double *out_real_ptr, double *out_imag_ptr)
{
	//tridiag_inv_worker_args args;

	int ii;
	int big_N; // total number of entries, N*nblocks
	int rc;
	long t;
    int block_ndx;
	pthread_attr_t attr;
	
//	printf("block size : %d \n", block_size);
//    printf("nblocks : %d, NUM_THREADS: %d, nblocks/NUM_THREADS: %d \n", nblocks, NUM_THREADS, (nblocks/ NUM_THREADS));
	pthread_t threads[NUM_THREADS];

	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    
//    printf("\n");
    // do all real values first
    for (int th_rep = 0; th_rep <= nblocks/NUM_THREADS; th_rep++) {
        for (int th_id = 0; th_id < NUM_THREADS; th_id++) {
            block_ndx = th_rep * NUM_THREADS + th_id;
//            printf("th_rep: %d, th_id: %d, block_ndx: %d \n", th_rep, th_id, block_ndx);
            if (block_ndx <= nblocks - 1) {
                thread_data_array[th_id].thread_id = th_id;
                thread_data_array[th_id].block_size = block_size;
                thread_data_array[th_id].subdiag_ptr = subdiag_ptr;// + th_id * block_size;
                thread_data_array[th_id].diagvals_ptr = diagvals_ptr;// + th_id * block_size;
                thread_data_array[th_id].supdiag_ptr = supdiag_ptr;// + th_id * block_size;
                thread_data_array[th_id].rhs_ptr = rhs_real_ptr + block_ndx * block_size;
                thread_data_array[th_id].out_ptr = out_real_ptr + block_ndx * block_size;
                rc = pthread_create(&threads[th_id], &attr, tridiag_inv_thr, (void *) &(thread_data_array[th_id]));
            }
        }
//        printf("\n");
        for (int th_id = 0; th_id < NUM_THREADS; th_id++) {
            block_ndx = th_rep * NUM_THREADS + th_id;
            if (block_ndx <= nblocks - 1) {
                rc = pthread_join(threads[th_id], NULL);
                if (rc) {
                    printf("ERROR; return code from pthread_join() is %d\n", rc);
                    exit(-1);
                }
//                printf("done with thread ndx %d ", block_ndx);
            }
        }
//        printf("\n");
//        printf("done with thread blocks %d - %d \n", th_rep*(nblocks/NUM_THREADS), block_ndx);
    }
    // do all complex values next
    for (int th_rep = 0; th_rep <= nblocks/NUM_THREADS; th_rep++) {
        for (int th_id = 0; th_id < NUM_THREADS; th_id++) {
            block_ndx = th_rep * NUM_THREADS + th_id;
            //            printf("th_rep: %d, th_id: %d, block_ndx: %d \n", th_rep, th_id, block_ndx);
            if (block_ndx <= nblocks - 1) {
                thread_data_array[th_id].thread_id = th_id;
                thread_data_array[th_id].block_size = block_size;
                thread_data_array[th_id].subdiag_ptr = subdiag_ptr;// + th_id * block_size;
                thread_data_array[th_id].diagvals_ptr = diagvals_ptr;// + th_id * block_size;
                thread_data_array[th_id].supdiag_ptr = supdiag_ptr;// + th_id * block_size;
                thread_data_array[th_id].rhs_ptr = rhs_imag_ptr + block_ndx * block_size;
                thread_data_array[th_id].out_ptr = out_imag_ptr + block_ndx * block_size;
                rc = pthread_create(&threads[th_id], &attr, tridiag_inv_thr, (void *) &(thread_data_array[th_id]));
            }
        }
        //        printf("\n");
        for (int th_id = 0; th_id < NUM_THREADS; th_id++) {
            block_ndx = th_rep * NUM_THREADS + th_id;
            if (block_ndx <= nblocks - 1) {
                rc = pthread_join(threads[th_id], NULL);
                if (rc) {
                    printf("ERROR; return code from pthread_join() is %d\n", rc);
                    exit(-1);
                }
                //                printf("done with thread ndx %d ", block_ndx);
            }
        }
        //        printf("\n");
        //        printf("done with thread blocks %d - %d \n", th_rep*(nblocks/NUM_THREADS), block_ndx);
    }
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
    double *rhs_imag;

    size_t N;                   /* size of tridiag matrix */
    size_t M;                   /* numcols of rhs matrix */
    double *x_real;                  /* output */
    double *x_imag;                  /* output */
    
	if (nrhs != 4 ) { // hard coding :(
		tridiag_inv_mex_help();
		//Call(mxu_arg, (nrhs, prhs))
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
    
    if (mxIsComplex(prhs[3])) {
        //printf("rhs is complex \n");
        rhs_imag = mxGetPi(prhs[3]);
        //out_imag_ptr = mxGetPi(plhs[0]);
        plhs[0] = mxCreateDoubleMatrix((mwSize)N, (mwSize) M, mxCOMPLEX);
    } else {
        rhs_imag = NULL;
        plhs[0] = mxCreateDoubleMatrix((mwSize)N, (mwSize) M, mxREAL);
        
    }
    x_real = mxGetPr(plhs[0]);
    x_imag = mxGetPi(plhs[0]); // NULL when rhs is real
    
   // check_types_and_sizes(sub, diag, sup, rhs, N, M);
    
	//Call(tridiag_inv_mex_thr, (nlhs, plhs, nrhs, prhs)); // why is using "call" better?
    tridiag_inv_mex_thr(sub, diag, sup, rhs, rhs_imag, N, M, x_real, x_imag);
	
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


