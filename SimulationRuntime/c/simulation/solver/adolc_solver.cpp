/*
 * This file is part of OpenModelica.
 *
 * Copyright (c) 1998-2014, Open Source Modelica Consortium (OSMC),
 * c/o Linköpings universitet, Department of Computer and Information Science,
 * SE-58183 Linköping, Sweden.
 *
 * All rights reserved.
 *
 * THIS PROGRAM IS PROVIDED UNDER THE TERMS OF THE BSD NEW LICENSE OR THE
 * GPL VERSION 3 LICENSE OR THE OSMC PUBLIC LICENSE (OSMC-PL) VERSION 1.2.
 * ANY USE, REPRODUCTION OR DISTRIBUTION OF THIS PROGRAM CONSTITUTES
 * RECIPIENT'S ACCEPTANCE OF THE OSMC PUBLIC LICENSE OR THE GPL VERSION 3,
 * ACCORDING TO RECIPIENTS CHOICE.
 *
 * The OpenModelica software and the OSMC (Open Source Modelica Consortium)
 * Public License (OSMC-PL) are obtained from OSMC, either from the above
 * address, from the URLs: http://www.openmodelica.org or
 * http://www.ida.liu.se/projects/OpenModelica, and in the OpenModelica
 * distribution. GNU version 3 is obtained from:
 * http://www.gnu.org/copyleft/gpl.html. The New BSD License is obtained from:
 * http://www.opensource.org/licenses/BSD-3-Clause.
 *
 * This program is distributed WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE, EXCEPT AS
 * EXPRESSLY SET FORTH IN THE BY RECIPIENT SELECTED SUBSIDIARY LICENSE
 * CONDITIONS OF OSMC-PL.
 *
 */

/*! \file adolc_solver.c
 */

#include <string.h>
#include <vector>
#include <adolc/edfclasses.h>
#include <adolc/taping.h>
#include <adolc/interfaces.h>
#include <adolc/adolc_fatalerror.h>
#include <adolc/adalloc.h>
#include <adolc/drivers/drivers.h>
#include <adolc/tapedoc/asciitapes.h>

#include "simulation/solver/adolc_solver.h"

extern "C" {

int dgetrf_(int *m, int *n, double *a, int *lda,
            int *ipiv, int *info);

int dgetrs_(char* trans, int *n, int *nrhs, double *a, int *lda, int *ipiv,
            double *b, int *ldb, int *info);

int dgesv_(int *n, int *nrhs, double *a, int *lda,
           int *ipiv, double *b, int *ldb, int *info);

int dgemv_(char* trans, int* m, int *n, double *alpha, double *A,
           int *lda, double *x, int *incx, double *beta, double *y, int *incy);

int dgemm_(char* transa, char* transb, int* m, int* n, int* k, double* alpha, 
           double* a, int* lda, double* b, int* ldb, double* beta, double* c, 
           int* ldc);
}

class LinearSolverEdf : public EDFobject_v2 {
public:
    LinearSolverEdf() : EDFobject_v2() {}
    virtual ~LinearSolverEdf() {}
    virtual int function(int iArrLen, int *iArr, int nin, int nout, int *insz, double **x, int *outsz, double **y, void* ctx);
    virtual int zos_forward(int iArrLen, int *iArr, int nin, int nout, int *insz, double **x, int *outsz, double **y, void* ctx);
    virtual int fos_forward(int iArrLen, int* iArr, int nin, int nout, int *insz, double **x, double **xp, int *outsz, double **y, double **yp, void *ctx);
    virtual int fov_forward(int iArrLen, int* iArr, int nin, int nout, int *insz, double **x, int ndir, double ***Xp, int *outsz, double **y, double ***Yp, void* ctx);
    virtual int fos_reverse(int iArrLen, int* iArr, int nout, int nin, int *outsz, double **up, int *insz, double **zp, double **x, double **y, void *ctx);
    virtual int fov_reverse(int iArrLen, int* iArr, int nout, int nin, int *outsz, int dir, double ***Up, int *insz, double ***Zp, double **x, double **y, void* ctx);
};

static std::vector<LinearSolverEdf> linSolEdfVec;

int LinearSolverEdf::function(int iArrLen, int *iArr, int nin, int nout, int *insz, double **x, int *outsz, double **y, void* ctx) {

  // nin == 2 and nout == 1
  // iArrlen == 2*nnz

  int nnz = insz[0];
  int nb  = insz[1];
  int nx  = outsz[0];
  int nrhs = 1;
  int info;

  double* A = (double*) calloc(nb*nx, sizeof(double));
  int *ipriv = (int*) calloc(nx, sizeof(int));

  double* b = y[0];
  for(int i=0; i<nb; ++i){
    b[i] = x[1][i];
    //printf("b[%d] = %f\n", i, b[i]);
  }

  for(int i = 0; i<nnz; ++i){
    A[iArr[2*i]+iArr[2*i+1]*nb] = x[0][i];
    //printf("A[%d, %d] = %f\n", iArr[2*i], iArr[2*i+1], A[iArr[2*i]+iArr[2*i+1]*nb]);
  }

  dgesv_(&nx, &nrhs, A, &nx, ipriv, b, &nb, &info);

  /*
  for(int i=0; i<nb; ++i){
      printf("sol[%d] = %f\n", i, b[i]);
  }
  */


  free(A);
  free(ipriv);
  if(info < 0){
    printf("function: Error solving linear system of equations. Argument %d illegal.\n", info);
    return -info;
  }else if (info > 0){
    printf("function: Error solving linear system of equations. System singular in row %d.\n", info);
    return -info;
  }

  return 0;
}

int LinearSolverEdf::zos_forward(int iArrLen, int *iArr, int nin, int nout, int *insz, double **x, int *outsz, double **y, void* ctx) {

  return this->function(iArrLen, iArr, nin, nout, insz, x, outsz, y, ctx);
}


int LinearSolverEdf::fos_forward(int iArrLen, int* iArr, int nin, int nout, int *insz, double **x, double **xp, int *outsz, double **y, double **yp, void *ctx) {


  // x = A^{-1}*b
  function(iArrLen, iArr, nin, nout, insz, x, outsz, y, ctx);

  //x^{.} = A^{-1}*(\dot b - \dot A*x)
  // \dot b -> xp[1]
  // \dot A -> xp[0]

  int nnz = insz[0];
  int nb  = insz[1];
  int nx  = outsz[0];
  int nrhs = 1;
  int info;

  double* A = (double*) calloc(nb*nx, sizeof(double));
  int *ipriv = (int*) calloc(nx, sizeof(int));

  for(int i = 0; i<nnz; ++i){
    A[iArr[2*i]+iArr[2*i+1]*nb] = xp[0][i];
  }

  double* b = yp[0];
  for(int i=0; i<nb; ++i){
    b[i] = xp[1][i];
  }

  char trans = 'N';
  double alpha = -1.0;
  double beta = 1.0;
  int incx = 1;
  dgemv_(&trans, &nb, &nx, &alpha, A, &nx, y[0], &incx, &beta, b, &incx);


  for(int i = 0; i<nnz; ++i){
    A[iArr[2*i]+iArr[2*i+1]*nb] = x[0][i];
  }

  dgesv_(&nx, &nrhs, A, &nx, ipriv, b, &nb, &info);

  free(A);
  free(ipriv);
  if(info < 0){
    printf("fos_forward: Error solving linear system of equations. Argument %d illegal.\n", info);
    return -info;
  }else if (info > 0){
    printf("fos_forward: Error solving linear system of equations. System singular in row %d.\n", info);
    return -info;
  }

  return 0;
}



int LinearSolverEdf::fov_forward(int iArrLen, int* iArr, int nin, int nout, int *insz, double **x, int ndir, double ***Xp, int *outsz, double **y, double ***Yp, void* ctx) {

  // x = A^{-1}*b
  function(iArrLen, iArr, nin, nout, insz, x, outsz, y, ctx);

  // \dot x = A^(-1)*(\dot B - \dot \underscore{A}*x)
  // \dot b -> xp[1]
  // \dot A -> xp[0]

  int nnz = insz[0];
  int nb  = insz[1];
  int nx  = outsz[0];
  int nrhs = ndir;
  int info;

  double* A = (double*) calloc(nb*nx, sizeof(double));
  int *ipriv = (int*) calloc(nx, sizeof(int));
  double** b = Yp[0];

  char trans = 'N';
  double alpha = -1.0;
  double beta = 1.0;
  int incx = 1;
  int incy = ndir;

  for (int k = 0; k < ndir; ++k){

    for(int i = 0; i<nnz; ++i){
      A[iArr[2*i]+iArr[2*i+1]*nb] = Xp[0][i][k];
    }

    for(int i=0; i<nb; ++i){
      b[i][k] = Xp[1][i][k];
    }

    dgemv_(&trans, &nb, &nx, &alpha, A, &nx, y[0], &incx, &beta, &b[0][k], &incy);
  }

  for(int i = 0; i<nnz; ++i){
    A[iArr[2*i]+iArr[2*i+1]*nb] = x[0][i];
  }

  dgesv_(&nx, &nrhs, A, &nx, ipriv, &b[0][0], &nb, &info);

  free(A);
  free(ipriv);
  if(info < 0){
    printf("fov_forward: Error solving linear system of equations. Argument %d illegal.\n", info);
    return -info;
  }else if (info > 0){
    printf("fov_forward: Error solving linear system of equations. System singular in row %d.\n", info);
    return -info;
  }

  return 0;
}


int LinearSolverEdf::fos_reverse(int iArrLen, int* iArr, int nout, int nin, int *outsz, double **up, int *insz, double **zp, double **x, double **y, void *ctx){
  return 0;
}
int LinearSolverEdf::fov_reverse(int iArrLen, int* iArr, int nout, int nin, int *outsz, int dir, double ***Up, int *insz, double ***Zp, double **x, double **y, void* ctx){
  return 0;
}

/* =============== */
/* Non-linear part */
/* =============== */

#include <simulation/solver/newtonIteration.h>

class NonLinearSolverEdf : public EDFobject_v2 {
protected:
    short trace1, trace2, trace3, nexttag;
public:
    DATA_NEWTON* data;
	double **J, **I2;
	int numIterVar;
	char *tmp;

    NonLinearSolverEdf(const char* nlsfbase, short tagstart);
    virtual ~NonLinearSolverEdf();
    virtual int function(int iArrLen, int *iArr, int nin, int nout, int *insz, double **x, int *outsz, double **y, void* ctx);
    virtual int zos_forward(int iArrLen, int *iArr, int nin, int nout, int *insz, double **x, int *outsz, double **y, void* ctx);
    virtual int fos_forward(int iArrLen, int* iArr, int nin, int nout, int *insz, double **x, double **xp, int *outsz, double **y, double **yp, void *ctx);
    virtual int fov_forward(int iArrLen, int* iArr, int nin, int nout, int *insz, double **x, int ndir, double ***Xp, int *outsz, double **y, double ***Yp, void* ctx);
    virtual int fos_reverse(int iArrLen, int* iArr, int nout, int nin, int *outsz, double **up, int *insz, double **zp, double **x, double **y, void *ctx);
    virtual int fov_reverse(int iArrLen, int* iArr, int nout, int nin, int *outsz, int dir, double ***Up, int *insz, double ***Zp, double **x, double **y, void* ctx);
    short get_next_tag() const { return nexttag; }
    friend int wrapper_fvec_newton_adolc(int* n, double* x, double* fvec, void* userdata, int fj);
};

static std::vector<NonLinearSolverEdf> nonLinSolEdfVec;

NonLinearSolverEdf::~NonLinearSolverEdf() {
	free(tmp);
	myfreeI2(numIterVar,I2);
	freeNewtonData(reinterpret_cast<void**>(&data));
}

template <typename Type>
static inline char* populate_dpp_with_contigdata(Type ***const pointer, char *const memory,
                                   int n, int m, Type *const data) {
    char* tmp;
    Type **tmp1; Type *tmp2;
    int i,j;
    tmp = (char*)memory;
    tmp1 = (Type**) memory;
    *pointer = tmp1;
    tmp = (char*)(tmp1+n);
    tmp2 = data;
    for (i=0;i<n;i++) {
        (*pointer)[i] = tmp2;
        tmp2 += m;
    }
    return tmp;
}

int wrapper_fvec_newton_adolc(int* n, double* x, double* fvec, void* userdata, int fj) {
	NonLinearSolverEdf* nledf = reinterpret_cast<NonLinearSolverEdf*>(userdata);
	DATA_NEWTON* dataLocal = nledf->data;
	if (!fj) {
		::zos_forward(nledf->trace1,*n,*n,0,x,fvec);
	} else {
		::fov_forward(nledf->trace1,*n,*n,*n,x,nledf->I2,fvec,nledf->J);
	}
	return 0;
}

NonLinearSolverEdf::NonLinearSolverEdf(const char* nlsfbase, short tagstart) : EDFobject_v2() {
    char *nlsfilename;
    size_t slen = strlen(nlsfbase);
    nlsfilename = new char[slen + 15];

    sprintf(nlsfilename,"%s_1_aat.txt",nlsfbase);
    trace1 = tagstart;
    tagstart = read_ascii_trace(nlsfilename,tagstart);

    /* //debug
    sprintf(nlsfilename,"%s_1_aat_w.txt",nlsfbase);
    write_ascii_trace(nlsfilename, trace1);
    printf("tag 1: %d\n", trace1);
    printTapeStats(stdout, trace1);
     */

    sprintf(nlsfilename,"%s_2_aat.txt",nlsfbase);
    trace2 = tagstart;
    tagstart = read_ascii_trace(nlsfilename,tagstart);

    /* //debug
    sprintf(nlsfilename,"%s_2_aat_w.txt",nlsfbase);
    write_ascii_trace(nlsfilename, trace2);
    printf("tag 2: %d\n", trace2);
    printTapeStats(stdout, trace2);
    */

    sprintf(nlsfilename,"%s_3_aat.txt",nlsfbase);
    trace3 = tagstart;
    tagstart = read_ascii_trace(nlsfilename,tagstart);
    nexttag = tagstart;

    /* //debug
    sprintf(nlsfilename,"%s_3_aat_w.txt",nlsfbase);
    write_ascii_trace(nlsfilename, trace3);
    printf("tag 3: %d\n", trace3);
    printTapeStats(stdout, trace3);
    printf("next tag: %d\n", nexttag);
    */

    delete [] nlsfilename;
}

int NonLinearSolverEdf::function(int iArrLen, int *iArr, int nin, int nout, int *insz, double **x, int *outsz, double **y, void* ctx) {
    // assumption: nin == 1, nout == 2
    //             insz[0] == sizeof input vars x
    //             outsz[0] == sizeof output vars y1
    //             outsz[1] == sizeof output vars y2
    int numouterparams, i;
    double *outerparams;
    numouterparams = alloc_copy_current_params(&outerparams);
    // solver takes input from x, jacobian from trace1 and gives output y2
    size_t stats1[STAT_SIZE];
    tapestats(trace1, stats1);
    int num_resid1 = stats1[NUM_DEPENDENTS];
    int num_indep1 = stats1[NUM_INDEPENDENTS]; // should be == outsz[1]
    int num_param1 = stats1[NUM_PARAM]; // should be numouterparams + insz[0]
    double* allparams1 = (double*) calloc(num_param1,sizeof(double));
    for (i = 0; i < numouterparams; i++)
        allparams1[i] = outerparams[i];
    for (i = 0; i < insz[0]; i++)
        allparams1[numouterparams+i] = x[0][i];
    set_param_vec(trace1, num_param1, allparams1);
    // call solver iteration with y2 and J until convergence
    _omc_newton(wrapper_fvec_newton_adolc,data,reinterpret_cast<void*>(this));
    if (data->info <0){
    	throw FatalError(data->info, "Nonlinear solver failed!", __func__, __FILE__, __LINE__);
    }
    for(i=0;i<outsz[1];i++) {
    	y[1][i] = data->x[i];
    }
    set_param_vec(trace3, num_param1, allparams1);
    ::zos_forward(trace3,outsz[1],outsz[0],0,y[1],y[0]);
    free(allparams1);
    free(outerparams);
    return 0;
}

int NonLinearSolverEdf::zos_forward(int iArrLen, int *iArr, int nin, int nout, int *insz, double **x, int *outsz, double **y, void* ctx) {

  return this->function(iArrLen, iArr, nin, nout, insz, x, outsz, y, ctx);
}


int NonLinearSolverEdf::fos_forward(int iArrLen, int* iArr, int nin, int nout, int *insz, double **x, double **xp, int *outsz, double **y, double **yp, void *ctx) {
    // do everything as in function
    // assumption: nin == 1, nout == 2
    //             insz[0] == sizeof input vars x
    //             outsz[0] == sizeof output vars y1
    //             outsz[1] == sizeof output vars y2
    int numouterparams, i;
    double *outerparams;
    numouterparams = alloc_copy_current_params(&outerparams);
    // solver takes input from x, jacobian from trace1 and gives output y2
    size_t stats1[STAT_SIZE];
    tapestats(trace1, stats1);
    int num_resid1 = stats1[NUM_DEPENDENTS];
    int num_indep1 = stats1[NUM_INDEPENDENTS]; // should be == outsz[1]
    int num_param1 = stats1[NUM_PARAM]; // should be numouterparams + insz[0]
    double* allparams1 = (double*) calloc(num_param1,sizeof(double));
    for (i = 0; i < numouterparams; i++)
        allparams1[i] = outerparams[i];
    for (i = 0; i < insz[0]; i++)
        allparams1[numouterparams+i] = x[0][i];
    set_param_vec(trace1, num_param1, allparams1);
    // call solver iteration with y and J until convergence
    _omc_newton(wrapper_fvec_newton_adolc,data,reinterpret_cast<void*>(this));
    if (data->info <0){
    	throw FatalError(data->info, "Nonlinear solver failed!", __func__, __FILE__, __LINE__);
    }
    for(i=0;i<outsz[1];i++) {
    	y[1][i] = data->x[i];
    }
    double **J1 = myalloc2(outsz[1],outsz[1]);
    jacobian(trace1,outsz[1],outsz[1],y[1],J1);
    double **J3 = myalloc2(outsz[1],outsz[0]);
    set_param_vec(trace3, num_param1, allparams1);
    jacobian(trace3,outsz[0],outsz[1],y[1],J3);
    free(allparams1);

    // finally compute jacobian or resid,y1 wrt y2 after convergence in J
    size_t stats2[STAT_SIZE];
    tapestats(trace2, stats2);
    int num_depen2 = stats2[NUM_DEPENDENTS];// should be == outsz[0] + outsz[1]
    int num_indep2 = stats2[NUM_INDEPENDENTS]; // should be == insz[0]
    int num_param2 = stats2[NUM_PARAM]; // should be numouterparams
    double* allparams2 = (double*) calloc(num_param2,sizeof(double));
    double* depen = (double*) calloc(num_depen2,sizeof(double));
    for (i = 0; i < numouterparams; i++)
        allparams2[i] = outerparams[i];
    for (i = 0; i < outsz[0]; i++)
        allparams2[numouterparams+i] = y[1][i];
    set_param_vec(trace2, num_param2, allparams2);
    // allocate J2 as directional deriv of resid w.r.t. x
    double *J2 = (double*) calloc(num_depen2,sizeof(double));
    ::fos_forward(trace2, num_depen2, num_indep2, 0, x[0], xp[0], depen, J2);
    for (i=0;i<outsz[0];i++) {
        	y[0][i] = depen[outsz[1]+i];
    }
    free(depen);

    // First outsz[1] rows of J2 contain dr/dx, next outsz[0] rows contain dy1/dx
    // solve yp[1] = - J1^{-1} * dr/dx using a linear solver

    int nrhs = 1;
    int info;

    /*
    double* A = (double*) malloc(outsz[1]*outsz[1]*sizeof(double));
    for(int i = 0; i<outsz[1]; ++i){
      for(int j = 0; j<outsz[1]; ++j){
    	  A[i+j*outsz[1]] = J1[j][i];
      }
      //printf("A[%d, %d] = %f\n", iArr[2*i], iArr[2*i+1], A[iArr[2*i]+iArr[2*i+1]*nb]);
    }
    */
    char trans = 'T';
    double* A = &J1[0][0];

    double* b = yp[1];
    for(int i=0; i<outsz[1]; ++i){
      b[i] = -J2[i];
      //printf("b[%d] = %f\n", i, b[i]);
    }
    int *ipriv = (int*) calloc(outsz[1], sizeof(int));
    
    dgetrf_(&outsz[1], &outsz[1], A, &outsz[1], ipriv, &info);
    if (info == 0) {
        dgetrs_(&trans, &outsz[1], &nrhs, A, &outsz[1], ipriv, b, &outsz[1], &info);
    }
    if (info < 0) {
        printf("fos_forward: Error solving linear system of equations. Argument %d illegal.\n", info);
        return -1;
    }
    //dgesv_(&outsz[1], &nrhs, A, &outsz[1], ipriv, b, &outsz[1], &info);

    free(ipriv);
    myfree2(J1);

    //free(A);
    // compute yp[0] = J3*yp[1] + dy1/dx <=> yp[0] = A*yp[1]+b

    A = &J3[0][0];
    b = yp[0];
    for(int i=0; i<outsz[0]; ++i){
      b[i] = J2[outsz[1]+i];
      //printf("b[%d] = %f\n", i, b[i]);
    }
    double alpha = 1.0;
    double beta = 1.0;
    int incx = 1;
    dgemv_(&trans, &outsz[1], &outsz[0], &alpha, A, &outsz[1], yp[1], &incx, &beta, b, &incx);

    myfree2(J3);
    free(J2);
    free(outerparams);
    free(allparams2);
    return 0;
}



int NonLinearSolverEdf::fov_forward(int iArrLen, int* iArr, int nin, int nout, int *insz, double **x, int ndir, double ***Xp, int *outsz, double **y, double ***Yp, void* ctx) {

    // almost same as fos_forward
    // do everything as in function
    // assumption: nin == 1, nout == 2
    //             insz[0] == sizeof input vars x
    //             outsz[0] == sizeof output vars y1
    //             outsz[1] == sizeof output vars y2
    int numouterparams, i;
    double *outerparams;
    numouterparams = alloc_copy_current_params(&outerparams);
    // solver takes input from x, jacobian from trace1 and gives output y2
    size_t stats1[STAT_SIZE];
    tapestats(trace1, stats1);
    int num_resid1 = stats1[NUM_DEPENDENTS];
    int num_indep1 = stats1[NUM_INDEPENDENTS]; // should be == outsz[1]
    int num_param1 = stats1[NUM_PARAM]; // should be numouterparams + insz[0]
    double* allparams1 = (double*) calloc(num_param1,sizeof(double));
    for (i = 0; i < numouterparams; i++)
        allparams1[i] = outerparams[i];
    for (i = 0; i < insz[0]; i++)
        allparams1[numouterparams+i] = x[0][i];
    set_param_vec(trace1, num_param1, allparams1);
    // call solver iteration with y and J until convergence
    _omc_newton(wrapper_fvec_newton_adolc,data,reinterpret_cast<void*>(this));
    if (data->info <0){
    	throw FatalError(data->info, "Nonlinear solver failed!", __func__, __FILE__, __LINE__);
    }
    for(i=0;i<outsz[1];i++) {
    	y[1][i] = data->x[i];
    }
    double **J1 = myalloc2(outsz[1],outsz[1]);
    jacobian(trace1,outsz[1],outsz[1],y[1],J1);
    double **J3 = myalloc2(outsz[0],outsz[1]);
    set_param_vec(trace3, num_param1, allparams1);
    jacobian(trace3,outsz[0],outsz[1],y[1],J3);
    free(allparams1);

    // finally compute jacobian or resid,y1 wrt y2 after convergence in J
    size_t stats2[STAT_SIZE];
    tapestats(trace2, stats2);
    int num_depen2 = stats2[NUM_DEPENDENTS];// should be == outsz[0] + outsz[1]
    int num_indep2 = stats2[NUM_INDEPENDENTS]; // should be == insz[0]
    int num_param2 = stats2[NUM_PARAM]; // should be numouterparams + outsz[1]
    double* allparams2 = (double*) calloc(num_param2,sizeof(double));
    double* depen = (double*) calloc(num_depen2,sizeof(double));
    for (i = 0; i < numouterparams; i++)
        allparams2[i] = outerparams[i];
    for (i = 0; i < outsz[1]; i++)
        allparams2[numouterparams+i] = y[1][i];
    set_param_vec(trace2, num_param2, allparams2);
    // use fov_forward for computing j2 instead
    // allocate J2 as directional deriv of resid w.r.t. x
    double **J2 = myalloc2(num_depen2,ndir);
    ::fov_forward(trace2, num_depen2, num_indep2, ndir, x[0], Xp[0], depen, J2);
    for (i=0;i<outsz[0];i++) {
        	y[0][i] = depen[outsz[1]+i];
    }
    free(depen);
    // First outsz[1] rows of J2 contain dr/dx, next outsz[0] rows contain dy1/dx
    // solve yp[1] = - J1^{-1} * dr/dx using a linear solver

    int nrhs = ndir;
    int info;
    char trans = 'T';
    double* A = &J1[0][0];

    double** b = Yp[1];
    for(int i=0; i<outsz[1]; ++i){
        for(int j=0; j<ndir; ++j) {
            b[i][j] = -J2[i][j];
      //printf("b[%d] = %f\n", i, b[i])
        };
    }
    int *ipriv = (int*) calloc(outsz[1], sizeof(int));
    
    dgetrf_(&outsz[1], &outsz[1], A, &outsz[1], ipriv, &info);
    if (info == 0) {
        dgetrs_(&trans, &outsz[1], &nrhs, A, &outsz[1], ipriv, &b[0][0], &outsz[1], &info);
    }
    if (info < 0) {
        printf("fos_forward: Error solving linear system of equations. Argument %d illegal.\n", info);
        return -1;
    }
    free(ipriv);
    myfree2(J1);

    //free(A);
    // compute yp[0] = J3*yp[1] + dy1/dx <=> yp[0] = A*yp[1]+b

    A = &J3[0][0];
    double alpha = 1.0;
    double beta = 1.0;
    int incx = 1;
    b = Yp[0];
    for(int k=0; k<ndir; ++k) {
        for(int i=0; i<outsz[0]; ++i){
            b[i][k] = J2[outsz[1]+i][k];
            //printf("b[%d] = %f\n", i, b[i]);
        }
    }
    dgemm_(&trans,&trans,&ndir,&outsz[1],&outsz[1],&alpha,&Yp[1][0][0],&outsz[1],A,&outsz[1],&beta,&b[0][0],&ndir);
    myfree2(J3);
    myfree2(J2);
    free(outerparams);
    free(allparams2);
    return 0;
}


int NonLinearSolverEdf::fos_reverse(int iArrLen, int* iArr, int nout, int nin, int *outsz, double **up, int *insz, double **zp, double **x, double **y, void *ctx){
  return 0;
}
int NonLinearSolverEdf::fov_reverse(int iArrLen, int* iArr, int nout, int nin, int *outsz, int dir, double ***Up, int *insz, double ***Zp, double **x, double **y, void* ctx){
  return 0;
}


unsigned int alloc_adolc_lin_sol(int nnz, int nb, int nx) {
    int insz[2], outsz[1];
    insz[0] = nnz;
    insz[1] = nb;
    outsz[0] = nx;
    linSolEdfVec.emplace_back();
    linSolEdfVec.back().allocate_mem(2,1,insz,outsz);
    return linSolEdfVec.back().get_index();
}

unsigned int alloc_adolc_nonlin_sol(char* fbase,int nx, int ny1, int ny2,short* usetag) {
    int insz[1], outsz[2];
    insz[0] = nx;
    outsz[0] = ny1;
    outsz[1] = ny2;
    nonLinSolEdfVec.emplace_back(fbase,*usetag);
    NonLinearSolverEdf& edf = nonLinSolEdfVec.back();
    edf.allocate_mem(1,2,insz,outsz);
    *usetag = edf.get_next_tag();
    edf.numIterVar = outsz[1];
    allocateNewtonData(edf.numIterVar,reinterpret_cast<void**>(&edf.data));
    edf.tmp = (char*) malloc(edf.numIterVar * sizeof(double*));
    populate_dpp_with_contigdata(&edf.J,edf.tmp,edf.numIterVar,edf.numIterVar,edf.data->fjac);
    edf.I2 = myallocI2(edf.numIterVar);
    edf.data->trans = 'T';
    return edf.get_index();
}

double *adolc_nonlin_sol_get_values_buffer(int index) {
	NonLinearSolverEdf& edf = nonLinSolEdfVec[index];
	return edf.data->x;
}
