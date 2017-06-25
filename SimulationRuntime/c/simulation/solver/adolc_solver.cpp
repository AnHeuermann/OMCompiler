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
#include <adolc/drivers/drivers.h>
#include <adolc/tapedoc/asciitapes.h>

extern "C" int dgesv_(int *n, int *nrhs, double *a, int *lda,
                  int *ipiv, double *b, int *ldb, int *info);

extern "C" int dgemv_(char* trans, int* m, int *n, double *alpha, double *A,
                  int *lda, double *x, int *incx, double *beta, double *y, int *incy);

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

class NonLinearSolverEdf : public EDFobject_v2 {
protected:
    short trace1, trace2, nexttag;
public:
    NonLinearSolverEdf(const char* nlsfbase, short tagstart);
    virtual ~NonLinearSolverEdf() {}
    virtual int function(int iArrLen, int *iArr, int nin, int nout, int *insz, double **x, int *outsz, double **y, void* ctx);
    virtual int zos_forward(int iArrLen, int *iArr, int nin, int nout, int *insz, double **x, int *outsz, double **y, void* ctx);
    virtual int fos_forward(int iArrLen, int* iArr, int nin, int nout, int *insz, double **x, double **xp, int *outsz, double **y, double **yp, void *ctx);
    virtual int fov_forward(int iArrLen, int* iArr, int nin, int nout, int *insz, double **x, int ndir, double ***Xp, int *outsz, double **y, double ***Yp, void* ctx);
    virtual int fos_reverse(int iArrLen, int* iArr, int nout, int nin, int *outsz, double **up, int *insz, double **zp, double **x, double **y, void *ctx);
    virtual int fov_reverse(int iArrLen, int* iArr, int nout, int nin, int *outsz, int dir, double ***Up, int *insz, double ***Zp, double **x, double **y, void* ctx);
    short get_next_tag() const { return nexttag; }
};

static std::vector<NonLinearSolverEdf> nonLinSolEdfVec;


NonLinearSolverEdf::NonLinearSolverEdf(const char* nlsfbase, short tagstart) : EDFobject_v2() {
    char *nlsfilename1, *nlsfilename2;
    size_t slen = strlen(nlsfbase);
    nlsfilename1 = new char[slen + 10];
    nlsfilename2 = new char[slen + 10];
    sprintf(nlsfilename1,"%s_1_aat.txt",nlsfbase);
    sprintf(nlsfilename2,"%s_2_aat.txt",nlsfbase);
    trace1 = tagstart;
    tagstart = read_ascii_trace(nlsfilename1,tagstart);
    trace2 = tagstart;
    tagstart = read_ascii_trace(nlsfilename2,tagstart);
    nexttag = tagstart;
}

int NonLinearSolverEdf::function(int iArrLen, int *iArr, int nin, int nout, int *insz, double **x, int *outsz, double **y, void* ctx) {
    // assumption: nin == 1, nout == 1
    //             insz[0] == sizeof input vars x
    //             outsz[0] == sizeof output vars y
    int numouterparams, i;
    double *outerparams;
    numouterparams = alloc_copy_current_params(&outerparams);
    // solver takes input from x, jacobian from trace1 and gives output y
    size_t stats1[STAT_SIZE];
    tapestats(trace1, stats1);
    int num_resid1 = stats1[NUM_DEPENDENTS];
    int num_indep1 = stats1[NUM_INDEPENDENTS]; // should be == outsz[0]
    int num_param1 = stats1[NUM_PARAM]; // should be numouterparams + insz[0]
    double* allparams1 = (double*) calloc(num_param1,sizeof(double));
    for (i = 0; i < numouterparams; i++)
        allparams1[i] = outerparams[i];
    for (i = 0; i < insz[0]; i++)
        allparams1[numouterparams+i] = x[0][i];
    set_param_vec(trace1, num_param1, allparams1);
    double **J;
    // allocate J as jacobian of resid w.r.t. y (sparse or dense) depending on
    // solver
    // set initial values in y[0]
    // compute jacobian resid wrt y sparse or dense
    jacobian(trace1, num_resid1, num_indep1, y[0], J);
    // or using sparse_jac()
    // call solver iteration with y and j until convergence
    free(allparams1);
    free(outerparams);
    return 0;
}

int NonLinearSolverEdf::zos_forward(int iArrLen, int *iArr, int nin, int nout, int *insz, double **x, int *outsz, double **y, void* ctx) {

  return this->function(iArrLen, iArr, nin, nout, insz, x, outsz, y, ctx);
}


int NonLinearSolverEdf::fos_forward(int iArrLen, int* iArr, int nin, int nout, int *insz, double **x, double **xp, int *outsz, double **y, double **yp, void *ctx) {
    int numouterparams, i;
    double *outerparams;
    // do everything as in function then
    // finally compute jacobian or resid wrt y after convergence in J

    size_t stats2[STAT_SIZE];
    tapestats(trace2, stats2);
    int num_resid2 = stats2[NUM_DEPENDENTS];
    int num_indep2 = stats2[NUM_INDEPENDENTS]; // should be == insz[0]
    int num_param2 = stats2[NUM_PARAM]; // should be numouterparams + outsz[0]
    double* allparams2 = (double*) calloc(num_param2,sizeof(double));
    double* resid = (double*) calloc(num_resid2,sizeof(double));
    for (i = 0; i < numouterparams; i++)
        allparams2[i] = outerparams[i];
    for (i = 0; i < outsz[0]; i++)
        allparams2[numouterparams+i] = y[0][i];
    set_param_vec(trace2, num_param2, allparams2);
    // allocate J2 as directional deriv of resid w.r.t. x
    double *J2;
    ::fos_forward(trace2, num_resid2, num_indep2, 0, x[0], xp[0], resid, J2);
    // solve yp[0] = - J^{-1} * J2 using a linear solver

    free(outerparams);
    free(allparams2);
    free(resid);
    return 0;
}



int NonLinearSolverEdf::fov_forward(int iArrLen, int* iArr, int nin, int nout, int *insz, double **x, int ndir, double ***Xp, int *outsz, double **y, double ***Yp, void* ctx) {

    // almost same as fos_forward
    // use fov_forward for computing j2 instead

  return 0;
}


int NonLinearSolverEdf::fos_reverse(int iArrLen, int* iArr, int nout, int nin, int *outsz, double **up, int *insz, double **zp, double **x, double **y, void *ctx){
  return 0;
}
int NonLinearSolverEdf::fov_reverse(int iArrLen, int* iArr, int nout, int nin, int *outsz, int dir, double ***Up, int *insz, double ***Zp, double **x, double **y, void* ctx){
  return 0;
}


#include "simulation/solver/adolc_solver.h"

unsigned int alloc_adolc_lin_sol(int nnz, int nb, int nx) {
    int insz[2], outsz[1];
    insz[0] = nnz;
    insz[1] = nb;
    outsz[0] = nx;
    linSolEdfVec.emplace_back();
    linSolEdfVec.back().allocate_mem(2,1,insz,outsz);
    return linSolEdfVec.back().get_index();
}

unsigned int alloc_adolc_nonlin_sol(char* fbase,int nx, int ny,short* usetag) {
    int insz[1], outsz[1];
    insz[0] = nx;
    outsz[0] = ny;
    nonLinSolEdfVec.emplace_back(fbase,*usetag);
    nonLinSolEdfVec.back().allocate_mem(1,1,insz,outsz);
    *usetag = nonLinSolEdfVec.back().get_next_tag();
    return nonLinSolEdfVec.back().get_index();
}
