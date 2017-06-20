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
public:
    NonLinearSolverEdf();
    virtual ~NonLinearSolverEdf() {}
    virtual int function(int iArrLen, int *iArr, int nin, int nout, int *insz, double **x, int *outsz, double **y, void* ctx);
    virtual int zos_forward(int iArrLen, int *iArr, int nin, int nout, int *insz, double **x, int *outsz, double **y, void* ctx);
    virtual int fos_forward(int iArrLen, int* iArr, int nin, int nout, int *insz, double **x, double **xp, int *outsz, double **y, double **yp, void *ctx);
    virtual int fov_forward(int iArrLen, int* iArr, int nin, int nout, int *insz, double **x, int ndir, double ***Xp, int *outsz, double **y, double ***Yp, void* ctx);
    virtual int fos_reverse(int iArrLen, int* iArr, int nout, int nin, int *outsz, double **up, int *insz, double **zp, double **x, double **y, void *ctx);
    virtual int fov_reverse(int iArrLen, int* iArr, int nout, int nin, int *outsz, int dir, double ***Up, int *insz, double ***Zp, double **x, double **y, void* ctx);
};

static std::vector<NonLinearSolverEdf> nonlinSolEdfVec;


NonLinearSolverEdf::NonLinearSolverEdf(){

}

int NonLinearSolverEdf::function(int iArrLen, int *iArr, int nin, int nout, int *insz, double **x, int *outsz, double **y, void* ctx) {


  return 0;
}

int NonLinearSolverEdf::zos_forward(int iArrLen, int *iArr, int nin, int nout, int *insz, double **x, int *outsz, double **y, void* ctx) {

  return this->function(iArrLen, iArr, nin, nout, insz, x, outsz, y, ctx);
}


int NonLinearSolverEdf::fos_forward(int iArrLen, int* iArr, int nin, int nout, int *insz, double **x, double **xp, int *outsz, double **y, double **yp, void *ctx) {


  return 0;
}



int NonLinearSolverEdf::fov_forward(int iArrLen, int* iArr, int nin, int nout, int *insz, double **x, int ndir, double ***Xp, int *outsz, double **y, double ***Yp, void* ctx) {

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

unsigned int alloc_adolc_nonlin_sol(int nnz, int nb, int nx) {
    int insz[2], outsz[1];
    insz[0] = nnz;
    insz[1] = nb;
    outsz[0] = nx;
    linSolEdfVec.emplace_back();
    linSolEdfVec.back().allocate_mem(2,1,insz,outsz);
    return linSolEdfVec.back().get_index();
}
