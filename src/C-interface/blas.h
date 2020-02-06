#ifndef __BLAS_H
#define __BLAS_H

#include <complex.h>

#include "../typed.h"

void C_SSCAL(
    const int *n,
    const float *a,
    float *x,
    const int *incx);
void C_DSCAL(
    const int *n,
    const double *a,
    double *x,
    const int *incx);
void C_CSCAL(
    const int *n,
    const float complex *a,
    float complex *x,
    const int *incx);
void C_ZSCAL(
    const int *n,
    const double complex *a,
    double complex *x,
    const int *incx);
void C_SGEMM(
    const char *transa,
    const char *transb,
    const int *m,
    const int *n,
    const int *k,
    const float *alpha,
    const float *a,
    const int *lda,
    const float *b,
    const int *ldb,
    const float *beta,
    float *c,
    const int *ldc);
void C_DGEMM(
    const char *transa,
    const char *transb,
    const int *m,
    const int *n,
    const int *k,
    const double *alpha,
    const double *a,
    const int *lda,
    const double *b,
    const int *ldb,
    const double *beta,
    double *c,
    const int *ldc);
void C_CGEMM(
    const char *transa,
    const char *transb,
    const int *m,
    const int *n,
    const int *k,
    const float complex *alpha,
    const float complex *a,
    const int *lda,
    const float complex *b,
    const int *ldb,
    const float complex *beta,
    float complex *c,
    const int *ldc);
void C_ZGEMM(
    const char *transa,
    const char *transb,
    const int *m,
    const int *n,
    const int *k,
    const double complex *alpha,
    const double complex *a,
    const int *lda,
    const double complex *b,
    const int *ldb,
    const double complex *beta,
    double complex *c,
    const int *ldc);
void C_SAXPY(
    const int *n,
    const float *alpha,
    const float *x,
    const int *incx,
    float *y,
    const int *incy);
void C_DAXPY(
    const int *n,
    const double *alpha,
    const double *x,
    const int *incx,
    double *y,
    const int *incy);
void C_CAXPY(
    const int *n,
    const float complex *alpha,
    const float complex *x,
    const int *incx,
    float complex *y,
    const int *incy);
void C_ZAXPY(
    const int *n,
    const double complex *alpha,
    const double complex *x,
    const int *incx,
    double complex *y,
    const int *incy);

void XSMM(
    C_SGEMM) (
    const char *transa,
    const char *transb,
    const int *m,
    const int *n,
    const int *k,
    const float *alpha,
    const float *a,
    const int *lda,
    const float *b,
    const int *ldb,
    const float *beta,
    float *c,
    const int *ldc);
void XSMM(
    C_DGEMM) (
    const char *transa,
    const char *transb,
    const int *m,
    const int *n,
    const int *k,
    const double *alpha,
    const double *a,
    const int *lda,
    const double *b,
    const int *ldb,
    const double *beta,
    double *c,
    const int *ldc);
#endif
