#ifndef __BML_GEMM_H
#define __BML_GEMM_H

#include <complex.h>

void bml_gemm_single_real(
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
void bml_gemm_double_real(
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
void bml_gemm_single_complex(
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
void bml_gemm_double_complex(
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

void bml_xsmm_gemm_single_real(
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
void bml_xsmm_gemm_double_real(
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
void bml_xsmm_gemm_single_complex(
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
void bml_xsmm_gemm_double_complex(
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

#endif
