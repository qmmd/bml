#ifndef __LAPACK_H
#define __LAPACK_H

#include <complex.h>

void C_SSYEV(
    const char *JOBZ,
    const char *UPLO,
    const int *N,
    float *A,
    const int *LDA,
    float *W,
    float *WORK,
    const int *LWORK,
    int *INFO);

void C_DSYEV(
    const char *JOBZ,
    const char *UPLO,
    const int *N,
    double *A,
    const int *LDA,
    double *W,
    double *WORK,
    const int *LWORK,
    int *INFO);

void C_SSYEVR(
    const char *JOBZ,
    const char *RANGE,
    const char *UPLO,
    const int *N,
    float *A,
    const int *LDA,
    const float *VL,
    const float *VU,
    const int *IL,
    const int *IU,
    const float *ABSTOL,
    int *M,
    float *W,
    float *Z,
    const int *LDZ,
    int *ISUPPZ,
    float *WORK,
    const int *LWORK,
    int *IWORK,
    const int *LIWORK,
    int *INFO);

void C_DSYEVR(
    const char *JOBZ,
    const char *RANGE,
    const char *UPLO,
    const int *N,
    double *A,
    const int *LDA,
    const double *VL,
    const double *VU,
    const int *IL,
    const int *IU,
    const double *ABSTOL,
    int *M,
    double *W,
    double *Z,
    const int *LDZ,
    int *ISUPPZ,
    double *WORK,
    const int *LWORK,
    int *IWORK,
    const int *LIWORK,
    int *INFO);

void C_CHEEVR(
    const char *JOBZ,
    const char *RANGE,
    const char *UPLO,
    const int *N,
    float complex * A,
    const int *LDA,
    const float *VL,
    const float *VU,
    const int *IL,
    const int *IU,
    const float *ABSTOL,
    int *M,
    float *W,
    float complex * Z,
    const int *LDZ,
    int *ISUPPZ,
    float complex * WORK,
    const int *LWORK,
    float *RWORK,
    int *LRWORK,
    int *IWORK,
    const int *LIWORK,
    int *INFO);

void C_ZHEEVR(
    const char *JOBZ,
    const char *RANGE,
    const char *UPLO,
    const int *N,
    double complex * A,
    const int *LDA,
    const double *VL,
    const double *VU,
    const int *IL,
    const int *IU,
    const double *ABSTOL,
    int *M,
    double *W,
    double complex * Z,
    const int *LDZ,
    int *ISUPPZ,
    double complex * WORK,
    const int *LWORK,
    double *RWORK,
    int *LRWORK,
    int *IWORK,
    const int *LIWORK,
    int *INFO);

void C_SGETRF(
    const int *M,
    const int *N,
    float *A,
    const int *LDA,
    int *IPIV,
    int *INFO);

void C_DGETRF(
    const int *M,
    const int *N,
    double *A,
    const int *LDA,
    int *IPIV,
    int *INFO);

void C_CGETRF(
    const int *M,
    const int *N,
    float complex * A,
    const int *LDA,
    int *IPIV,
    int *INFO);

void C_ZGETRF(
    const int *M,
    const int *N,
    double complex * A,
    const int *LDA,
    int *IPIV,
    int *INFO);

void C_SGETRI(
    const int *N,
    float *A,
    const int *LDA,
    int *IPIV,
    float *WORK,
    const int *LWORK,
    int *INFO);

void C_DGETRI(
    const int *N,
    double *A,
    const int *LDA,
    int *IPIV,
    double *WORK,
    const int *LWORK,
    int *INFO);

void C_CGETRI(
    const int *N,
    float *A,
    const int *LDA,
    int *IPIV,
    float complex * WORK,
    const int *LWORK,
    int *INFO);

void C_ZGETRI(
    const int *N,
    double *A,
    const int *LDA,
    int *IPIV,
    double complex * WORK,
    const int *LWORK,
    int *INFO);

#endif
