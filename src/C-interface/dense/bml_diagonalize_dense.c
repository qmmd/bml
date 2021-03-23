#include "../../macros.h"
#include "../bml_allocate.h"
#include "../bml_logger.h"
#include "../bml_types.h"
#include "bml_allocate_dense.h"
#include "bml_diagonalize_dense.h"
#include "bml_types_dense.h"
#include "../bml_utilities.h"
#include <float.h>

#ifdef BML_USE_MAGMA
#include "magma_v2.h"
#ifdef BML_USE_CUSOLVER
#include <cuda_runtime.h>
#include <cusolverDn.h>
#endif
#else
#include "../lapack.h"
#endif

#include <string.h>

/** \page diagonalize
 *
 * Note: We can't generify these functions easily since the API
 * differs between the real and complex types. rwork and lrwork are
 * only used in the complex cases. We opted instead to explicitly
 * implement the four versions.
 */

void
bml_diagonalize_dense_single_real(
    bml_matrix_dense_t * A,
    void *eigenvalues,
    bml_matrix_dense_t * eigenvectors)
{
    int info;
    float *A_matrix;
    float *typed_eigenvalues = (float *) eigenvalues;

//    void mkl_thread_free_buffers(void);

#ifdef BML_USE_MAGMA
    int nb = magma_get_ssytrd_nb(A->N);
    float *evecs;
    magma_int_t ret = magma_smalloc(&evecs, A->N * A->ld);
    assert(ret == MAGMA_SUCCESS);

    float *evals;
    evals = malloc(A->N * sizeof(float));
    float *work;
    int lwork = 2 * A->N + A->N * nb;
    int tmp = 1 + 6 * A->N + 2 * A->N * A->N;
    if (tmp > lwork)
        lwork = tmp;
    work = malloc(lwork * sizeof(float));
    int liwork = 3 + 5 * A->N;
    int *iwork;
    iwork = malloc(liwork * sizeof(int));
    float *wa;
    int ldwa = A->ld;
    wa = malloc(A->N * ldwa * sizeof(float));

    //copy matrix into evecs
    magmablas_slacpy(MagmaFull, A->N, A->N, A->matrix, A->ld, evecs, A->ld,
                     A->queue);

    magma_ssyevd_gpu(MagmaVec, MagmaUpper, A->N, evecs, A->ld, evals,
                     wa, ldwa, work, lwork, iwork, liwork, &info);
    if (info != 0)
        LOG_ERROR("ERROR in magma_ssyevd_gpu");

#else
    int lwork = 3 * A->N;
    float *evecs = calloc(A->N * A->N, sizeof(float));
    float *evals = calloc(A->N, sizeof(float));
    float *work = calloc(lwork, sizeof(float));
    memcpy(evecs, A->matrix, A->N * A->N * sizeof(float));

#ifdef NOBLAS
    LOG_ERROR("No BLAS library");
#else
    C_SSYEV("V", "U", &A->N, evecs, &A->N, evals, work, &lwork, &info);
#endif

#endif
    // mkl_free_buffers();

    A_matrix = (float *) eigenvectors->matrix;
#ifdef BML_USE_MAGMA
    magmablas_stranspose(A->N, A->N, evecs, A->ld,
                         A_matrix, eigenvectors->ld, A->queue);
    for (int i = 0; i < A->N; i++)
        typed_eigenvalues[i] = (float) evals[i];
#else
    for (int i = 0; i < A->N; i++)
    {
        typed_eigenvalues[i] = (float) evals[i];
        for (int j = 0; j < A->N; j++)
        {
            A_matrix[ROWMAJOR(i, j, A->N, A->N)] =
                evecs[COLMAJOR(i, j, A->N, A->N)];
        }
    }
#endif

#ifdef BML_USE_MAGMA
    magma_free(evecs);
    magma_free(wa);
#else
    free(evecs);
#endif
    free(evals);
    free(work);

//    free(lwork);
//    mkl_thread_free_buffers();
}

void
bml_diagonalize_dense_double_real(
    bml_matrix_dense_t * A,
    void *eigenvalues,
    bml_matrix_dense_t * eigenvectors)
{
    int info;
    double *A_matrix;
    double *typed_eigenvalues = (double *) eigenvalues;

#ifdef BML_USE_MAGMA
    //copy matrix into evecs
    double *evecs;
    magma_int_t ret = magma_dmalloc(&evecs, A->N * A->ld);
    assert(ret == MAGMA_SUCCESS);

    magmablas_dlacpy(MagmaFull, A->N, A->N, A->matrix, A->ld, evecs, A->ld,
                     A->queue);

    magma_queue_sync(A->queue);

#ifdef BML_USE_CUSOLVER
    // create cusolver/cublas handle
    cusolverDnHandle_t cusolverH = NULL;
    cusolverStatus_t cusolver_status = cusolverDnCreate(&cusolverH);
    assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);

    // allocate memory for eigenvalues
    double *d_W = NULL;
    cudaError_t cudaStat = cudaMalloc((void **) &d_W, sizeof(double) * A->N);
    assert(cudaSuccess == cudaStat);

    // compute eigenvalues and eigenvectors
    cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_VECTOR;
    cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;

    // query working space of syevd
    int lwork = 0;
    cusolver_status =
        cusolverDnDsyevd_bufferSize(cusolverH, jobz, uplo, A->N, evecs, A->ld,
                                    d_W, &lwork);
    assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

    double *d_work = NULL;
    cudaStat = cudaMalloc((void **) &d_work, sizeof(double) * lwork);
    assert(cudaSuccess == cudaStat);

    // solve
    int *devInfo = NULL;
    cudaStat = cudaMalloc((void **) &devInfo, sizeof(int));
    assert(cudaSuccess == cudaStat);

    cusolver_status =
        cusolverDnDsyevd(cusolverH, jobz, uplo, A->N, evecs, A->ld, d_W,
                         d_work, lwork, devInfo);
    cudaStat = cudaDeviceSynchronize();
    assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
    assert(cudaSuccess == cudaStat);

    // copy eigenvalues to CPU
    cudaStat =
        cudaMemcpy(typed_eigenvalues, d_W, sizeof(double) * A->N,
                   cudaMemcpyDeviceToHost);
    assert(cudaSuccess == cudaStat);

    // free resources
    cudaFree(d_W);
    cudaFree(devInfo);
    cudaFree(d_work);

    if (cusolverH)
        cusolverDnDestroy(cusolverH);
#else // MAGMA
    int nb = magma_get_ssytrd_nb(A->N);

    double *evals = malloc(A->N * sizeof(double));
    int lwork = 2 * A->N + A->N * nb;
    int tmp = 1 + 6 * A->N + 2 * A->N * A->N;
    if (tmp > lwork)
        lwork = tmp;
    double *work = malloc(lwork * sizeof(double));
    int liwork = 3 + 5 * A->N;
    int *iwork = malloc(liwork * sizeof(int));
    int ldwa = A->N;
    double *wa = malloc(A->N * ldwa * sizeof(double));

    magma_dsyevd_gpu(MagmaVec, MagmaUpper, A->N, evecs, A->ld, evals,
                     wa, ldwa, work, lwork, iwork, liwork, &info);
    if (info != 0)
        LOG_ERROR("ERROR in magma_dsyevd_gpu");

    free(wa);
    free(work);

    for (int i = 0; i < A->N; i++)
        typed_eigenvalues[i] = (double) evals[i];
    free(evals);

    //verify norm eigenvectors
    //for(int i=0;i<A->N;i++)
    //{
    //    double norm = magma_dnrm2(A->N, evecs+A->ld*i, 1, A->queue);
    //    printf("norm = %le\n", norm);
    //}
    magma_queue_sync(A->queue);
#endif

    // transpose eigenvactors matrix on GPU
    A_matrix = (double *) eigenvectors->matrix;
    magmablas_dtranspose(A->N, A->N, evecs, A->ld,
                         A_matrix, eigenvectors->ld, A->queue);
    magma_queue_sync(eigenvectors->queue);

    //verify norm eigenvectors transposed
    //for(int i=0;i<A->N;i++)
    //{
    //    double norm = magma_dnrm2(A->N, evecs+i, A->ld, A->queue);
    //    printf("norm transposed vector = %le\n", norm);
    //}
    magma_free(evecs);

#else // CPU code
    int lwork = 3 * A->N;
    double *evecs = calloc(A->N * A->N, sizeof(double));
    double *evals = calloc(A->N, sizeof(double));
    double *work = calloc(lwork, sizeof(double));
    memcpy(evecs, A->matrix, A->N * A->N * sizeof(double));

#ifdef NOBLAS
    LOG_ERROR("No BLAS library");
#else
    C_DSYEV("V", "U", &A->N, evecs, &A->N, evals, work, &lwork, &info);
#endif

    A_matrix = (double *) eigenvectors->matrix;
    for (int i = 0; i < A->N; i++)
    {
        typed_eigenvalues[i] = (double) evals[i];
        for (int j = 0; j < A->N; j++)
        {
            A_matrix[ROWMAJOR(i, j, A->N, A->N)] =
                evecs[COLMAJOR(i, j, A->N, A->N)];
        }
    }
    free(evecs);
    free(evals);
    free(work);
#endif
}

void
bml_diagonalize_dense_single_complex(
    bml_matrix_dense_t * A,
    void *eigenvalues,
    bml_matrix_dense_t * eigenvectors)
{
    int info;
    float _Complex *typed_eigenvalues = (float _Complex *) eigenvalues;

#ifdef BML_USE_MAGMA
    int nb = magma_get_ssytrd_nb(A->N);
    magmaFloatComplex *evecs;
    magma_int_t ret = magma_cmalloc(&evecs, A->N * A->ld);
    assert(ret == MAGMA_SUCCESS);

    float *evals = malloc(A->N * sizeof(float));
    int lwork = 2 * A->N + A->N * nb;
    int tmp = 1 + 6 * A->N + 2 * A->N * A->N;
    if (tmp > lwork)
        lwork = tmp;
    magmaFloatComplex *work = malloc(lwork * sizeof(magmaFloatComplex));
    int liwork = 3 + 5 * A->N;
    int *iwork = malloc(liwork * sizeof(int));
    int ldwa = A->ld;
    magmaFloatComplex *wa = malloc(A->N * ldwa * sizeof(magmaFloatComplex));
    int lrwork = 1 + 5 * A->N + 2 * A->N * A->N;
    float *rwork = malloc(lrwork * sizeof(float));

    //copy matrix into evecs
    magmablas_clacpy(MagmaFull, A->N, A->N, A->matrix, A->ld, evecs, A->ld,
                     A->queue);

    magma_queue_sync(A->queue);

    magma_cheevd_gpu(MagmaVec, MagmaUpper, A->N, evecs, A->ld, evals,
                     wa, ldwa,
                     work, lwork, rwork, lrwork, iwork, liwork, &info);
    if (info != 0)
        LOG_ERROR("ERROR in magma_cheevd_gpu");

    magmaFloatComplex *A_matrix = (magmaFloatComplex *) eigenvectors->matrix;
    magmablas_ctranspose(A->N, A->N, evecs, A->ld,
                         A_matrix, eigenvectors->ld, A->queue);
    for (int i = 0; i < A->N; i++)
        typed_eigenvalues[i] = (double) evals[i];

    magma_free(evecs);
    free(wa);
#else
    int *isuppz = calloc(2 * A->N, sizeof(int));
    int lwork = 2 * A->N;
    int liwork = 10 * A->N;
    int lrwork = 24 * A->N;
    int *iwork = calloc(liwork, sizeof(int));
    int M;
    float *evals = calloc(A->N, sizeof(float));
    float *rwork = calloc(lrwork, sizeof(float));
    float abstol = 0;
    float complex *A_copy = calloc(A->N * A->N, sizeof(float complex));
    float complex *A_matrix;
    float complex *evecs = calloc(A->N * A->N, sizeof(float complex));
    float complex *work = calloc(lwork, sizeof(float complex));

    memcpy(A_copy, A->matrix, A->N * A->N * sizeof(float complex));
#ifdef NOBLAS
    LOG_ERROR("No BLAS library");
#else
    {
        float vl = -FLT_MAX;
        float vu = FLT_MAX;
        int il = 1;
        int iu = A->N;
        C_CHEEVR("V", "A", "U", &A->N, A_copy, &A->N, &vl, &vu, &il, &iu,
	          &abstol, &M, evals, evecs, &A->N, isuppz, work, &lwork,
		  rwork, &lrwork, iwork, &liwork, &info);
    }
#endif
    A_matrix = (float complex *) eigenvectors->matrix;
    for (int i = 0; i < A->N; i++)
    {
        typed_eigenvalues[i] = (double) evals[i];
        for (int j = 0; j < A->N; j++)
        {
            A_matrix[ROWMAJOR(i, j, A->N, A->N)] =
                evecs[COLMAJOR(i, j, A->N, A->N)];
        }
    }

    free(A_copy);
    free(evecs);
    free(isuppz);
#endif
    free(work);
    free(rwork);
    free(iwork);
    free(evals);
}

void
bml_diagonalize_dense_double_complex(
    bml_matrix_dense_t * A,
    void *eigenvalues,
    bml_matrix_dense_t * eigenvectors)
{
    int info;
    double _Complex *typed_eigenvalues = (double _Complex *) eigenvalues;

#ifdef BML_USE_MAGMA
    int nb = magma_get_ssytrd_nb(A->N);
    magmaDoubleComplex *evecs;
    magma_int_t ret = magma_zmalloc(&evecs, A->N * A->ld);
    assert(ret == MAGMA_SUCCESS);

    double *evals = malloc(A->N * sizeof(double));
    int lwork = 2 * A->N + A->N * nb;
    int tmp = 1 + 6 * A->N + 2 * A->N * A->N;
    if (tmp > lwork)
        lwork = tmp;
    magmaDoubleComplex *work = malloc(lwork * sizeof(magmaDoubleComplex));
    int liwork = 3 + 5 * A->N;
    int *iwork = malloc(liwork * sizeof(int));
    int ldwa = A->ld;
    magmaDoubleComplex *wa = malloc(A->N * ldwa * sizeof(magmaDoubleComplex));
    int lrwork = 1 + 5 * A->N + 2 * A->N * A->N;
    double *rwork = malloc(lrwork * sizeof(double));

    //copy matrix into evecs
    magmablas_zlacpy(MagmaFull, A->N, A->N, A->matrix, A->ld, evecs, A->ld,
                     A->queue);

    magma_zheevd_gpu(MagmaVec, MagmaUpper, A->N, evecs, A->ld, evals,
                     wa, ldwa,
                     work, lwork, rwork, lrwork, iwork, liwork, &info);
    if (info != 0)
        LOG_ERROR("ERROR in magma_csyevd_gpu");

    magmaDoubleComplex *A_matrix =
        (magmaDoubleComplex *) eigenvectors->matrix;
    magmablas_ztranspose(A->N, A->N, evecs, A->ld, A_matrix, eigenvectors->ld,
                         A->queue);
    for (int i = 0; i < A->N; i++)
        typed_eigenvalues[i] = (double) evals[i];

    magma_free(evecs);
    free(wa);
#else
    int *isuppz = calloc(2 * A->N, sizeof(int));
    int liwork = 10 * A->N;
    int lrwork = 24 * A->N;
    int lwork = 2 * A->N;
    int *iwork = calloc(liwork, sizeof(int));
    int M;
    double *evals = calloc(A->N, sizeof(double));
    double *rwork = calloc(lrwork, sizeof(double));
    double abstol = 0;
    double complex *A_copy = calloc(A->N * A->N, sizeof(double complex));
    double complex *A_matrix;
    double complex *evecs = calloc(A->N * A->N, sizeof(double complex));
    double complex *work = calloc(lwork, sizeof(double complex));

    memcpy(A_copy, A->matrix, A->N * A->N * sizeof(double complex));
#ifdef NOBLAS
    LOG_ERROR("No BLAS library");
#else
    {
        double vl = -DBL_MAX;
        double vu = DBL_MAX;
        int il = 1;
        int iu = A->N;

        C_ZHEEVR("V", "A", "U", &A->N, A_copy, &A->N, &vl, &vu, &il, &iu,
	          &abstol, &M, evals, evecs, &A->N, isuppz, work, &lwork,
		  rwork, &lrwork, iwork, &liwork, &info);
    }
#endif
    A_matrix = (double complex *) eigenvectors->matrix;
    for (int i = 0; i < A->N; i++)
    {
        typed_eigenvalues[i] = (double) evals[i];
        for (int j = 0; j < A->N; j++)
        {
            A_matrix[ROWMAJOR(i, j, A->N, A->N)] =
                evecs[COLMAJOR(i, j, A->N, A->N)];
        }
    }

    free(A_copy);
    free(evecs);
    free(isuppz);
#endif
    free(work);
    free(rwork);
    free(iwork);
    free(evals);
}

void
bml_diagonalize_dense(
    bml_matrix_dense_t * A,
    void *eigenvalues,
    bml_matrix_t * eigenvectors)
{
    switch (A->matrix_precision)
    {
        case single_real:
            bml_diagonalize_dense_single_real(A, eigenvalues, eigenvectors);
            break;
        case double_real:
            bml_diagonalize_dense_double_real(A, eigenvalues, eigenvectors);
            break;
        case single_complex:
            bml_diagonalize_dense_single_complex(A, eigenvalues,
                                                 eigenvectors);
            break;
        case double_complex:
            bml_diagonalize_dense_double_complex(A, eigenvalues,
                                                 eigenvectors);
            break;
        default:
            LOG_ERROR("unknown precision\n");
            break;
    }
}
