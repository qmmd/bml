#include "../../macros.h"
#include "../bml_allocate.h"
#include "../bml_logger.h"
#include "../bml_types.h"
#include "../lapack.h"
#include "bml_allocate_dense.h"
#include "bml_diagonalize_dense.h"
#include "bml_types_dense.h"
#include "../bml_utilities.h"

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
    const bml_matrix_dense_t * A,
    void *eigenvalues,
    bml_matrix_dense_t * eigenvectors)
{
    float *evecs = calloc(A->N * A->N, sizeof(float));
    float *evals = calloc(A->N, sizeof(float));
    int lwork = 3 * A->N;
    float *work = calloc(lwork, sizeof(float));
    int info;
    float *A_matrix;
    float *typed_eigenvalues = (float *) eigenvalues;

//    void mkl_thread_free_buffers(void);

    memcpy(evecs, A->matrix, A->N * A->N * sizeof(float));
#ifdef NOBLAS
    LOG_ERROR("No BLAS library");
#else
    C_SSYEV("V", "U", &A->N, evecs, &A->N, evals, work, &lwork, &info);
#endif
    // mkl_free_buffers();

    A_matrix = (float *) eigenvectors->matrix;
    for (int i = 0; i < A->N; i++)
    {
        typed_eigenvalues[i] = (float) evals[i];
        for (int j = 0; j < A->N; j++)
        {
            A_matrix[LINEARINDEX(i, j, A->N, A->N)] =
                evecs[COLMAJOR(i, j, A->N, A->N)];
        }
    }

    free(evals);
    free(evecs);
    free(work);
//    free(lwork);
//    mkl_thread_free_buffers();
}

void
bml_diagonalize_dense_double_real(
    const bml_matrix_dense_t * A,
    void *eigenvalues,
    bml_matrix_dense_t * eigenvectors)
{
    double *evecs = calloc(A->N * A->N, sizeof(double));
    double *evals = calloc(A->N, sizeof(double));
    int lwork = 3 * A->N;
    double *work = calloc(lwork, sizeof(double));
    int info;
    double *A_matrix;
    double *typed_eigenvalues = (double *) eigenvalues;

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
            A_matrix[LINEARINDEX(i, j, A->N, A->N)] =
                evecs[COLMAJOR(i, j, A->N, A->N)];
        }
    }

    free(evals);
    free(evecs);
    free(work);
//    free(lwork);
//    mkl_thread_free_buffers();
}

void
bml_diagonalize_dense_single_complex(
    const bml_matrix_dense_t * A,
    void *eigenvalues,
    bml_matrix_dense_t * eigenvectors)
{
    int *isuppz = calloc(2 * A->N, sizeof(int));
    int lwork = 2 * A->N;
    int liwork = 10 * A->N;
    int lrwork = 24 * A->N;
    int *iwork = calloc(liwork, sizeof(int));
    int M;
    int info;
    float *evals = calloc(A->N, sizeof(float));
    float *rwork = calloc(lrwork, sizeof(float));
    float abstol = 0;
    float complex *A_copy = calloc(A->N * A->N, sizeof(float complex));
    float complex *A_matrix;
    float complex *evecs = calloc(A->N * A->N, sizeof(float complex));
    float complex *work = calloc(lwork, sizeof(float complex));
    float complex *typed_eigenvalues = (float complex *) eigenvalues;

    memcpy(A_copy, A->matrix, A->N * A->N * sizeof(float complex));
#ifdef NOBLAS
    LOG_ERROR("No BLAS library");
#else
    C_CHEEVR("V", "A", "U", &A->N, A_copy, &A->N, NULL, NULL, NULL, NULL,
             &abstol, &M, evals, evecs, &A->N, isuppz, work, &lwork, rwork,
             &lrwork, iwork, &liwork, &info);
#endif
    A_matrix = (float complex *) eigenvectors->matrix;
    for (int i = 0; i < A->N; i++)
    {
        typed_eigenvalues[i] = (double) evals[i];
        for (int j = 0; j < A->N; j++)
        {
            A_matrix[LINEARINDEX(i, j, A->N, A->N)] =
                evecs[COLMAJOR(i, j, A->N, A->N)];
        }
    }

    free(A_copy);
    free(evecs);
    free(isuppz);
    free(work);
    free(rwork);
    free(iwork);
    free(evals);
}

void
bml_diagonalize_dense_double_complex(
    const bml_matrix_dense_t * A,
    void *eigenvalues,
    bml_matrix_dense_t * eigenvectors)
{
    int *isuppz = calloc(2 * A->N, sizeof(int));
    int liwork = 10 * A->N;
    int lrwork = 24 * A->N;
    int lwork = 2 * A->N;
    int *iwork = calloc(liwork, sizeof(int));
    int M;
    int info;
    double *evals = calloc(A->N, sizeof(double));
    double *rwork = calloc(lrwork, sizeof(double));
    double abstol = 0;
    double complex *A_copy = calloc(A->N * A->N, sizeof(double complex));
    double complex *A_matrix;
    double complex *evecs = calloc(A->N * A->N, sizeof(double complex));
    double complex *work = calloc(lwork, sizeof(double complex));
    double complex *typed_eigenvalues = (double complex *) eigenvalues;

    memcpy(A_copy, A->matrix, A->N * A->N * sizeof(double complex));
#ifdef NOBLAS
    LOG_ERROR("No BLAS library");
#else
    C_ZHEEVR("V", "A", "U", &A->N, A_copy, &A->N, NULL, NULL, NULL, NULL,
             &abstol, &M, evals, evecs, &A->N, isuppz, work, &lwork, rwork,
             &lrwork, iwork, &liwork, &info);
#endif
    A_matrix = (double complex *) eigenvectors->matrix;
    for (int i = 0; i < A->N; i++)
    {
        typed_eigenvalues[i] = (double) evals[i];
        for (int j = 0; j < A->N; j++)
        {
            A_matrix[LINEARINDEX(i, j, A->N, A->N)] =
                evecs[COLMAJOR(i, j, A->N, A->N)];
        }
    }

    free(A_copy);
    free(evecs);
    free(isuppz);
    free(work);
    free(rwork);
    free(iwork);
    free(evals);
}

void
bml_diagonalize_dense(
    const bml_matrix_dense_t * A,
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
