#define _POSIX_C_SOURCE 200809L

#include "bml.h"
#include "../macros.h"
#include "../typed.h"

#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

// print content of file in stdout
void TYPED_FUNC(
    filecontent) (
    char *fname)
{
    FILE *fp = fopen(fname, "r");
    int ch;

    if (fp == NULL)
    {
        fprintf(stderr, "File %s doesn't exist\n", fname);
        return;
    }

    // Read contents of file
    while ((ch = fgetc(fp)) != EOF)
    {
        // print current character
        putc(ch, stdout);
    }

    fclose(fp);
}

int TYPED_FUNC(
    test_print) (
    const int N,
    const bml_matrix_type_t matrix_type,
    const bml_matrix_precision_t matrix_precision,
    const int M)
{
#ifndef BML_COMPLEX
    if (matrix_precision == single_complex
        || matrix_precision == double_complex)
    {
        printf("[FIXME] Skipping unsupported test\n");
        return 0;
    }
#endif

    //generate random matrix
    bml_matrix_t *A = NULL;
    A = bml_random_matrix(matrix_type, matrix_precision, N, M, sequential);
    REAL_T *A_dense = bml_export_to_dense(A, dense_row_major);
    bml_print_dense_matrix(N, matrix_precision, dense_row_major, A_dense, 0,
                           N, 0, N);
    /* Create unique filename (in case we run tests in parallel). */
    char *filename = strdup(tmpnam(NULL));
    fprintf(stdout, "Filename used for this test: %s\n", filename);
    int fd = open(filename, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
    if (fd < 0)
    {
        fprintf(stderr, "Failed to open %s\n", filename);
        return -1;
    }

    /* Flush stdout before redirecting to file. */
    fflush(stdout);
    int original_stdout = dup(fileno(stdout));
    if (dup2(fd, fileno(stdout)) < 0)
    {
        fprintf(stderr, "Failed to duplicate stdout\n");
        return -1;
    }

    /* Assumes matrix is at least 2x2 */
    const int up = 2;
    bml_print_bml_matrix(A, 0, up, 0, up);
    /* Flush stdout before switching back. */
    fflush(stdout);
    close(fd);
    /* Close file and re-instate stdout. */
    if (dup2(original_stdout, fileno(stdout)) < 0)
    {
        fprintf(stderr, "Failed to re-activate stdout\n");
        return -1;
    }

    printf("FILE content:\n");
    TYPED_FUNC(filecontent) (filename);
    //now read file just written
    REAL_T *data = calloc(up * up, sizeof(REAL_T));
    if (data == NULL)
    {
        fprintf(stderr, "calloc failed!\n");
        return -1;
    }

    FILE *fp2 = fopen(filename, "r");
    if (fp2 == NULL)
    {
        fprintf(stderr, "Failed to open %s\n", filename);
        return -1;
    }

    //FILE *fp3 = fopen("test_output.dat", "w");

    float realp;
    float imagp;
    char sign;
    for (int i = 0; i < up; i++)
    {
        for (int j = 0; j < up; j++)
        {
            switch (matrix_precision)
            {
                case single_real:
                case double_real:
                    fscanf(fp2, "%f", &realp);
                    printf("Read: %f\n", realp);
                    data[ROWMAJOR(i, j, up, up)] = realp;
                    break;
                case single_complex:
                case double_complex:
                    //read complex number in 3 parts, discarding 'i'
                    fscanf(fp2, "%f %c %fi", &realp, &sign, &imagp);
                    if (sign == '-')
                        imagp *= -1;
                    REAL_T tmp = realp + imagp * _Complex_I;
                    printf("realp %f\n", REAL_PART(tmp));
                    printf("imagp %f\n", IMAGINARY_PART(tmp));
                    data[ROWMAJOR(i, j, up, up)] = tmp;
                    break;
                default:
                    fprintf(stderr, "Unknown precision\n");
                    break;
            }
        }
    }
    if (fclose(fp2) == EOF)
    {
        fprintf(stderr, "ERROR closing file2\n");
        return -1;
    }

    //compare data just read with original data
    const double tol = 1.e-3;
    for (int i = 0; i < up; i++)
    {
        for (int j = 0; j < up; j++)
        {
            REAL_T val1 = data[ROWMAJOR(i, j, up, up)];
            REAL_T val2 = A_dense[ROWMAJOR(i, j, N, M)];
            double diff1 = REAL_PART(val2 - val1);
            double diff2 = IMAGINARY_PART(val2 - val1);
            double diff = sqrt(diff1 * diff1 + diff2 * diff2);
//            if (diff > tol)
            {
                printf("real parts= %f and %f\n", REAL_PART(val1),
                       REAL_PART(val2));
                if (matrix_precision == single_complex
                    || matrix_precision == double_complex)
                {
                    printf("imag part= %f and %f\n",
                           IMAGINARY_PART(val1), IMAGINARY_PART(val2));
                }
                printf("i=%d, j=%d, diff=%lf\n", i, j, diff);
                //fclose(fp3);
                if (diff > tol)
                {
                    fprintf(stderr, "test failed!!!\n");
                    return -1;
                }
            }
        }
    }
    //fprintf(fp3, "Test Successful!\n");
    //fclose(fp3);

    bml_deallocate(&A);
    free(data);
    if (remove(filename) != 0)
    {
        fprintf(stderr, "Failed removing file %s\n", filename);
        return -1;
    }
    free(filename);
    bml_free_memory(A_dense);
    return 0;
}
