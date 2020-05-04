#define _POSIX_C_SOURCE 200809L

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

#include "bml_test.h"

const int NUM_TESTS = 28;

typedef struct
{
    char *test_name;
    char *test_description;
    const test_function_t tester;
} tests_t;

const char *test_name[] = {
    "add",
    "adjacency",
    "adjungate_triangle",
    "allocate",
    "bml_gemm",
    "import_export",
    "convert",
    "copy",
    "diagonalize",
    "get_element",
    "get_set_diagonal",
    "get_sparsity",
    "inverse",
    "io_matrix",
    "multiply",
    "multiply_banded",
    "multiply_x2",
    "norm",
    "normalize",
    "print",
    "scale",
    "set_element",
    "set_row",
    "submatrix",
    "threshold",
    "trace",
    "trace_mult",
    "transpose"
};

const char *test_description[] = {
    "Add two bml matrices",
    "Adjacency CSR arrays for metis",
    "Adjungate triangle (conjugate transpose) of bml matrices",
    "Allocate bml matrices",
    "Internal GEMM implmentation",
    "Convert by import/export of bml matrices",
    "Convert bml matrix",
    "Copy bml matrices",
    "Diagonalize a bml matrix",
    "Get an element from a bml matrix",
    "Set the diagonal elements of bml matrices",
    "Get the sparsity",
    "Matrix inverse",
    "Read and write an mtx matrix",
    "Multiply two bml matrices",
    "Multiply two banded bml matrices",
    "Multiply two identical matrices",
    "Norm of bml matrix",
    "Normalize bml matrices",
    "Print bml matrix to stdout",
    "Scale bml matrices",
    "Set a single element of a bml matrix",
    "Set the elements of a row in a bml matrix",
    "Submatrix bml matrices",
    "Threshold bml matrices",
    "Trace of bml matrices",
    "Trace from multiplication of two bml matrices",
    "Transpose of bml matrices"
};

const test_function_t testers[] = {
    test_add,
    test_adjacency,
    test_adjungate_triangle,
    test_allocate,
    test_bml_gemm,
    test_import_export,
    test_convert,
    test_copy,
    test_diagonalize,
    test_get_element,
    test_get_set_diagonal,
    test_get_sparsity,
    test_inverse,
    test_io_matrix,
    test_multiply,
    test_multiply_banded,
    test_multiply_x2,
    test_norm,
    test_normalize,
    test_print,
    test_scale,
    test_set_element,
    test_set_row,
    test_submatrix,
    test_threshold,
    test_trace,
    test_trace_mult,
    test_transpose
};

void
print_usage(
    void)
{
    printf("Usage:\n");
    printf("\n");
    printf("-h | --help           This help\n");
    printf("-n | --testname TEST  Run test TEST\n");
    printf("-t | --type T         Choose the matrix type\n");
    printf("-p | --precision P    Choose matrix precision\n");
    printf("-l | --list           List all available tests\n");
    printf("-N | --N N            Test N x N matrices\n");
    printf("-M | --M M            Pass an extra parameter M to the test\n");
    printf("\n");
    printf("Recognized types:\n");
    printf("\n");
    printf("  dense\n");
    printf("  ellpack\n");
    printf("  ellsort\n");
    printf("  ellblock\n");
    printf("  csr\n");
    printf("\n");
    printf("Recognized precisions:\n");
    printf("\n");
    printf("  single_real\n");
    printf("  double_real,\n");
    printf("  single_complex\n");
    printf("  double_complex\n");
    printf("\n");

    int max_width = 0;
    for (int i = 0; i < NUM_TESTS; i++)
    {
        if (strlen(test_name[i]) > max_width)
        {
            max_width = strlen(test_name[i]);
        }
    }
    char desc_format[100];
    snprintf(desc_format, 100, "%%%ds    %%s\n", max_width);

    printf("Available tests:\n");
    printf("\n");
    for (int i = 0; i < NUM_TESTS; i++)
    {
        printf(desc_format, test_name[i], test_description[i]);
    }
}

int
main(
    int argc,
    char **argv)
{
    bml_init(&argc, &argv);

    int N = 13;
    int M = -1;
    char *test = NULL;
    int test_index = -1;
    int test_result;
    bml_matrix_type_t matrix_type = dense;
    bml_matrix_precision_t precision = single_real;

    const char *short_options = "hn:t:p:lN:M:";
    const struct option long_options[] = {
        {"help", no_argument, NULL, 'h'},
        {"testname", required_argument, NULL, 'n'},
        {"type", required_argument, NULL, 't'},
        {"precision", required_argument, NULL, 'p'},
        {"list", no_argument, NULL, 'l'},
        {"N", required_argument, NULL, 'N'},
        {"M", required_argument, NULL, 'M'},
        {NULL, 0, NULL, 0}
    };
    int c;

    while ((c =
            getopt_long(argc, argv, short_options, long_options, NULL)) != -1)
    {
        switch (c)
        {
            case 'h':
                print_usage();
                return 0;
                break;
            case 'n':
            {
                if (test)
                {
                    free(test);
                }
                test = strdup(optarg);
                for (test_index = 0; test_index < NUM_TESTS; test_index++)
                {
                    if (strcasecmp(test, test_name[test_index]) == 0)
                    {
                        break;
                    }
                }
                if (test_index == NUM_TESTS)
                {
                    fprintf(stderr, "unknown test %s\n", test);
                    return 1;
                }
                break;
            }
            case 't':
                if (strcasecmp(optarg, "dense") == 0)
                {
                    matrix_type = dense;
                }
                else if (strcasecmp(optarg, "ellpack") == 0)
                {
                    matrix_type = ellpack;
                }
                else if (strcasecmp(optarg, "ellsort") == 0)
                {
                    matrix_type = ellsort;
                }
                else if (strcasecmp(optarg, "ellblock") == 0)
                {
                    matrix_type = ellblock;
                }
                else if (strcasecmp(optarg, "csr") == 0)
                {
                    matrix_type = csr;
                }
                else
                {
                    fprintf(stderr, "unknown matrix type %s\n", optarg);
                }
                break;
            case 'p':
                if (strcasecmp(optarg, "single_real") == 0)
                {
                    precision = single_real;
                }
                else if (strcasecmp(optarg, "double_real") == 0)
                {
                    precision = double_real;
                }
                else if (strcasecmp(optarg, "single_complex") == 0)
                {
                    precision = single_complex;
                }
                else if (strcasecmp(optarg, "double_complex") == 0)
                {
                    precision = double_complex;
                }
                else
                {
                    fprintf(stderr, "unknow matrix precision %s\n", optarg);
                }
                break;
            case 'l':
                for (int i = 0; i < NUM_TESTS; i++)
                {
                    printf("%s", test_name[i]);
                    if (i < NUM_TESTS - 1)
                    {
                        printf(", ");
                    }
                }
                printf("\n");
                break;
            case 'N':
                N = strtol(optarg, NULL, 10);
                break;
            case 'M':
                M = strtol(optarg, NULL, 10);
                break;
            default:
                fprintf(stderr, "unknown option\n");
                return 1;
                break;
        }
    }

    if (!test)
    {
        fprintf(stderr, "missing test\n");
        return 1;
    }

    if (M < 0)
    {
        M = N;
    }

    fprintf(stderr, "%s\n", test);
    fprintf(stderr, "N = %d\n", N);
    free(test);

    test_result = testers[test_index] (N, matrix_type, precision, M);

    bml_shutdown();

    return test_result;
}
