#ifndef __TYPED_H
#define __TYPED_H

/* Fix the generated function name. */
#if defined(SINGLE_REAL)
#define FUNC_SUFFIX single_real
#elif defined(DOUBLE_REAL)
#define FUNC_SUFFIX double_real
#elif defined(SINGLE_COMPLEX)
#define FUNC_SUFFIX single_complex
#elif defined(DOUBLE_COMPLEX)
#define FUNC_SUFFIX double_complex
#else
#error Unknown precision type
#endif

/* Define numeric types. */
#if defined(SINGLE_REAL) || (defined(SINGLE_COMPLEX) && ! defined(BML_COMPLEX))
#define REAL_T float
#define MAGMA_T float
#define MATRIX_PRECISION single_real
#define BLAS_PREFIX S
#define MAGMA_PREFIX s
#define REAL_PART(x) (x)
#define IMAGINARY_PART(x) (0)
#define COMPLEX_CONJUGATE(x) (x)
#define ABS(x) (fabsf(x))
#define is_above_threshold(x, t) (fabsf(x) > (float) (t))
#elif defined(DOUBLE_REAL) || (defined(DOUBLE_COMPLEX) && ! defined(BML_COMPLEX))
#define REAL_T double
#define MAGMA_T double
#define MATRIX_PRECISION double_real
#define BLAS_PREFIX D
#define MAGMA_PREFIX d
#define REAL_PART(x) (x)
#define IMAGINARY_PART(x) (0)
#define COMPLEX_CONJUGATE(x) (x)
#define ABS(x) (fabs(x))
#define is_above_threshold(x, t) (fabs(x) > (double) (t))
#elif defined(SINGLE_COMPLEX)
#define REAL_T float _Complex
#define MAGMA_T magmaFloatComplex
#define MATRIX_PRECISION single_complex
#define BLAS_PREFIX C
#define MAGMA_PREFIX c
#define REAL_PART(x) (crealf(x))
#define IMAGINARY_PART(x) (cimagf(x))
#define COMPLEX_CONJUGATE(x) (conjf(x))
#define ABS(x) (cabsf(x))
#define is_above_threshold(x, t) ((creal(x)*creal(x)+cimag(x)+cimag(x)) > (float)(t*t))
#elif defined(DOUBLE_COMPLEX)
#define REAL_T double _Complex
#define MAGMA_T magmaDoubleComplex
#define MATRIX_PRECISION double_complex
#define BLAS_PREFIX Z
#define MAGMA_PREFIX z
#define REAL_PART(x) (creal(x))
#define IMAGINARY_PART(x) (cimag(x))
#define COMPLEX_CONJUGATE(x) (conj(x))
#define ABS(x) (cabs(x))
#define is_above_threshold(x, t) ((creal(x)*creal(x)+cimag(x)+cimag(x)) > (double)(t*t))
#else
#error Unknown precision type
#endif

#define CONCAT2_(a, b) a ## _ ## b
#define CONCAT_(a, b) CONCAT2_(a, b)

#define CONCAT2(a, b) a ## b
#define CONCAT(a, b) CONCAT2(a, b)

#define TYPED_FUNC(a) CONCAT_(a, FUNC_SUFFIX)
#define C_BLAS(a) CONCAT_(C, CONCAT(BLAS_PREFIX , a))
#define MAGMACOMPLEX(a) CONCAT_(MAGMA, CONCAT_(BLAS_PREFIX, a))
#define MAGMA(a) CONCAT_(magma, CONCAT(MAGMA_PREFIX , a))
#define MAGMAGPU(a) CONCAT_(magma, CONCAT(MAGMA_PREFIX , CONCAT_(a, gpu)))
#define MAGMABLAS(a) CONCAT_(magmablas, CONCAT(MAGMA_PREFIX , a))
#endif
