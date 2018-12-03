#ifndef SIMD_H
#define SIMD_H

#include <x86intrin.h>

typedef __m256d avx_vector;
typedef __m128d sse_vector;

#define avx_load(a,b) \
	a = _mm256_loadu_pd((double*)&b)

#define sse_load(a,b) \
	a = _mm_loadu_pd((double*)&b)

#define avx_store(a,b) \
	_mm256_storeu_pd((double*)&a,b)

#define sse_store(a,b) \
	_mm_storeu_pd((double*)&a,b)

#define avx_mul(a,b,c) \
	a = _mm256_mul_pd(b,c)

#define sse_mul(a,b,c) \
	a = _mm_mul_pd(b,c)

#define avx_mul_assign(a,b) \
	a = _mm256_mul_pd(a,b)

#define sse_mul_assign(a,b) \
	a = _mm_mul_pd(a,b)

#define avx_add(a,b,c) \
	a = _mm256_add_pd(b,c)

#define sse_add(a,b,c) \
	a = _mm_add_pd(b,c)

#define avx_add_assign(a,b) \
	a = _mm256_add_pd(a,b)

#define sse_add_assign(a,b) \
	a = _mm_add_pd(a,b)

#define avx_sub(a,b,c) \
	a = _mm256_sub_pd(b,c)

#define sse_sub(a,b,c) \
	a = _mm_sub_pd(b,c)

#define avx_sub_assign(a,b) \
	a = _mm256_sub_pd(a,b)

#define sse_sub_assign(a,b) \
	a = _mm_sub_pd(a,b)

#define avx_div(a,b,c) \
	a = _mm256_div_pd(b,c)

#define sse_div(a,b,c) \
	a = _mm_div_pd(b,c)

#define avx_set_dbl(a,b) \
	a = _mm256_set1_pd(b)

#define sse_set_dbl(a,b) \
	a = _mm_set1_pd(b)

#define avx_set_zero(a) \
	a = _mm256_setzero_pd()

#define sse_set_zero(a) \
	a = _mm_setzero_pd()

#define avx_sum(a,b) \
	do { \
		double *v = (double*)&b; \
		a = v[0] + v[1] + v[2] + v[3]; \
	} while(0)

#define sse_sum(a,b) \
	do { \
		double *v = (double*)&b; \
		a = v[0] + v[1]; \
	} while(0)

#endif
