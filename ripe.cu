/*
* haval-256 kernel implementation.
*
* ==========================(LICENSE BEGIN)============================
*
* Copyright (c) 2014 djm34
*               2016 tpruvot
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the
* "Software"), to deal in the Software without restriction, including
* without limitation the rights to use, copy, modify, merge, publish,
* distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to
* the following conditions:
*
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*
* ===========================(LICENSE END)=============================
*/

#ifdef __INTELLISENSE__
#define atomicExch(p,y) y
#define __byte_perm(x,y,z) x
#endif

#include <stdio.h>
#include <memory.h>

#include "cuda_helper.h"

static __constant__ const uint32_t c_IV[5] = { 0x67452301u, 0xEFCDAB89u, 0x98BADCFEu, 0x10325476u, 0xC3D2E1F0u };
static __constant__ const uint32_t c_K1[4] = { 0, 0x5A827999, 0x6ED9EBA1, 0x8F1BBCDC};
static __constant__ const uint32_t c_K2[4] = { 0x50A28BE6, 0x5C4DD124, 0x6D703EF3, 0 };
//#define sK11   SPH_C32(0x00000000)
//#define sK12   SPH_C32(0x5A827999)
//#define sK13   SPH_C32(0x6ED9EBA1)
//#define sK14   SPH_C32(0x8F1BBCDC)

//#define sK21   SPH_C32(0x50A28BE6)
//#define sK22   SPH_C32(0x5C4DD124)
//#define sK23   SPH_C32(0x6D703EF3)
//#define sK24   SPH_C32(0x00000000)

/*
* Round functions for RIPEMD-160.
*/
//#define F1(x, y, z)   xor3x(x, y, z)

__device__ __forceinline__
static uint32_t ROTATE(const uint32_t x, const uint32_t r) {
	if (r == 8)
		return __byte_perm(x, 0, 0x2103);
	else
		return ROTL32(x, r);
}

__device__ __forceinline__
uint32_t F1(const uint32_t a, const uint32_t b, const uint32_t c) {
	uint32_t result;
#if __CUDA_ARCH__ >= 500 && CUDA_VERSION >= 7050
	asm volatile ("lop3.b32 %0, %1, %2, %3, 0x96;" : "=r"(result) : "r"(a), "r"(b), "r"(c));
#else
	result = a^b^c;
#endif
	return result;
}
//#define F2(x, y, z)   ((x & (y ^ z)) ^ z)
__device__ __forceinline__
uint32_t F2(const uint32_t a, const uint32_t b, const uint32_t c) {
	uint32_t result;
#if __CUDA_ARCH__ >= 500 && CUDA_VERSION >= 7050
	asm volatile ("lop3.b32 %0, %1, %2, %3, 0xCA;" : "=r"(result) : "r"(a), "r"(b), "r"(c)); //0xCA=((F0∧(CC⊻AA))⊻AA)
#else
	result = ((a & (b ^ c)) ^ c);
#endif
	return result;
}
//#define F3(x, y, z)   ((x | ~y) ^ z)
__device__ __forceinline__
uint32_t F3(const uint32_t x, const uint32_t y, const uint32_t z) {
	uint32_t result;
#if __CUDA_ARCH__ >= 500 && CUDA_VERSION >= 7050
	asm volatile ("lop3.b32 %0, %1, %2, %3, 0x59;" : "=r"(result) : "r"(x), "r"(y), "r"(z)); //0x59=((F0∨(¬CC))⊻AA)
#else
	result = ((x | ~y) ^ z);
#endif
	return result;
}
//#define F4(x, y, z)   (y ^ ((x ^ y) & z))
__device__ __forceinline__
uint32_t F4(const uint32_t x, const uint32_t y, const uint32_t z) {
	uint32_t result;
#if __CUDA_ARCH__ >= 500 && CUDA_VERSION >= 7050
	asm volatile ("lop3.b32 %0, %1, %2, %3, 0xE4;" : "=r"(result) : "r"(x), "r"(y), "r"(z)); //0xE4=(CC⊻((F0⊻CC)∧AA))
#else
	result = (y ^ ((x ^ y) & z));
#endif
	return result;
}
//#define F5(x, y, z)   (x ^ (y | ~z))
__device__ __forceinline__
uint32_t F5(const uint32_t x, const uint32_t y, const uint32_t z) {
	uint32_t result;
#if __CUDA_ARCH__ >= 500 && CUDA_VERSION >= 7050
	asm volatile ("lop3.b32 %0, %1, %2, %3, 0x2D;" : "=r"(result) : "r"(x), "r"(y), "r"(z)); //0x2D=(F0⊻(CC∨(¬AA)))
#else
	result = (x ^ (y | ~z));
#endif
	return result;
}

typedef unsigned int sph_u32;
//#define F1(x, y, z)   ((x) ^ (y) ^ (z))
//#define F2(x, y, z)   ((((y) ^ (z)) & (x)) ^ (z))
//#define F3(x, y, z)   (((x) | ~(y)) ^ (z))
//#define F4(x, y, z)   ((((x) ^ (y)) & (z)) ^ (y))
//#define F5(x, y, z)   ((x) ^ ((y) | ~(z)))



//#define ROTL    SPH_ROTL32
//#define ROTL(x, n) SPH_T32(((x) << (n)) | ((x) >> (32 - (n))))



//#define sRR(a, b, c, d, f, s, r, k)   do { \
//		a = ROTL(SPH_T32(a + f(b, c, d) + r + k), s); \
//	} while (0)

#define RR(a, b, c, d, f, s, r, k) { \
	a = ROTATE((a + k + r + f(b, c, d)), s); \
	c = ROTL32(c, 10); \
}

//#define sROUND1(a, b, c, d, f, s, r, k)  \
//	sRR(a ## 1, b ## 1, c ## 1, d ## 1, f, s, r, sK1 ## k)

#define ROUND1(a, b, c, d, f, s, r, k) \
	RR(a[0], b[0], c[0], d[0], f, s, r, c_K1[k])


//#define sROUND2(a, b, c, d, f, s, r, k)  \
	sRR(a ## 2, b ## 2, c ## 2, d ## 2, f, s, r, sK2 ## k)

#define ROUND2(a, b, c, d, f, s, r, k) \
	RR(a[1], b[1], c[1], d[1], f, s, r, c_K2[k])



#define RIPEMD128_ROUND_BODY(in, h) { \
	uint32_t A[2], B[2], C[2], D[2], E[2]; \
	uint32_t tmp; \
 \
	A[0] = A[1] = h[0]; \
	B[0] = B[1] = h[1]; \
	C[0] = C[1] = h[2]; \
	D[0] = D[1] = h[3]; \
 \
		ROUND1(A, B, C, D, F1, 11, in[ 0],  0); \
		ROUND1(D, A, B, C, F1, 14, in[ 1],  0); \
		ROUND1(C, D, A, B, F1, 15, in[ 2],  0); \
		ROUND1(B, C, D, A, F1, 12, in[ 3],  0); \
		ROUND1(A, B, C, D, F1,  5, in[ 4],  0); \
		ROUND1(D, A, B, C, F1,  8, in[ 5],  0); \
		ROUND1(C, D, A, B, F1,  7, in[ 6],  0); \
		ROUND1(B, C, D, A, F1,  9, in[ 7],  0); \
		ROUND1(A, B, C, D, F1, 11, in[ 8],  0); \
		ROUND1(D, A, B, C, F1, 13, in[ 9],  0); \
		ROUND1(C, D, A, B, F1, 14, in[10],  0); \
		ROUND1(B, C, D, A, F1, 15, in[11],  0); \
		ROUND1(A, B, C, D, F1,  6, in[12],  0); \
		ROUND1(D, A, B, C, F1,  7, in[13],  0); \
		ROUND1(C, D, A, B, F1,  9, in[14],  0); \
		ROUND1(B, C, D, A, F1,  8, in[15],  0); \
 \
		ROUND1(A, B, C, D, F2,  7, in[ 7],  1); \
		ROUND1(D, A, B, C, F2,  6, in[ 4],  1); \
		ROUND1(C, D, A, B, F2,  8, in[13],  1); \
		ROUND1(B, C, D, A, F2, 13, in[ 1],  1); \
		ROUND1(A, B, C, D, F2, 11, in[10],  1); \
		ROUND1(D, A, B, C, F2,  9, in[ 6],  1); \
		ROUND1(C, D, A, B, F2,  7, in[15],  1); \
		ROUND1(B, C, D, A, F2, 15, in[ 3],  1); \
		ROUND1(A, B, C, D, F2,  7, in[12],  1); \
		ROUND1(D, A, B, C, F2, 12, in[ 0],  1); \
		ROUND1(C, D, A, B, F2, 15, in[ 9],  1); \
		ROUND1(B, C, D, A, F2,  9, in[ 5],  1); \
		ROUND1(A, B, C, D, F2, 11, in[ 2],  1); \
		ROUND1(D, A, B, C, F2,  7, in[14],  1); \
		ROUND1(C, D, A, B, F2, 13, in[11],  1); \
		ROUND1(B, C, D, A, F2, 12, in[ 8],  1); \
 \
		ROUND1(A, B, C, D, F3, 11, in[ 3],  2); \
		ROUND1(D, A, B, C, F3, 13, in[10],  2); \
		ROUND1(C, D, A, B, F3,  6, in[14],  2); \
		ROUND1(B, C, D, A, F3,  7, in[ 4],  2); \
		ROUND1(A, B, C, D, F3, 14, in[ 9],  2); \
		ROUND1(D, A, B, C, F3,  9, in[15],  2); \
		ROUND1(C, D, A, B, F3, 13, in[ 8],  2); \
		ROUND1(B, C, D, A, F3, 15, in[ 1],  2); \
		ROUND1(A, B, C, D, F3, 14, in[ 2],  2); \
		ROUND1(D, A, B, C, F3,  8, in[ 7],  2); \
		ROUND1(C, D, A, B, F3, 13, in[ 0],  2); \
		ROUND1(B, C, D, A, F3,  6, in[ 6],  2); \
		ROUND1(A, B, C, D, F3,  5, in[13],  2); \
		ROUND1(D, A, B, C, F3, 12, in[11],  2); \
		ROUND1(C, D, A, B, F3,  7, in[ 5],  2); \
		ROUND1(B, C, D, A, F3,  5, in[12],  2); \
 \
		ROUND1(A, B, C, D, F4, 11, in[ 1],  3); \
		ROUND1(D, A, B, C, F4, 12, in[ 9],  3); \
		ROUND1(C, D, A, B, F4, 14, in[11],  3); \
		ROUND1(B, C, D, A, F4, 15, in[10],  3); \
		ROUND1(A, B, C, D, F4, 14, in[ 0],  3); \
		ROUND1(D, A, B, C, F4, 15, in[ 8],  3); \
		ROUND1(C, D, A, B, F4,  9, in[12],  3); \
		ROUND1(B, C, D, A, F4,  8, in[ 4],  3); \
		ROUND1(A, B, C, D, F4,  9, in[13],  3); \
		ROUND1(D, A, B, C, F4, 14, in[ 3],  3); \
		ROUND1(C, D, A, B, F4,  5, in[ 7],  3); \
		ROUND1(B, C, D, A, F4,  6, in[15],  3); \
		ROUND1(A, B, C, D, F4,  8, in[14],  3); \
		ROUND1(D, A, B, C, F4,  6, in[ 5],  3); \
		ROUND1(C, D, A, B, F4,  5, in[ 6],  3); \
		ROUND1(B, C, D, A, F4, 12, in[ 2],  3); \
 \
		ROUND2(A, B, C, D, F4,  8, in[ 5],  0); \
		ROUND2(D, A, B, C, F4,  9, in[14],  0); \
		ROUND2(C, D, A, B, F4,  9, in[ 7],  0); \
		ROUND2(B, C, D, A, F4, 11, in[ 0],  0); \
		ROUND2(A, B, C, D, F4, 13, in[ 9],  0); \
		ROUND2(D, A, B, C, F4, 15, in[ 2],  0); \
		ROUND2(C, D, A, B, F4, 15, in[11],  0); \
		ROUND2(B, C, D, A, F4,  5, in[ 4],  0); \
		ROUND2(A, B, C, D, F4,  7, in[13],  0); \
		ROUND2(D, A, B, C, F4,  7, in[ 6],  0); \
		ROUND2(C, D, A, B, F4,  8, in[15],  0); \
		ROUND2(B, C, D, A, F4, 11, in[ 8],  0); \
		ROUND2(A, B, C, D, F4, 14, in[ 1],  0); \
		ROUND2(D, A, B, C, F4, 14, in[10],  0); \
		ROUND2(C, D, A, B, F4, 12, in[ 3],  0); \
		ROUND2(B, C, D, A, F4,  6, in[12],  0); \
 \
		ROUND2(A, B, C, D, F3,  9, in[ 6],  1); \
		ROUND2(D, A, B, C, F3, 13, in[11],  1); \
		ROUND2(C, D, A, B, F3, 15, in[ 3],  1); \
		ROUND2(B, C, D, A, F3,  7, in[ 7],  1); \
		ROUND2(A, B, C, D, F3, 12, in[ 0],  1); \
		ROUND2(D, A, B, C, F3,  8, in[13],  1); \
		ROUND2(C, D, A, B, F3,  9, in[ 5],  1); \
		ROUND2(B, C, D, A, F3, 11, in[10],  1); \
		ROUND2(A, B, C, D, F3,  7, in[14],  1); \
		ROUND2(D, A, B, C, F3,  7, in[15],  1); \
		ROUND2(C, D, A, B, F3, 12, in[ 8],  1); \
		ROUND2(B, C, D, A, F3,  7, in[12],  1); \
		ROUND2(A, B, C, D, F3,  6, in[ 4],  1); \
		ROUND2(D, A, B, C, F3, 15, in[ 9],  1); \
		ROUND2(C, D, A, B, F3, 13, in[ 1],  1); \
		ROUND2(B, C, D, A, F3, 11, in[ 2],  1); \
 \
		ROUND2(A, B, C, D, F2,  9, in[15],  2); \
		ROUND2(D, A, B, C, F2,  7, in[ 5],  2); \
		ROUND2(C, D, A, B, F2, 15, in[ 1],  2); \
		ROUND2(B, C, D, A, F2, 11, in[ 3],  2); \
		ROUND2(A, B, C, D, F2,  8, in[ 7],  2); \
		ROUND2(D, A, B, C, F2,  6, in[14],  2); \
		ROUND2(C, D, A, B, F2,  6, in[ 6],  2); \
		ROUND2(B, C, D, A, F2, 14, in[ 9],  2); \
		ROUND2(A, B, C, D, F2, 12, in[11],  2); \
		ROUND2(D, A, B, C, F2, 13, in[ 8],  2); \
		ROUND2(C, D, A, B, F2,  5, in[12],  2); \
		ROUND2(B, C, D, A, F2, 14, in[ 2],  2); \
		ROUND2(A, B, C, D, F2, 13, in[10],  2); \
		ROUND2(D, A, B, C, F2, 13, in[ 0],  2); \
		ROUND2(C, D, A, B, F2,  7, in[ 4],  2); \
		ROUND2(B, C, D, A, F2,  5, in[13],  2); \
 \
		ROUND2(A, B, C, D, F1, 15, in[ 8],  3); \
		ROUND2(D, A, B, C, F1,  5, in[ 6],  3); \
		ROUND2(C, D, A, B, F1,  8, in[ 4],  3); \
		ROUND2(B, C, D, A, F1, 11, in[ 1],  3); \
		ROUND2(A, B, C, D, F1, 14, in[ 3],  3); \
		ROUND2(D, A, B, C, F1, 14, in[11],  3); \
		ROUND2(C, D, A, B, F1,  6, in[15],  3); \
		ROUND2(B, C, D, A, F1, 14, in[ 0],  3); \
		ROUND2(A, B, C, D, F1,  6, in[ 5],  3); \
		ROUND2(D, A, B, C, F1,  9, in[12],  3); \
		ROUND2(C, D, A, B, F1, 12, in[ 2],  3); \
		ROUND2(B, C, D, A, F1,  9, in[13],  3); \
		ROUND2(A, B, C, D, F1, 12, in[ 9],  3); \
		ROUND2(D, A, B, C, F1,  5, in[ 7],  3); \
		ROUND2(C, D, A, B, F1, 15, in[10],  3); \
		ROUND2(B, C, D, A, F1,  8, in[14],  3); \
 \
	tmp  = h[1] + C[0] + D[1]; \
	h[1] = h[2] + D[0] + A[1]; \
	h[2] = h[3] + A[0] + B[1]; \
	h[3] = h[4] + B[0] + C[1]; \
	h[0] = tmp; \
	} 

// END OF RIPEMD MACROS----------------------------------------------------------------------


__global__ /* __launch_bounds__(256, 6) */
void x17_ripemd160_gpu_hash_64(const uint32_t threads, uint64_t *g_hash, const int outlen)
{
	const uint32_t thread = (blockDim.x * blockIdx.x + threadIdx.x);
	if (thread < threads)
	{
		const uint64_t hashPosition = thread * 8U;
		uint64_t *pHash = &g_hash[hashPosition];
		uint32_t h[5];
		
			uint64_t dat[8];
#pragma unroll 5
			for (int i = 0; i<5; i++)
				h[i] = c_IV[i];

#pragma unroll
		for (int i = 0; i<8; i++) {
			dat[i] = pHash[i];
		}

		///////// input big /////////////////////

		uint32_t buf[32];
		RIPEMD128_ROUND_BODY(dat, h);



		if (outlen == 512) {
			pHash[4] = 0; //hash.h8[4];
			pHash[5] = 0; //hash.h8[5];
			pHash[6] = 0; //hash.h8[6];
			pHash[7] = 0; //hash.h8[7];
		}
	}
}

__host__
void x17_ripemd160_cpu_init(int thr_id, uint32_t threads)
{
}

__host__
void x17_ripemd160_cpu_hash_64(int thr_id, uint32_t threads, uint32_t startNounce, uint32_t *d_hash, const int outlen)
{
	const uint32_t threadsperblock = 256;

	dim3 grid((threads + threadsperblock - 1) / threadsperblock);
	dim3 block(threadsperblock);

	x17_ripemd160_gpu_hash_64 << <grid, block >> > (threads, (uint64_t*)d_hash, outlen);
}
