#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<xmmintrin.h>

void transpose(double *__restrict copy_to, const int lda, const double *__restrict copy_from);
void populate(const int lda, double *__restrict copy_to, const int M, const int K, const double *__restrict copy_from);
const char* dgemm_desc = "My dgemm.";
void basic_dgemm2(const int lda, const int M, const int N, const int K,
                 const double *__restrict A, const double *__restrict B, double *__restrict C, const double *Ap, const double *Bp);

#ifndef BLOCK_SIZE
#define BLOCK_SIZE ((int) 256)
#define MAT_SIZE ((int)BLOCK_SIZE * BLOCK_SIZE)
#endif
/*
  A is M-by-K
  B is K-by-N
  C is M-by-N

  lda is the leading dimension of the matrix (the M of square_dgemm).
*/
void basic_dgemm2(const int lda, const int M, const int N, const int K,
                 const double *__restrict A, const double *__restrict B, double *__restrict C, const double *Ap, const double *Bp)
{
    int i, j, k, idx;
    double cij;
    double ctmp[2] __attribute__((aligned(16)));
    __m128d c;
    for (j = 0; j < N; ++j) {
        for (i = 0; i < M; ++i) {
            idx = j*lda+i;
            cij = C[idx];
            c = _mm_set_pd(cij, 0.);
            for (k = 0; k < K; k=k+2) {
                c = _mm_add_pd(c, _mm_mul_pd(_mm_load_pd(A+i*BLOCK_SIZE+k), _mm_load_pd(B+j*BLOCK_SIZE+k)));
            }        
            _mm_store_pd(ctmp, c);
            cij = ctmp[1] + ctmp[0];
            idx = j*lda+i;
            C[idx] = cij;
        }
    }
}
void basic_dgemm(const int lda, const int M, const int N, const int K,
                 const double *__restrict A, const double *__restrict B, double *__restrict C)
{
    int i, j, k, idx;
    double a1, a2, b1, b2, cij;
    double ctmp[2] __attribute__((aligned(16)));
    const int klim = (K&1 == 1)? K-1 : K;
    __m128d c;
    for (j = 0; j < N; ++j) {
        for (i = 0; i < M; ++i) {
            idx = j*lda+i;
            cij = C[idx];
            c = _mm_set_pd(cij, 0.);
            for (k = 0; k < klim; k=k+2) {
					 c = _mm_add_pd(c, _mm_mul_pd(_mm_set_pd(A[i*lda+k], A[i*lda+k+1]), _mm_set_pd(B[j*lda+k],B[j*lda+k+1])));
            }
            if (K&1==1) {
					 c = _mm_add_pd(c, _mm_mul_pd(_mm_set_pd(A[i*lda+K-1], 0.), _mm_set_pd(B[j*lda+K-1], 0.)));
            }
            _mm_store_pd(ctmp, c);
            cij = ctmp[1] + ctmp[0];
            idx = j*lda+i;
            C[idx] = cij;
        }
    }
}

void do_block(const int lda,
              const double *A, const double *B, double *C,
              const int i, const int j, const int k)
{
  
    const int M = (i+BLOCK_SIZE > lda? lda-i : BLOCK_SIZE);
    const int N = (j+BLOCK_SIZE > lda? lda-j : BLOCK_SIZE);
    const int K = (k+BLOCK_SIZE > lda? lda-k : BLOCK_SIZE);
     
    if (lda > 190) {
        static double copy_A[MAT_SIZE] __attribute__((aligned(16)));
        static double copy_B[MAT_SIZE] __attribute__((aligned(16)));
        populate(lda, copy_A, M, K, A + k + i*lda);
        populate(lda, copy_B, N, K, B + k + j*lda);
        basic_dgemm2(lda, M, N, K, copy_A, copy_B, C + i + j*lda, A + k + i*lda, B + k + j*lda);
    } else {
        basic_dgemm(lda, M, N, K,
                A + k + i*lda, B + k + j*lda, C + i + j*lda);
    }
}

void transpose(double *__restrict copy_to, const int lda, const double *__restrict copy_from) {
    int c, r, idx;
    double from;
    for (c = 0; c < lda; c++) {
        for (r = 0; r < lda; r++) {
              idx = r * lda + c;
              from = copy_from[c * lda + r];
              copy_to[idx] = from;
        }
    }
}
 
void populate(const int lda, double *__restrict copy_to, const int M, const int K, const double *__restrict copy_from) {
    int i, j;
    double from;
    
    for (i = 0; i < BLOCK_SIZE; ++i) {
        for (j = 0; j < BLOCK_SIZE; ++j) {
            if (i < K && j < M) {
                from = copy_from[j*lda + i];
                copy_to[j*BLOCK_SIZE + i] = from;
            } else {
                copy_to[j*BLOCK_SIZE + i] = 0.;
            }
        }
    }
}

void square_dgemm(const int M, const double *__restrict A, const double *__restrict B, double *__restrict C)
{
    const int leftover = M%BLOCK_SIZE;
    const int n_blocks = M / BLOCK_SIZE + (leftover? 1 : 0);
    // copy matrix
    const int lda = BLOCK_SIZE * n_blocks;
    const int tot_num_elems = M * M;
    double transposed_A[tot_num_elems];
    transpose(transposed_A, M, A); 
     int bi, bj, bk;
    for (bi = 0; bi < n_blocks; ++bi) {
        const int i = bi * BLOCK_SIZE;
        for (bj = 0; bj < n_blocks; ++bj) {
            const int j = bj * BLOCK_SIZE;
            for (bk = 0; bk < n_blocks; ++bk) {
                const int k = bk * BLOCK_SIZE;
                do_block(M, transposed_A, B, C, j, i, k);
            }
        }
    }
}

