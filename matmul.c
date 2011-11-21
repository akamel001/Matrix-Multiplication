/*
  Ideally, you won't need to change this file.  You may want to change
  a few settings to speed debugging runs, but remember to change back
  to the original settings during final testing.

  The output format: "Size: %u\tmflop/s: %g\n"

  These hands have touched the file:
    David Bindel
    David Garmire
    Jason Riedy
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <sys/types.h>
#include <sys/resource.h>
#include <unistd.h>
#include <time.h>

#include "timing.h"

#ifndef COMPILER
#  define COMPILER "unknown"
#endif
#ifndef FLAGS
#  define FLAGS "unknown"
#endif

/*
  Your function _MUST_ have the following signature:
*/
extern const char* dgemm_desc;
extern void square_dgemm();

/*
  We try to run enough iterations to get reasonable timings.  The matrices
  are multiplied at least MIN_RUNS times.  If that doesn't take MIN_SECS
  seconds, then we double the number of iterations and try again.

  You may want to modify these to speed debugging...
*/
#define MIN_RUNS 4
/* #define MIN_SECS 1.0 */
#define MIN_SECS 0.25

/*
  Note the strange sizes...  You'll see some interesting effects
  around some of the powers-of-two.
*/
const int test_sizes[] = {
    31, 32, 96, 97, 127, 128, 129, 191, 192, 229,
#if defined(DEBUG_RUN)
# define MAX_SIZE 229u
#else
    255, 256, 257, 319, 320, 321, 417, 479, 480, 511, 512, 639, 640,
    767, 768, 769,
# define MAX_SIZE 769u
#endif
};

#define N_SIZES (sizeof (test_sizes) / sizeof (int))


/* --
 * Initialize A to random numbers (A is MAX_SIZE * MAX_SIZE)
 */
void matrix_init(double *A)
{
    int i;
    for (i = 0; i < MAX_SIZE*MAX_SIZE; ++i) 
        A[i] = drand48();
}


/* --
 * Zero out C (which is MAX_SIZE * MAX_SIZE)
 */
void matrix_clear(double *C)
{
    memset(C, 0, MAX_SIZE * MAX_SIZE * sizeof(double));
}


/* --
 * Check that C = A*B to within roundoff error.
 *
 * We use the fact that dot products satisfy the error bound
 *
 *   float(sum a_i * b_i) = sum a_i * b_i * (1 + delta_i)
 *
 * where delta_i <= n * epsilon.  In order to check your matrix
 * multiply, we compute each element in turn and make sure that
 * your product is within three times the given error bound.
 * We make it three times because there are three sources of
 * error:
 *
 *  - the roundoff error in your multiply
 *  - the roundoff error in our multiply
 *  - the roundoff error in computing the error bound
 *
 *  That last source of error is not so significant, but that's a
 *  story for another day.
 */
void validate_dgemm(const int M, const double *A, const double *B, double *C)
{
    int i, j, k;

    matrix_clear(C);
    square_dgemm(M, A, B, C);

    for (i = 0; i < M; ++i) {
        for (j = 0; j < M; ++j) {
            double dotprod = 0;
            double errorbound = 0;
            double err;
            for (k = 0; k < M; ++k) {
                double prod = A[k*M + i] * B[j*M + k];
                dotprod += prod;
                errorbound += fabs(prod);
            }
            errorbound *= (M * DBL_EPSILON);
            err = fabs(C[j*M + i] - dotprod);
            if (err > 3*errorbound) {
                printf("Matrix multiply failed.\n");
                printf("C(%d,%d) should be %lg, was %lg\n", i, j,
                       dotprod, C[j*M + i]);
                printf("Error of %lg, acceptable limit %lg\n",
                       err, 3*errorbound);
                exit(-1);
            }
        }
    }
}


/* --
 * Compute a MFlop/s rate for C += A*B.
 *
 * The code runs the multiplication repeatedly in a loop MIN_RUNS times,
 * then doubles the loop time if it did not take MIN_SECS to perform the
 * run.  This helps us get around the limits of timer resolution.
 */
double time_dgemm(const int M, const double *A, const double *B, double *C)
{
    timing_t start, finish;
    double mflops, mflop_s;
    double secs = -1.0;
    int num_iterations = MIN_RUNS;
    int i;
    while (secs < MIN_SECS) {
        matrix_clear(C);
        get_time(&start);
        for (i = 0; i < num_iterations; ++i) {
            square_dgemm(M, A, B, C);
        }
        get_time(&finish);
        secs = timespec_diff(start, finish);
        mflops  = 2.0 * num_iterations * M * M * M / 1.0e6;
        mflop_s = mflops/secs;
        num_iterations *= 2;
    }
    return mflop_s;
}


int main()
{
    int i;
    double mflop_s;
    double* A __attribute__((aligned(16))) = (double*) malloc(MAX_SIZE * MAX_SIZE * sizeof(double));
    double* B __attribute__((aligned(16))) = (double*) malloc(MAX_SIZE * MAX_SIZE * sizeof(double));
    double* C __attribute__((aligned(16))) = (double*) malloc(MAX_SIZE * MAX_SIZE * sizeof(double));

    matrix_init(A);
    matrix_init(B);

    printf("Compiler:\t%s\nOptions:\t%s\nDescription:\t%s\n\n",
           COMPILER, FLAGS, dgemm_desc);

    for (i = 0; i < N_SIZES; ++i) {
        const int M = test_sizes[i];
        validate_dgemm(M, A, B, C);
        mflop_s = time_dgemm(M, A, B, C);
        printf("Size: %u\tmflop/s: %lg\n", M, mflop_s);
    }

    free(C);
    free(B);
    free(A);

    return 0;
}

