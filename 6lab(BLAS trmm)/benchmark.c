#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cblas.h>
#include <math.h>

#include "trmm.h"

#define N 2000

double get_time()
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}

void fill_matrix(double *M, int n)
{
    for(int i = 0; i < n*n; i++)
        M[i] = (double)rand() / RAND_MAX;
}

int main()
{
    srand(0);

    double *A  = malloc(N*N*sizeof(double));
    double *B  = malloc(N*N*sizeof(double));
    double *B2 = malloc(N*N*sizeof(double));
    double *B0 = malloc(N*N*sizeof(double));

    fill_matrix(A, N);
    fill_matrix(B0, N);

    printf("Matrix size %d x %d\n\n", N, N);

    double geom = 1.0;

    for(int run = 0; run < 10; run++)
    {
        for(int i = 0; i < N*N; i++)
        {
            B[i]  = B0[i];
            B2[i] = B0[i];
        }

        double start = get_time();

        dtrmm(0,0,0,0,N,N,1.0,A,N,B,N);

        double end = get_time();
        double my_time = end - start;

        start = get_time();

        cblas_dtrmm(CblasRowMajor,CblasLeft,CblasLower,CblasNoTrans,CblasNonUnit,N,N,1.0,A,N,B2,N);

        end = get_time();
        double blas_time = end - start;

        double perf = (blas_time / my_time) * 100.0;

        geom *= perf;

        printf("Run %d\n", run + 1);
        printf("my_trmm:    %f sec\n", my_time);
        printf("openblas:   %f sec\n", blas_time);
        printf("performance %.2f %%\n\n", perf);
    }

    double geo_mean = pow(geom, 1.0/10.0);

    printf("Average geometric performance: %.2f %%\n", geo_mean);

    free(A);
    free(B);
    free(B2);
    free(B0);

    return 0;
}
