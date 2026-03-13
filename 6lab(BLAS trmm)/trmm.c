#include "trmm.h"
#include <stdlib.h>

void dtrmm(
    int side,
    int uplo,
    int trans,
    int diag,
    int m,
    int n,
    double alpha,
    double *A,
    int lda,
    double *B,
    int ldb
)
{
    if(side == 0 && uplo == 0 && trans == 0)
    {
        double *temp = malloc(n * sizeof(double));

        for(int i = 0; i < m; i++)
        {
            int k0 = 0;

            double a;

            if(diag && i == 0)
                a = 1.0;
            else
                a = A[i*lda + 0];

            a *= alpha;

            double *B0 = &B[0*ldb];

            for(int j = 0; j < n; j++)
                temp[j] = a * B0[j];

            for(int k = 1; k <= i; k++)
            {
                if(diag && i == k)
                    a = 1.0;
                else
                    a = A[i*lda + k];

                a *= alpha;

                double *Bk = &B[k*ldb];

                for(int j = 0; j < n; j++)
                    temp[j] += a * Bk[j];
            }

            double *Bi = &B[i*ldb];

            for(int j = 0; j < n; j++)
                Bi[j] = temp[j];
        }

        free(temp);
    }
}


void strmm(
    int side,
    int uplo,
    int trans,
    int diag,
    int m,
    int n,
    float alpha,
    float *restrict A,
    int lda,
    float *restrict B,
    int ldb
)
{
    if(side == 0 && uplo == 0 && trans == 0)
    {
        float *temp = malloc(n * sizeof(float));

        for(int i = 0; i < m; i++)
        {
            for(int j = 0; j < n; j++)
                temp[j] = 0.0f;

            float *Bi = &B[i*ldb];

            for(int k = 0; k <= i; k++)
            {
                float a;

                if(diag && i == k)
                    a = 1.0f;
                else
                    a = A[i*lda + k];

                a *= alpha;

                float *Bk = &B[k*ldb];

                for(int j = 0; j < n; j++)
                    temp[j] += a * Bk[j];
            }

            for(int j = 0; j < n; j++)
                Bi[j] = temp[j];
        }

        free(temp);
    }
}