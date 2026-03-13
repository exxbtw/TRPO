#ifndef TRMM_H
#define TRMM_H

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
);

void strmm(
    int side,
    int uplo,
    int trans,
    int diag,
    int m,
    int n,
    float alpha,
    float *A,
    int lda,
    float *B,
    int ldb
);

#endif