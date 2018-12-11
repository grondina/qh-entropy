#include <assert.h>
#include <cblas.h>
#include <lapacke.h>
#include <lapacke_utils.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "data.h"
#include "util.h"

static void print_matrix(int N, int M, double (*A)[M])
{
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            char *s = (j == M-1) ? "\n" : "  ";
            printf("% -8.4f%s", A[i][j], s);
        }
    }
}


static double det3x3(double *A)
{
        return (A[0]*A[4]*A[8] + A[1]*A[5]*A[6] + A[2]*A[3]*A[7] - A[2]*A[4]*A[6] - A[1]*A[3]*A[8] - A[0]*A[5]*A[7]);
}

int cmpint(const void *p1, const void *p2)
{
    int a = *((int *)p1);
    int b = *((int *)p2);
    return a - b;
}

static double distance(double x1, double y1, double z1,
                       double x2, double y2, double z2)
{
    return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2));
}

double gyration(struct molecule *mol)
{
    double sum = 0;
    double xi, yi, zi;
    double xj, yj, zj;

    for (int i = 0; i < mol->m; ++i) {

        xi = mol->R[i][0];
        yi = mol->R[i][1];
        zi = mol->R[i][2];

        for (int j = 0; j < mol->m; ++j) {

            xj = mol->R[j][0];
            yj = mol->R[j][1];
            zj = mol->R[j][2];

            double dist = distance(xi, yi, zi, xj, yj, zj);
            sum += dist*dist;
        }
    }

    return sum/(2.0 * mol->m * mol->m);
}

void kabsch(struct molecule *mol, struct molecule *ref)
{
    if (mol == NULL)
        return;

    if (ref == NULL)
        return;

    assert(mol->n == ref->n);
    int N = mol->n;

    double ( *P)[3] = mol->R;
    double ( *Q)[3] = ref->R;
    double ( *H)[3] = malloc(sizeof(double[3][3]));
    double ( *U)[3] = malloc(sizeof(double[3][3]));
    double ( *S)[3] = malloc(sizeof(double[3][3]));
    double ( *V)[3] = malloc(sizeof(double[3][3]));
    double ( *X)[3] = malloc(sizeof(double[3][3]));
    double ( *M)[3] = malloc(sizeof(double[3][3]));
    double ( *R)[3] = malloc(sizeof(double[3][3]));
    double ( *Y)[3] = malloc(sizeof(double[N][3]));
    double (*Pt)[N] = malloc(sizeof(double[3][N]));
    double (*Vt)[3] = malloc(sizeof(double[3][3]));
    double (*Ut)[3] = malloc(sizeof(double[3][3]));

    assert(H  != NULL);
    assert(U  != NULL);
    assert(S  != NULL);
    assert(V  != NULL);
    assert(X  != NULL);
    assert(M  != NULL);
    assert(R  != NULL);
    assert(Y  != NULL);
    assert(Pt != NULL);
    assert(Vt != NULL);
    assert(Ut != NULL);

    double *pP  = &( P[0][0]);
    double *pQ  = &( Q[0][0]);
    double *pH  = &( H[0][0]);
    double *pU  = &( U[0][0]);
    double *pS  = &( S[0][0]);
    double *pV  = &( V[0][0]);
    double *pX  = &( X[0][0]);
    double *pM  = &( M[0][0]);
    double *pR  = &( R[0][0]);
    double *pY  = &( Y[0][0]);
    double *pPt = &(Pt[0][0]);
    double *pVt = &(Vt[0][0]);
    double *pUt = &(Ut[0][0]);

    printf("-------------------------------\n");

    printf("P = \n");
    print_matrix(N, 3, P);

    printf("Q = \n");
    print_matrix(N, 3, Q);

    /* Get Pt = trans(P) */
    LAPACKE_dge_trans(LAPACK_ROW_MAJOR, N, 3, pP, 3, pPt, N);

    printf("Pt = \n");
    print_matrix(3, N, Pt);

    /* Get H = Pt*Q */
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, N, 1.0, pPt, N, pQ, 3, 0, pH, 3);
    printf("H = \n");
    print_matrix(3, 3, H);

    /* Perform SVD: H = U*S*Vt */
    double superb[10];
    LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A', 3, 3, pH, 3, pS, pU, 3, pVt, 3, superb);

    /* Get V = trans(Vt) */
    LAPACKE_dge_trans(LAPACK_ROW_MAJOR, 3, 3, pVt, 3, pV, 3);

    /* Get Ut = trans(U) */
    LAPACKE_dge_trans(LAPACK_ROW_MAJOR, 3, 3, pUt, 3, pU, 3);

    /* Get X = V*Ut */
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0, pV, 3, pUt, 3, 0, pX, 3);

    /* Get d = det(V*Ut) */
    double d = det3x3(pX);

    /* Build middle matrix M */
    M[0][0] = 1;
    M[0][1] = 0;
    M[0][2] = 0;
    M[1][0] = 0;
    M[1][1] = 1;
    M[1][2] = 0;
    M[2][0] = 0;
    M[2][1] = 0;
    M[2][2] = d;

    /* Calculate X = V*M */
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0, pV, 3, pM, 3, 0, pX, 3);

    /* Calculate R = X*Ut */
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0, pX, 3, pUt, 3, 0, pR, 3);

    /* Apply R to P: R*Pt -> Y */
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, N, 3, 1.0, pR, 3, pPt, N, 0, pY, N);

    /* Copy rotated coordinates back Y -> P */
    memcpy(mol->R, pY, (N * 3 * sizeof(double)));

    /* Free allocated memory */
    free(H);
    free(S);
    free(V);
    free(X);
    free(M);
    free(R);
    free(Y);
    free(Pt);
    free(Vt);
    free(Ut);
}

void remove_com(struct molecule *mol)
{
    double xcom = 0;
    double ycom = 0;
    double zcom = 0;
    double mass = 0;

    for (int j = 0; j < mol->m; ++j) {
        xcom += (mol->R[j][0] * mol->mass[j]);
        ycom += (mol->R[j][1] * mol->mass[j]);
        zcom += (mol->R[j][2] * mol->mass[j]);
        mass += mol->mass[j];
    }
    xcom /= mass;
    ycom /= mass;
    zcom /= mass;

    for (int j = 0; j < mol->m; ++j) {
        mol->R[j][0] -= xcom;
        mol->R[j][1] -= ycom;
        mol->R[j][2] -= zcom;
    }
}
