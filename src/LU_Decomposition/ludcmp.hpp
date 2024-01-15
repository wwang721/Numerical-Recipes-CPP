#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

/*
 * LU decomposition of a matrix.
 *  @param a: given a matrix a[0..n-1][0..n-1], this routine replaces it by the LU decomposition of a rowwise permutation of itself. 
 *  @param n: is input as the dimension of the matrix a.
 *  @param indx: indx[0..n-1] is an output vector that records the row permutation effected by the partial pivoting;
 *  @returns d: +1 or -1 depending on whether the number of row interchanges was even or odd, respectively.
 * 
 * Note that a is also output. This routine is used in combination with lubksb to solve linear equations or invert a matrix.
 * Adapted from Numerical Recipes in C, 2nd Ed., Chapt. 2.3, p. 46.
 */
#define TINY 1.0e-20
int ludcmp(std::vector<std::vector<double>> &a, int n, std::vector<int> &indx) {
    int i, imax, j, k, d;
    double big, dum, sum, temp;
    double *vv = (double *)malloc(n * sizeof(double));

    d = 1;
    for (i = 0; i < n; i++) {
        big = 0.0;
        for (j = 0; j < n; j++)
            if ((temp = fabs(a[i][j])) > big) big = temp;

        if (big == 0.0) {
            fprintf(stderr, "Singular matrix in routine ludcmp");
            exit(1);
        }
        vv[i] = 1.0 / big;
    }

    for (j = 0; j < n; j++) {
        for (i = 0; i < j; i++) {
            sum = a[i][j];
            for (k = 0; k < i; k++) sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
        }
        big = 0.0;
        for (i = j; i < n; i++) {
            sum = a[i][j];
            for (k = 0; k < j; k++) sum -= a[i][k] * a[k][j];
            a[i][j] = sum;
            if ((dum = vv[i] * fabs(sum)) >= big) {
                big = dum;
                imax = i;
            }
        }
        if (j != imax) {
        for (k = 0; k < n; k++) {
            dum = a[imax][k];
            a[imax][k] = a[j][k];
            a[j][k] = dum;
        }
        d = -d;
        vv[imax] = vv[j];
        }
        indx[j] = imax;

        if (a[j][j] == 0.0) a[j][j]=TINY; 
        // If the pivot element is zero the matrix is singular (at least to the precision of the algorithm).
        // For some applications on singular matrices, it is desirable to substitute TINY for zero

        if (j != n - 1) {
            dum = (1.0 / a[j][j]);
            for (i = j + 1; i < n; i++) a[i][j] *= dum;
        }
    }
    free(vv);

    return d;
}
#undef TINY


/*
 * Solves the set of n linear equations A X = B. (lubksb stands for LU Back Substitution)
 *  @param a: here a[0..n-1][0..n-1] is input, not as the matrix A but rather as its LU decomposition, determined by the routine ludcmp.
 *  @param n: is input as the dimension of A and the vector B.
 *  @param indx: indx[0..n-1] is input as the permutation vector returned by ludcmp.
 *  @param b: b[0..n-1] is input as the right-hand side vector B, and returns with the solution vector X. 
 
 * a, n, and indx are not modified by this routine and can be left in place for successive calls with different right-hand sides b.
 * This routine takes into account the possibility that b will begin with many zero elements, so it is efficient for use in matrix inversion. 
 * Adapted from Numerical Recipes in C, 2nd Ed., Chapt. 2.3, p. 47.
 */
void lubksb(std::vector<std::vector<double>> &a, int n, std::vector<int> &indx, std::vector<double> &b) {
    int i, ii = 0, ip, j;
    double sum;

    for (i = 0; i < n; i++) {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];

        if (ii)
            for (j = ii - 1; j < i; j++) sum -= a[i][j] * b[j];
        else if (sum)
            ii = i + 1;

        b[i] = sum;
    }

    for (i = n - 1; i >= 0; i--) {
        sum = b[i];
        for (j = i + 1; j < n; j++) sum -= a[i][j] * b[j];
        b[i] = sum / a[i][i];
    }
}