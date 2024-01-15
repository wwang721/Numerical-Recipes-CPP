#include <iostream>
#include <vector>
#include "ludcmp.hpp"

int main() {
    // Create a 3x3 matrix
    int n = 3;
    std::vector<std::vector<double>> matrix_A(n, std::vector<double>(n));

    matrix_A[0][0] = 1;
    matrix_A[0][1] = 2;
    matrix_A[1][1] = 2;
    matrix_A[1][2] = 1;
    matrix_A[2][0] = 1;
    matrix_A[2][1] = 1;
    matrix_A[2][2] = 2; 


    std::vector<int> indx(n);
    int d = ludcmp(matrix_A, n, indx);  // LU decomposition result is stored in matrix_A and indx

    /*
     * To solve linear equations Ax = b, we only need to do the following:
     *
     * std::vector<double> b {1, 2, 3};
     * lubksb(matrix_A, n, indx, b);    // solution x stored in b, but note A here must be the LU decomposition of A.
     */

    // We can also use lubksb to get the inverse matrix
    std::vector<std::vector<double>> inverse_matrix_A(n, std::vector<double>(n));
    // Using LU decomposition (ludcmp + lubksb's) to get matrix inversion is a ~O(N^3) operation, 
    // so this method is equally efficient as Gauss-Jordan elimination method.
    for(int j = 0; j < n; ++j) {
        std::vector<double> b (n);
        b[j] = 1.;
        lubksb(matrix_A, n, indx, b);
        // Matrix A and indx are not modified by this routine and can be left in place for successive calls with different right-hand sides b.
        for(int i = 0; i < n; i++)
            inverse_matrix_A[i][j] = b[i];
    }

    // Print the matrix
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << inverse_matrix_A[i][j] << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}
