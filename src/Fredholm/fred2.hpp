#include <vector>
#include "ludcmp.hpp"

// Fredholm equation of the second kind
// Solve (1-omk)f=g and the result is stored in f
// Adapted from Numerical Recipes in C, 2nd Ed., Chapt. 18.1, p. 793.
void fred2(int n, std::vector<std::vector<double>> &omk, std::vector<double> &g, std::vector<double> &f) {

    std::vector<int> indx(n);

    for (int i=0; i<n; ++i) 
    { 
        for (int j=0; j<n; ++j) 
            omk[i][j]=(double)(i == j) - omk[i][j];
        f[i]=g[i];
    }

    int d = ludcmp(omk, n, indx); 
    lubksb(omk, n, indx, f);

}