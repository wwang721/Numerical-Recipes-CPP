# Numerical Recipes in C++

Adapt some utilities from the famous book [*"Numerical Recipes in C"*](http://kfes-16.karlov.mff.cuni.cz/~standa/nc/www.library.cornell.edu/nr/cbookcpdf.html) (Cambridge University Press make it online now!) to C++.


## Installation

1. CMake version >= 3.0 needed ![CMake](https://img.shields.io/badge/CMake-v3.0-brightgreen.svg)

2. Navigate to the [build](/build) folder by `cd build`, then build the source code by `cmake ..` and `make`

3. Now you can find the executable file with .exe extension in the `build` folder

## Contents

1. [**LU decomposition**](https://courses.physics.illinois.edu/cs357/sp2020/notes/ref-9-linsys.html) ($\mathbf{A=LU}$) to solve linear equations $\mathbf{Ax=b}$ or calculate the inverse $\mathbf{A}^{-1}$. See header file [ludcmp.hpp](/src/LU_Decomposition) [Chapt. 2.3].

2. [**Fredholm integral equation**](https://en.wikipedia.org/wiki/Fredholm_integral_equation) of the second kind $f(t)=\lambda\int_a^bK (t,s)f(s)\mathrm{d}s+g(t)$ can be tranformed into a matrix representation $(1-\lambda\tilde{\mathbf{K}})\cdot \mathbf{f}=\mathbf{g}$. See header file [fred2.hpp](/src/Fredholm) [Chapt. 18.1].



## License

This project is licensed under a MIT License ([to see why](https://choosealicense.com/)). Please see the [LICENSE](/LICENSE) file for details.