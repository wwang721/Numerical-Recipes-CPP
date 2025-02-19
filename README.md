# Numerical Recipes in C++

Adapt some utilities from the famous book [*"Numerical Recipes in C"*](http://kfes-16.karlov.mff.cuni.cz/~standa/nc/www.library.cornell.edu/nr/cbookcpdf.html) (Cambridge University Press make it online now!) to C++.


## Installation

1. CMake version >= 3.0 needed ![CMake](https://img.shields.io/badge/CMake-v3.0-brightgreen.svg)

2. Navigate to the [build](/build) folder by `cd build`, then build the source code by `cmake ..` and `make`

3. Now you can find the executable file with .exe extension in the `build` folder

## Contents

1. [**LU decomposition**](https://courses.physics.illinois.edu/cs357/sp2020/notes/ref-9-linsys.html) ($\mathbf{A=LU}$) to solve linear equations $\mathbf{Ax=b}$ or calculate the inverse $\mathbf{A}^{-1}$. See header file [ludcmp.hpp](https://github.com/wwang721/Numerical-Recipes-CPP/tree/main/src/LU_Decomposition/ludcmp.hpp) [Chapt. 2.3].

2. [**Fredholm integral equation**](https://en.wikipedia.org/wiki/Fredholm_integral_equation) of the second kind $f(t)=\lambda\int_a^bK (t,s)f(s)\mathrm{d}s+g(t)$ can be tranformed into a matrix representation $(1-\lambda\tilde{\mathbf{K}})\cdot \mathbf{f}=\mathbf{g}$. See header file [fred2.hpp](https://github.com/wwang721/Numerical-Recipes-CPP/tree/main/src/Fredholm/fred2.hpp) [Chapt. 18.1].

* Details are in the book's chapters, denoted by square brackets.

## Project information & Citation
[![DOI](https://zenodo.org/badge/743327087.svg)](https://zenodo.org/doi/10.5281/zenodo.11116966)

## License

This project is licensed under the [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0) license (CC-BY-4.0), allowing re-distribution and re-use of the licensed work on the condition that the creator is appropriately credited. Please see the [LICENSE](/LICENSE) file for more details.
