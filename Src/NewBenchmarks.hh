#ifndef __NEWBENCHMARKS_HH_
#define __NEWBENCHMARKS_HH_

// Mostly identical to Benchmarks.hh

#include <string.h>
#include <gmpxx.h>
#include <fstream>
#include "helpers.hh"
#include "Matrix.hh"

// Runs CG on the given matrix and plots residual after each iteration.
void CGTest(Matrix<mpf_class> M, string matrixname="unknown matrix",
    string identifier="", bool quire=false, bool plot=false,
    double relativeTolerance=1e-5, bool scale=false, bool stochastic=false);

// Runs Trisolve on the given matrix and reports final residual
void trisolveTest(Matrix<mpf_class> M, string matrixname="unknkown matrix",
    string identifier="", bool cholesky=false, bool scale=false);

void fftTest();

void fftTest2D();

void convolveImage(string imgFilename, string kernelFilename, string output);

#endif