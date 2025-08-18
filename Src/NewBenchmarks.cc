#include <vector>
#include <cmath>
#include <stdio.h>
#include <iostream>
#include <CImg.h>
#include <fenv.h>
#include <dirent.h>
#include <string.h>
#include <boost/filesystem.hpp>
#include "FFT.hh"
#include "Quire.hh"

using namespace std;
using namespace cimg_library;
using half_float::half;

// For now, only works with 32 bit floats.
// TODO: add options for other precisions
void CGTest(Matrix<mpf_class> M, string matrixname="unknown matrix", double relativeTolerance=1e-5) {
    if (!(M.isSymmetric())) {fprintf(stderr, "Please input symmetric matrix for CG test."); return;}

    int n = M.nCols();

    vector<mpf_class> xM = vector<mpf_class>(n, 1/sqrt(mpf_class(n)));
    vector<mpf_class> bM = matVec(M, xM);

    Matrix<double    > D;
    Matrix<float     > F;

    D.set(M);
    F.set(M);

    vector<double    > xD(n);
    vector<float     > xF(n);

    vector<double    > bD(n);
    vector<float     > bF(n);

    downcast(bD, bM);
    downcast(bF, bM);

    double tolerance = M.vectorNorm(bM).get_d()*relativeTolerance;
    cout << "Aiming for relative error of " << relativeTolerance << endl;

    int f, d;
    puts("starting CG step");
    f = F.conjugateGradientSolver(tolerance, F, bF, xF);
    d = D.conjugateGradientSolver(tolerance, D, bD, xD);

    cout << '\t' << "double iterations = "  << d << endl;
    cout << '\t' << "float iterations = "   << f << endl;

    vector<mpf_class> bmD(n), bmF(n);
    upcast(bmD,matVec(D, xD)); 
    upcast(bmF,matVec(F, xF));
    mpf_class rD = M.vectorNorm(M.vectorCombination(1.0, bM, -1.0, bmD));
    mpf_class rF = M.vectorNorm(M.vectorCombination(1.0, bM, -1.0, bmF));

    puts("\nFinal Residuals:");
    cout << '\t' << "double relative residual = "   << rD/M.vectorNorm(bM) << endl;
    cout << '\t' << "float relative residual = "    << rF/M.vectorNorm(bM) << endl;
}

void trisolveTest(Matrix<mpf_class> M) {
    if (!(M.isSquare())) {fprintf(stderr, "Please input square matrix for Tri-Solve test."); return;}
    cout << "Running direct solve benchmark..." << endl;
    int n = M.nCols();

    Matrix<double    > D;
    Matrix<float     > F;

    vector<mpf_class> xM = vector<mpf_class>(n, 1/sqrt(mpf_class(n)));
    vector<mpf_class> bM = matVec(M, xM);

    D.set(M);
    F.set(M);

    vector<double    > xD(n);
    vector<float     > xF(n);

    vector<double> bD(n);
    vector<float> bF(n);

    downcast(bD, bM);
    downcast(bF, bM);

    D.triSolve(xD, bD);
    F.triSolve(xF, bF);

    vector<mpf_class> bmD(n), bmF(n);
    upcast(bmD,matVec(D, xD));
    upcast(bmF,matVec(F, xF));

    mpf_class rD = M.vectorNorm(M.vectorCombination(1.0, bM, -1.0, bmD));
    mpf_class rF = M.vectorNorm(M.vectorCombination(1.0, bM, -1.0, bmF));

    cout << "Final residuals:" << endl;
    cout << "double relative residual = "   << rD/M.vectorNorm(bM) << endl;
    cout << "float relative residual = "    << rF/M.vectorNorm(bM) << endl;
}

int main(int argc, char* argv[]) {
    if (argc != 1) {
        cout << "Usage: ./NewBenchmarks" << endl;
        return EXIT_FAILURE;
    }
    cout << "Enter test ID (0 for CG, 1 for Trisolve):" << endl;
    int testID = getInteger();
    while (testID < 0 || testID > 1) {
        cout << "Please enter a valid test ID" << endl;
        testID = getInteger();
    }
    Matrix<mpf_class> systemM;
    string root = "./MatrixMarket/";
    string filename;
    cout << "Enter mtx file name (see MatrixMarket):" << endl;
    cin >> filename;
    while (!systemM.loadMPF((root + filename).c_str())) {
        cout << "Try again" << endl;
        cin >> filename;
    }
    switch (testID) {
        case 0:
            CGTest(systemM);
            break;
        case 1:
            trisolveTest(systemM);
            break;
    }
    return EXIT_SUCCESS;
}