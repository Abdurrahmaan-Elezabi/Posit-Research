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
#include "bfloat16.hh"

using namespace std;
using namespace cimg_library;
using half_float::half;

// TODO: half doesn't work
void CGTest(Matrix<mpf_class> M, string matrixname="unknown matrix",
    string identifier="", bool quire=false, bool plot=false,
    double relativeTolerance=1e-5, bool scale=false, bool stochastic=false) {
    if (!(M.isSymmetric())) {fprintf(stderr, "Please input symmetric matrix for CG test."); return;}

    string infoFilename = "plots/CG" + identifier;
    string plotfile;

    if (plot) plotfile = "plots/" + (matrixname + identifier) + ".csv";
    ofstream infoFile;
    
    int n = M.nCols();

    vector<mpf_class> xM = vector<mpf_class>(n, 1/sqrt(mpf_class(n)));
    vector<mpf_class> bM = matVec(M, xM);

    if (scale) {
        M.scaleNorm(M, bM);
    }

    Matrix<double    > D;
    Matrix<float     > F;
    //Matrix<Posit32gmp> P;
    //Matrix<half      > H;
    //Matrix<bfloat16  > B;

    D.set(M);
    F.set(M);
    //P.set(M);
    //H.set(M);
    //B.set(M);

    vector<double    > xD(n);
    vector<float     > xF(n);
    //vector<Posit32gmp> xP(n);
    //vector<half      > xH(n);
    //vector<bfloat16  > xB(n);

    vector<double    > bD(n);
    vector<float     > bF(n);
    //vector<Posit32gmp> bP(n);
    //vector<half      > bH(n);
    //vector<bfloat16  > bB(n);

    downcast(bD, bM);
    downcast(bF, bM);
    //downcast(bP, bM);
    //downcast(bH, bM);
    //downcast(bB, bM);

    double tolerance = M.vectorNorm(bM).get_d()*relativeTolerance;
    cout << "Aiming for relative error of " << relativeTolerance << endl;

    Posit32::clearCounter();
    int f, d, p, h, b;
    puts("starting CG step");
    
    if (quire) {
        f = F.conjugateGradientSolverQuire(tolerance, F, bF, xF, plotfile);
        d = D.conjugateGradientSolverQuire(tolerance, D, bD, xD, plotfile);
        //p = conjugateGradientSolverQ(tolerance, P, bP, xP, plotfile, "", false);
        //b = B.conjugateGradientSolverQuire(tolerance, B, bB, xB, plotfile);
        //h = H.conjugateGradientSolverQuire(tolerance, H, bH, xH, plotfile);
    } else if (stochastic) {
        f = F.conjugateGradientSolverStochastic(tolerance, F, bF, xF, plotfile);
        d = D.conjugateGradientSolverStochastic(tolerance, D, bD, xD, plotfile);
        //p = P.conjugateGradientSolverStochastic(tolerance, P, bP, xP, plotfile);
        //b = B.conjugateGradientSolverStochastic(tolerance, B, bB, xB, plotfile);
        //h = H.conjugateGradientSolverStochastic(tolerance, H, bH, xH, plotfile);
    } else {
        f = F.conjugateGradientSolver(tolerance, F, bF, xF, plotfile);
        d = D.conjugateGradientSolver(tolerance, D, bD, xD, plotfile);
        //p = P.conjugateGradientSolver(tolerance, P, bP, xP, plotfile);
        //b = B.conjugateGradientSolver(tolerance, B, bB, xB, plotfile);
        //h = H.conjugateGradientSolver(tolerance, H, bH, xH, plotfile);
    }
    

    cout << "Scale? " << scale << endl;
    cout << "Max entry: " << M.getMax() << endl;
    cout << "Min entry: " << M.getMin() << endl;
    // Using minus here made more sense to me when talking about range.
    // I was also getting issues with division by zero. Should I change it back?
    cout << "Range: " << M.getMax() - M.getMin() << endl;
    cout << "double iterations = "   << d << endl;
    cout << "float iterations = "    << f << endl;
    //cout << "posit iterations = "    << p << endl;
    //cout << "half iterations = "     << h << endl;
    //cout << "bfloat16 iterations = " << b << endl;

    infoFile.open(infoFilename, ofstream::app);
    infoFile << matrixname << endl;
    infoFile << "relative error tolerance: " << relativeTolerance << endl;
    infoFile << "Scale?: " << scale << endl;
    infoFile << "Max/min entry: " << M.getMax() << " : " << M.getMin() << "Range: " << M.getMax() - M.getMin() << endl;
    infoFile << "double iterations = "   << d << endl;
    infoFile << "float iterations = "    << f << endl;
    //infoFile << "posit iterations = "    << p << endl;
    //infoFile << "half iterations = "     << h << endl;
    //infoFile << "bfloat16 iterations = " << b << endl;
    infoFile.close();

    vector<mpf_class> bmD(n), bmF(n), bmP(n), bmH(n), bmB(n);

    upcast(bmD,matVec(D, xD)); 
    upcast(bmF,matVec(F, xF));
    //upcast(bmP,matVec(P, xP));
    //upcast(bmH,matVec(H, xH));
    //upcast(bmB,matVec(B, xB));

    mpf_class rD = M.vectorNorm(M.vectorCombination(1.0, bM, -1.0, bmD));
    mpf_class rF = M.vectorNorm(M.vectorCombination(1.0, bM, -1.0, bmF));
    //mpf_class rP = M.vectorNorm(M.vectorCombination(1.0, bM, -1.0, bmP));
    //mpf_class rH = M.vectorNorm(M.vectorCombination(1.0, bM, -1.0, bmH));
    //mpf_class rB = M.vectorNorm(M.vectorCombination(1.0, bM, -1.0, bmB));

    puts("\nFinal Residuals:");
    cout << '\t' << "double relative residual = "   << rD/M.vectorNorm(bM) << endl;
    cout << '\t' << "float relative residual = "    << rF/M.vectorNorm(bM) << endl;
    //cout << '\t' << "posit relative residual = "    << rP/M.vectorNorm(bM) << endl;
    //cout << '\t' << "half relative residual = "     << rH/M.vectorNorm(bM) << endl;
    //cout << '\t' << "bfloat16 relative residual = " << rB/M.vectorNorm(bM) << endl;

    infoFile.open(infoFilename, ofstream::app);
    infoFile << endl;
    infoFile << '\t' << "double relative residual = "   << rD/M.vectorNorm(bM) << endl;
    infoFile << '\t' << "float relative residual = "    << rF/M.vectorNorm(bM) << endl;
    //infoFile << '\t' << "posit relative residual = "    << rP/M.vectorNorm(bM) << endl;
    //infoFile << '\t' << "half relative residual = "     << rH/M.vectorNorm(bM) << endl;
    //infoFile << '\t' << "bfloat16 relative residual = " << rB/M.vectorNorm(bM) << endl;
    infoFile << endl;
    infoFile.close();
}

// TODO: implement half after fixing issues
void trisolveTest(Matrix<mpf_class> M, string matrixname="unknown matrix",
    string identifier="", bool cholesky=false, bool scale=false) {
    if (!(M.isSquare())) {fprintf(stderr, "Please input square matrix for Tri-Solve test."); return;}
    cout << "Running direct solve benchmark..." << endl;
    int n = M.nCols();

    string infoFilename = "plots/TriSolve" + identifier;
    ofstream infoFile;

    Matrix<double    > D;
    Matrix<float     > F;
    Matrix<Posit32gmp> P;
    Matrix<bfloat16  > B;

    vector<mpf_class> xM = vector<mpf_class>(n, 1/sqrt(mpf_class(n)));
    vector<mpf_class> bM = matVec(M, xM);

    D.set(M);
    F.set(M);
    P.set(M);
    B.set(M);

    vector<double    > xD(n);
    vector<float     > xF(n);
    vector<Posit32gmp> xP(n);
    vector<bfloat16  > xB(n);

    vector<double> bD(n);
    vector<float> bF(n);
    vector<Posit32gmp> bP(n);
    vector<bfloat16> bB(n);

    downcast(bD, bM);
    downcast(bF, bM);
    downcast(bP, bM);
    downcast(bB, bM);

    if (scale) {
        M.diagScaleAvg(3, M, bM);
        D.diagScaleAvg(3, D, bD);
        F.diagScaleAvg(3, F, bF);
        P.diagScaleAvg(3, P, bP);
        B.diagScaleAvg(3, B, bB);
    }

    Posit32::clearCounter();

    if (cholesky) {
        D.symmetricTriSolve(xD, bD);
        F.symmetricTriSolve(xF, bF);
        P.symmetricTriSolve(xP, bP);
        B.symmetricTriSolve(xB, bB);
    } else {
        D.triSolve(xD, bD);
        F.triSolve(xF, bF);
        P.triSolve(xP, bP);
        B.triSolve(xB, bB);
    }

    infoFile.open(infoFilename, ofstream::app);
    infoFile << matrixname << endl;
    infoFile << "Scale?: " << scale << endl;
    infoFile << "Cholesky?: " << cholesky << endl;
    // infoFile << "Average posit advantage: " << Posit32::distillAdvantage() << endl; // has a bug
    infoFile.close();

    vector<mpf_class> bmD(n), bmF(n), bmP(n), bmB(n);
    upcast(bmD,matVec(D, xD));
    upcast(bmF,matVec(F, xF));
    upcast(bmP,matVec(P, xP));
    upcast(bmB,matVec(B, xB));

    mpf_class rD = M.vectorNorm(M.vectorCombination(1.0, bM, -1.0, bmD));
    mpf_class rF = M.vectorNorm(M.vectorCombination(1.0, bM, -1.0, bmF));
    mpf_class rP = M.vectorNorm(M.vectorCombination(1.0, bM, -1.0, bmP));
    mpf_class rB = M.vectorNorm(M.vectorCombination(1.0, bM, -1.0, bmB));

    cout << "Final residuals:" << endl;
    cout << "Scale?: " << scale << endl;
    cout << "Cholesky?: " << cholesky << endl;
    cout << "double relative residual = "   << rD/M.vectorNorm(bM) << endl;
    cout << "float relative residual = "    << rF/M.vectorNorm(bM) << endl;
    cout << "posit relative residual = "    << rP/M.vectorNorm(bM) << endl;
    cout << "bfloat16 relative residual = " << rB/M.vectorNorm(bM) << endl;

    infoFile.open(infoFilename, ofstream::app);
    infoFile << '\t' << "double relative residual = "   << rD/M.vectorNorm(bM) << endl;
    infoFile << '\t' << "float relative residual = "    << rF/M.vectorNorm(bM) << endl;
    infoFile << '\t' << "posit relative residual = "    << rP/M.vectorNorm(bM) << endl;
    infoFile << '\t' << "bfloat16 relative residual = " << rB/M.vectorNorm(bM) << endl;
    infoFile.close();
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
    bool scale;
    bool plot;
    bool quire;
    bool stochastic;
    string identifier = "";
    string matrixname;
    double tolerance = 1e-7;

    cout << "Enter mtx file name (see MatrixMarket):" << endl;
    cin >> filename;
    while (!systemM.loadMPF((root + filename).c_str())) {
        cout << "Try again" << endl;
        cin >> filename;
    }
    cout << "Enter matrix name:" << endl;
    cin >> matrixname;
    cout << "Enter identifier:" << endl;
    cin >> identifier;
    cout << "Scale? (y/n)" << endl;
    scale = yesNo();
    cout << "Quire? (y/n)" << endl;
    quire = yesNo();
    cout << "Stochastic? (y/n)" << endl;
    stochastic = yesNo();

    switch (testID) {
        case 0:
            //cout << "Plot? (y/n)" << endl;
            //plot = yesNo();
            // No reason not to plot?
            plot = true;
            CGTest(systemM, matrixname, identifier, quire, plot, tolerance, scale, stochastic);
            break;
        case 1:
            cout << "Cholesky? (y/n)" << endl;
            bool cholesky = yesNo();
            trisolveTest(systemM, matrixname, identifier, cholesky, scale);
            break;
    }
    return EXIT_SUCCESS;
}