#include <math.h>
#include <iostream>
#include <bitset>
#include <vector>
#include "FFT.hh"
// #include "gmpxx.h"
using namespace std;

/*Convenience Functions*/
int mod(int m, int n)  { return m % n >= 0 ? (m % n) : ((m % n) + n); }

int powerOfTwo(uint n) {return (n & (n-1)) == 0;}

// *****************************

template <class T>
void partition(vector<T> &out, vector<T> &in, int odd) {
    for (int i=0;i<out.size();++i) out[i]=in[2*i+odd];
}

