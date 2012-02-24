#include <iostream>
#include <list>
#include <string>
#include <vector>
#include <cmath>
#include  <cassert>
#include <cfloat>
#include <algorithm>
#include <complex>
#include <ctime>
#include "Array.h"

// There is some mismatch between fftw++ and Array about what namespace these are in, so this is a kludge.
// #ifndef __FractalGravity_libs_kludge__
// #define __FractalGravity_libs_kludge__
// namespace Array {
//   inline void free0(void *p) { ::free0(p); }
//   template<typename T> inline void newAlign(T *&v, size_t len, size_t align) { ::newAlign(v, len, align); }
//   template<typename T> inline void deleteAlign(T *v, size_t len) { ::deleteAlign(v, len); }
// }
// #endif

#include "fftw++.h"
using namespace std;
using namespace Array;
using namespace fftwpp;
