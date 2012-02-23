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

// These shouldn't be necessary, but there's some inconsistency between Array and fftw++.
namespace Array {
  void free0(void *p) { ::free0(p); }
  template<class T> void newAlign(T *&v, size_t len, size_t align) { ::newAlign(v, len, align); }
  template<class T> void deleteAlign(T *v, size_t len) { ::deleteAlign(v, len); }
}

#include "fftw++.h"
using namespace std;
using namespace Array;
using namespace fftwpp;
