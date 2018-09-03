#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <vector>
#include <cmath>
#include <cassert>
#include <cfloat>
#include <numeric>
#include <algorithm>
#include <complex>
#include <ctime>
#include <random>
#include <climits>
#include <cerrno>
#include <sys/stat.h>
int main()
{
  vector <int>a(3);
  a[0]=1;
  a[1]=2;
  a[2]=3;
  vector <int>b;
  cout << " a, b 1 " << a.begin() << " " << b.begin() << endl;
  a.swap(b);
  cout << " a, b 2 " << a.begin() << " " << b.begin() << endl;
}
