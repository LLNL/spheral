#include <cstdlib>
#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <vector>
#include <cmath>
#include <cassert>
#include <cfloat>
#include <algorithm>
#include <complex>
#include <ctime>
#include <climits>
using namespace std;
int fractal_force_wrapper(string* pfm,bool* pf);
int main()
{
  string* p_fractal_memory= 0;
  bool* p_fractal=0;
  cout << " pfm aa " << p_fractal_memory << endl;
  cout << " pf aa " << p_fractal << endl;
  int result=fractal_force_wrapper(p_fractal_memory,p_fractal);
  cout << " pfm a " << p_fractal_memory << endl;
  cout << " pf a " << p_fractal << endl;
  return result;
}
int fractal_force_wrapper(string* pfm,bool* pf)
{
  using namespace std;
  cout << " pfm b " << pfm << endl;
  cout << " pf b " << pf << endl;
  pfm=new string;
  pf=new bool;
  cout << " pfm c " << pfm << endl;
  cout << " pf c " << pf << endl;
}
