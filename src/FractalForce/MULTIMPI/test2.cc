#include <cstdlib>
#include <iostream>
#include <fstream>
using namespace std;
int main()
{
  int i=1;
  double x=1.0;
  for(int j=0;j<100;j++)
    {
      cout << " " << j << " " << i << " " << x << endl;
      i*=2;
      x*=2.0;
    }
  return 1;
}
