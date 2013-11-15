#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  template <class T> int shortest_vector(vector<T>& veca,vector<T>& vecb,vector<T>& vecc)
  {
    int sizea=veca.size();
    int sizeb=vecb.size();
    int sizec=vecc.size();
    int sizeabc=min(min(sizea,sizeb),sizec);
    if(sizea == sizeabc) return 0;
    if(sizeb == sizeabc) return 1;
    if(sizec == sizeabc) return 2;
    return -1;
  }
  template int shortest_vector(vector<Point*>& veca,vector<Point*>& vecb,vector<Point*>& vecc);
}
