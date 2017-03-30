#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  template <typename T> void clean_vector(T& vec)
  {
    vector<T> Brian_of_Nazareth;
    vec.swap(Brian_of_Nazareth);
  }
}
