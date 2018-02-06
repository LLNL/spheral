#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  template <typename T> void clean_vector(vector<T>& vec)
  {
    vector<T> Brian_of_Nazareth;
    vec.swap(Brian_of_Nazareth);
  }
}
namespace FractalSpace
{
  template void clean_vector(vector <bool>& vec);
  template void clean_vector(vector <double>& vec);
  template void clean_vector(vector <int>& vec);
  template void clean_vector(vector <Point*>& vec);
  template void clean_vector(vector <Particle*>& vec);
}
namespace FractalSpace
{
  template <typename T> void clean_deque(deque<T>& deq)
  {
    deque<T> Brian_of_Nazareth;
    deq.swap(Brian_of_Nazareth);
  }
}
namespace FractalSpace
{
  template void clean_deque(deque <bool>& deq);
  template void clean_deque(deque <double>& deq);
  template void clean_deque(deque <int>& deq);
  template void clean_deque(deque <Point*>& deq);
  template void clean_deque(deque <Particle*>& deq);
}
