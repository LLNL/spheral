#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  template <class GO_AWAY> void really_clear(vector <GO_AWAY>& die)
  {
    die.clear();
  }
}
// namespace FractalSpace
// {
//   template <class GO_AWAY> void really_clear(vector <GO_AWAY>& die)
//   {
//     die.clear();
//     die.shrink_to_fit();
//   }
// }
// namespace FractalSpace
// {
//   template <class GO_AWAY> void really_resize(vector <GO_AWAY>& die,int how_big)
//   {
//     die.resize(how_big);
//     die.shrink_to_fit();
//   }
// }
// namespace FractalSpace
// {
//   template <class GO_AWAY> void really_resize(vector <vector <GO_AWAY> >& die,int how_big)
//   {
//     die.resize(how_big);
//     die.shrink_to_fit();
//     for(vector <GO_AWAY>& di : die)
//       di.shrink_to_fit();
//   }
// }
// namespace FractalSpace
// {
//   template <class GO_AWAY> void really_resize(vector <vector <vector <GO_AWAY> > >& die,int how_big)
//   {
//     die.resize(how_big);
//     die.shrink_to_fit();
//     for(vector <vector <GO_AWAY> >& di : die)
//       {
// 	di.shrink_to_fit();
// 	for(vector <GO_AWAY>& d : di)
// 	  d.shrink_to_fit();
//       }
//   }
// }
// namespace FractalSpace
// {
//   template <class GO_AWAY> void multi_shrink(vector <vector <GO_AWAY> >& die)
//   {
//     die.shrink_to_fit();
//     for(vector <GO_AWAY>& di : die)
//       di.shrink_to_fit();
//   }
// }
// namespace FractalSpace
// {
//   template <class GO_AWAY> void multi_shrink(vector <vector <vector <GO_AWAY> > >& die)
//   {
//     die.shrink_to_fit();
//     for(vector <vector <GO_AWAY> >& di : die)
//       {
// 	di.shrink_to_fit();
// 	for(vector <GO_AWAY>& d : di)
// 	  d.shrink_to_fit();
//       }
//   }
// }
