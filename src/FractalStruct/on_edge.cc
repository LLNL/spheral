#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  bool on_edge(vector <int>& pos,vector <int>& Box)
  {
    if(pos[0] == Box[0] || pos[0] == Box[1])
      return(pos[1] >= Box[2] && pos[1] <= Box[3] && pos[2] >= Box[4] && pos[2] <= Box[5]);
    if(pos[1] == Box[2] || pos[1] == Box[3])
      return(pos[0] >= Box[0] && pos[0] <= Box[1] && pos[2] >= Box[4] && pos[2] <= Box[5]);
    if(pos[2] == Box[4] || pos[2] == Box[5])
      return(pos[0] >= Box[0] && pos[0] <= Box[1] && pos[1] >= Box[2] && pos[1] <= Box[3]);
    return false;
  }
  bool on_edge(array <int,3>& pos,vector <int>& Box)
  {
    if(pos[0] == Box[0] || pos[0] == Box[1])
      return(pos[1] >= Box[2] && pos[1] <= Box[3] && pos[2] >= Box[4] && pos[2] <= Box[5]);
    if(pos[1] == Box[2] || pos[1] == Box[3])
      return(pos[0] >= Box[0] && pos[0] <= Box[1] && pos[2] >= Box[4] && pos[2] <= Box[5]);
    if(pos[2] == Box[4] || pos[2] == Box[5])
      return(pos[0] >= Box[0] && pos[0] <= Box[1] && pos[1] >= Box[2] && pos[1] <= Box[3]);
    return false;
  }
  bool on_edge(Point* p,vector <int>& Box)
  {
    array<int,3>pos(p->get_pos_point_a());
    if(pos[0] == Box[0] || pos[0] == Box[1])
      return(pos[1] >= Box[2] && pos[1] <= Box[3] && pos[2] >= Box[4] && pos[2] <= Box[5]);
    if(pos[1] == Box[2] || pos[1] == Box[3])
      return(pos[0] >= Box[0] && pos[0] <= Box[1] && pos[2] >= Box[4] && pos[2] <= Box[5]);
    if(pos[2] == Box[4] || pos[2] == Box[5])
      return(pos[0] >= Box[0] && pos[0] <= Box[1] && pos[1] >= Box[2] && pos[1] <= Box[3]);
    return false;
  }
}
