#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  extern int node_number;
  extern bool node_master;
  extern bool node_client;
  extern int nodes_total;
  extern int nodes_zero_level_total;
  extern int nodes_zero_level_1;
  extern bool node_zero_level;
  extern bool nodes_multi;
  void node_start(Fractal& fractal,Fractal_Memory& mem)
  {
    if(!FractalSpace::nodes_multi)
      return;
    int length=mem.grid_length;
    double points=mem.node_points_max;
    double padd=1+max(1,mem.padding);
    int n0=pow(points+0.01,1.0/3.0)-padd*2.0;
    while(length % n0 != 0)
      n0--;
    int n1=length/n0;
    FractalSpace::nodes_zero_level_1=n1;
    FractalSpace::nodes_multi=FractalSpace::nodes_total > 1;
    FractalSpace::nodes_zero_level_total=n1*n1*n1;
    if(FractalSpace::nodes_zero_level_total > FractalSpace::nodes_total)
      {
	node_end();
	assert(0);
      }
    FractalSpace::node_zero_level=node_client && FractalSpace::node_number <= FractalSpace::nodes_zero_level_total;
  }
}
