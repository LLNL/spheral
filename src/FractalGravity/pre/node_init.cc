#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  int node_number=0;
  int nodes_total=1;
  int node_points=0;
  int node_particles=0;
  int nodes_zero_level_1;
  int nodes_zero_level_total=0;
  bool node_master=false;
  bool node_client=true;
  bool node_zero_level=false;
  bool nodes_multi=false;
    //
  void node_init()
  {
    mpi_init();
    int spam=0;
    mpi_node_number(spam);
    FractalSpace::node_number=spam;
    FractalSpace::node_master=FractalSpace::node_number==0;
    FractalSpace::node_client=FractalSpace::node_number!=0;
    mpi_nodes_total(spam);
    FractalSpace::nodes_total=spam;
    FractalSpace::nodes_multi=spam > 1;
    cout << " not yet " << endl;
  }
}
