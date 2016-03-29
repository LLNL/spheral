#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void sor_solver(Group& group, Fractal& fractal)
  {
    const double pi = 4.0*atan(1.0);
    double g_c=4.0*pi/(double)(fractal.get_grid_length()*fractal.get_grid_length())*pow(4.0,-group.get_level());
    group.set_force_const(g_c);
    //
    int n_tot=group.list_points.size();
    int n_inside=0;
    vector < vector <Point*> > list_edge(6);
    vector <Point*> ud(6);
    int n_xy=0;
    int n_xz=0;
    int n_yz=0;
    //
    for(vector <Point*>::iterator point_itr=group.list_points.begin();point_itr != group.list_points.end();++point_itr)
      {
	Point& point=**point_itr;
	point.get_point_ud(ud);
	if(!point.get_inside())
	  {
	    for(int ni=0;ni <6;ni++)
	      if(ud[ni] != 0 && ud[ni]->get_inside()) list_edge[ni].push_back(ud[ni]);
	  }
	else
	  n_inside++;
	if(ud[0]==0) {n_yz++;}
	if(ud[2]==0) {n_xz++;}
	if(ud[4]==0) {n_xy++;}
      }
    double rj=(cos(pi*(double)n_yz/(double)n_tot)+
	       cos(pi*(double)n_xz/(double)n_tot)+
	       cos(pi*(double)n_xy/(double)n_tot))/3.0;
    group.set_rjac(rj);
    //
    int dir=5;
    if(n_yz >= max(n_xz,n_xy))
      dir=1;
    else if(n_xz >= max(n_yz,n_xy))
      dir=3;
    //
    dir=1;
    //
    sor(group,fractal,list_edge[dir],dir);
  }
}
