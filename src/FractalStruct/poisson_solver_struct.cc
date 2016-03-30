#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void poisson_solver_struct(Fractal& fractal,Fractal_Memory& mem,const int& level)
  {
    static int _COUNTER=0;
    int RANK=-1;
    MPI_Comm_rank(MPI_COMM_WORLD,&RANK);
    int spacing=Misc::pow(2,fractal.get_level_max()-level);
    cerr << "POISSON A " << RANK << " " << level << " " << " " << _COUNTER << "\n";
    for(int ni=0;ni<2;ni++)
      {
	cerr << "POISSON B " << RANK << " " << level << " " << ni << " " << _COUNTER << "\n";
	vector <vector <Point*> >hypre_points;
	vector <vector <Point*> >SPoints;
	vector <vector <int> >SBoxes;
	cerr << " HYPRE RES B " <<  RANK << " " << level << " " << ni << " " << _COUNTER << "\n";
	hypre_points_struct(mem,mem.all_groups[level],hypre_points,ni == 1,level);
	cerr << " HYPRE RES C " <<  RANK << " " << level << " " << ni << " " << _COUNTER << "\n";
	hypre_points_boxes(hypre_points,spacing,SBoxes,SPoints);
	cerr << " HYPRE RES D " <<  RANK << " " << level << " " << ni << " " << _COUNTER << "\n";
	hypre_points.clear();
	hypre_world_create(mem,level,SBoxes,ni == 1);
	cerr << " HYPRE RES E " <<  RANK << " " << level << " " << ni << " " << _COUNTER << "\n";
	if(mem.p_mess->IAmAHypreNode)
	  {
	    hypre_solve_struct(mem,level,SBoxes,SPoints);
	    if(ni == 100)
	      add_buffer_values(mem,level,SBoxes,SPoints);
	    cerr << " HYPRE RES F " <<  RANK << " " << level << " " << ni << " " << _COUNTER << "\n";
	  }
	cerr << " HYPRE RES G " <<  RANK << " " << level << " " << ni << " " << _COUNTER << "\n";
	mem.p_mess->HypreGroupFree();
	cerr << " HYPRE RES H " <<  RANK << " " << level << " " << ni << " " << _COUNTER << "\n";
	SBoxes.clear();
	SPoints.clear();
	_COUNTER++;
      }
    cerr << " HYPRE RES I " <<  RANK << " " << level << " " << " " << _COUNTER << "\n";
  }
}
