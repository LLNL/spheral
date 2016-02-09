#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void poisson_solver_struct(Fractal& fractal,Fractal_Memory& mem,const int& level)
  {
    int spacing=Misc::pow(2,fractal.get_level_max()-level);
    for(int ni=0;ni<2;ni++)
      {
	vector <vector <Point*> >hypre_points;
	vector <vector <Point*> >SPoints;
	vector <vector <int> >SBoxes;
	hypre_points_struct(mem,mem.all_groups[level],hypre_points,ni == 1,level);
	hypre_points_boxes(hypre_points,spacing,SBoxes,SPoints);
	hypre_points.clear();
	hypre_world_create(mem,level,SBoxes,ni == 1);
	hypre_solve_struct(mem,level,SBoxes,SPoints);
	//
	hypre_world_destroy();
	SBoxes.clear();
	SPoints.clear();
      }
  }
}
