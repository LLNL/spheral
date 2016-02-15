#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void poisson_solver_struct(Fractal& fractal,Fractal_Memory& mem,const int& level)
  {
    int RANK=-1;
    MPI_Comm_rank(MPI_COMM_WORLD,&RANK);
//     bool Ranky = RANK == 21;
//     Ranky=false;
    int spacing=Misc::pow(2,fractal.get_level_max()-level);
    for(int ni=0;ni<2;ni++)
      {
	vector <vector <Point*> >hypre_points;
	vector <vector <Point*> >SPoints;
	vector <vector <int> >SBoxes;
	mem.p_mess->Full_Stop_Do_Not_Argue();
	cerr << " structAA " << RANK << " " << ni << " " << level  << " " << ni << " " << spacing << endl;
	hypre_points_struct(mem,mem.all_groups[level],hypre_points,ni == 1,level);
	mem.p_mess->Full_Stop_Do_Not_Argue();
	cerr << " structBA " << RANK << " " << ni << " " << level  << " " << ni << " " << spacing << endl;
	hypre_points_boxes(hypre_points,spacing,SBoxes,SPoints);
	mem.p_mess->Full_Stop_Do_Not_Argue();
	cerr << " boxA " << RANK << " " << ni << " " << level << " " << SBoxes.size() << " " << ni << " " << spacing << endl;
	hypre_points.clear();
	cerr << " createAA " << RANK << " " << ni << " " << level  << " " << ni << " " << spacing << endl;
	hypre_world_create(mem,level,SBoxes,ni == 1);
	mem.p_mess->Full_Stop_Do_Not_Argue();
	cerr << " createBA " << RANK << " " << ni << " " << level  << " " << ni << " " << spacing << endl;
	if(mem.p_mess->IAmAHypreNode)
	  hypre_solve_struct(mem,level,SBoxes,SPoints);
	mem.p_mess->Full_Stop_Do_Not_Argue();
	cerr << " runaA " << RANK << " " << ni << " " << level << " " << ni << " " << spacing << endl;
	mem.p_mess->HypreGroupFree();
	cerr << " freeA " << RANK << " " << ni << " " << level  << " " << ni << " " << spacing << endl;
	mem.p_mess->Full_Stop_Do_Not_Argue();
	SBoxes.clear();
	SPoints.clear();
	assert(false);
      }
  }
}
// 	mem.p_mess->Full_Stop_Do_Not_Argue();
