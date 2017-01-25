#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_clever_boxes(Fractal_Memory& mem,vector <vector<Point*> >& hypre_points,int spacing,
			  int VOLMIN,double FILLFACTOR,
			  vector < vector<int> >& SBoxes,vector < vector<Point*> >& SPoints)
  {
    int niMAX=3;
    for(int ni=0;ni<niMAX;ni++)
      {
	SBoxes.clear();
	SPoints.clear();
	hypre_points_boxes(mem,hypre_points,spacing,VOLMIN,FILLFACTOR,SBoxes,SPoints);
	// hypre_points=std::move(SPoints);
	hypre_points=SPoints;
      }
    hypre_points_zero(SPoints);
  }
}
