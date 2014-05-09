#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_dump(int level,vector <Point*>& hypre_points,ofstream& FH)
  {
    FH << " HYPRE DUMP " << level << "\n";
    for(vector<Point*>::const_iterator point_itr=hypre_points.begin();point_itr !=hypre_points.end();++point_itr)
      {
	FH << " EDUMP ";
	(*point_itr)->dumpp(FH);
      }
  }
}
