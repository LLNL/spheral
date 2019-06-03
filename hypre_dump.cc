#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_dump(int level,vector <Point*>& hypre_points,ofstream& FH)
  {
    FH << " HYPRE DUMP " << level << "\n";
    for(auto p : hypre_points)
      {
	FH << " EDUMP ";
	p->dumpp(FH);
      }
  }
}
