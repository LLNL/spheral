#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_points_zero(vector<vector<Point*>>& SPoints)
  {
    Point* pFAKE=0;
    for(auto &SP : SPoints)
      for(int S=0;S<SP.size();S++)
	{
	  Point* p=SP[S];
	  if(p != 0 && p->get_really_passive())
	    {
	      delete p;
	      SP[S]=pFAKE;
	    }
	}
  }
  // void hypre_points_unzero(vector<vector<int>>& SBoxes,vector<vector<Point*>>& SPoints)
  // {
  //   int npa=0;
  //   for(vector<int> &SB : SBoxes)
  //     {
  // 	int npb=0;
  // 	for(int nz=SB[4];nz<=SB[5];nz+=spacing)
  // 	  {
  // 	    for(int ny=SB[2];ny<=SB[3];ny+=spacing)
  // 	      {
  // 		for(int nx=SB[0];nx<=SB[1];nx+=spacing)
  // 		  {
  // 		    if(SPoints[npa][npb]==0)
  // 		      {
  // 			Point* pp=new Point;
  // 			pp->set_really_passive(true);
  // 			pp->set_pos_point(nx,ny,nz);
  // 			SPoints[npa][npb]=pp;
  // 		      }
  // 		    npb++;
  // 		  }
  // 	      }
  // 	  }
  // 	npa++;
  //     }
}
