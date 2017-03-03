#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void any_overlaps(Fractal_Memory& mem,int spacing,int VOLMIN,double FILLFACTOR,vector<vector<int>>& SBoxes,vector<vector<Point*>>& SPoints)
  {
    bool overlaps=false;
    std::map<array<int,4>,Point*,point_comp4> dupes;
    vector<int>pos(3);
    array<int,4>ar4;
    vector<bool>STrouble(SBoxes.size(),false);
    int nPa=0;
    for(vector<int> &SB : SBoxes)
      {
	int nPb=0;
	for(int nz=SB[4];nz<=SB[5];nz+=spacing)
	  {
	    for(int ny=SB[2];ny<=SB[3];ny+=spacing)
	      {
		for(int nx=SB[0];nx<=SB[1];nx+=spacing)
		  {
		    Point* pp=SPoints[nPa][nPb];
		    if(pp == 0)
		      ar4={{nx,ny,nz,nPa}};
		    else
		      {
			pp->get_pos_point(pos);
			assert(pos[0] == nx);
			assert(pos[1] == ny);
			assert(pos[2] == nz);
			std::move(pos.begin(),pos.end(),ar4.begin());
			ar4[3]=nPa;
		      }
		    std::pair<std::map<array<int,4>,Point*>::iterator,bool> ret;
		    ret=dupes.insert(std::pair<array<int,4>,Point*>(ar4,pp));
		    if(!ret.second)
		      {
			overlaps=true;
			array<int,4>arprevious=ret.first->first;
			STrouble[arprevious[3]]=true;
			STrouble[nPa]=true;
		      }
		    nPb++;  
		  }
	      }
	  }
	nPa++;
      }
    if(overlaps)
      clean_overlaps(mem,spacing,VOLMIN,FILLFACTOR,STrouble,SBoxes,SPoints);
  }
}
