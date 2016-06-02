#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void remove_dupe_points(int spacing,vector<vector<Point*>>& hypre_points,vector<vector<int>>& SBoxes,vector<vector<Point*>>& SPoints)
  {
    bool no_overlaps=true;
    std::map<array<int,3>,Point*,point_comp2> dupes;
    vector<int>pos(3);
    array<int,3>ar3;
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
		      ar3={{nx,ny,nz}};
		    else
		      {
			pp->get_pos_point(pos);
			assert(pos[0] == nx);
			assert(pos[1] == ny);
			assert(pos[2] == nz);
			std::move(pos.begin(),pos.end(),ar3.begin());
		      }
		    std::pair<std::map<array<int,3>,Point*>::iterator,bool> ret;
		    ret=dupes.insert(std::pair<array<int,3>,Point*>(ar3,pp));
		    if(!ret.second && pp != 0)
		      {
			no_overlaps=false;
			dupes.erase(ar3);
			ret=dupes.insert(std::pair<array<int,3>,Point*>(ar3,pp));
		      }
		    nPb++;  
		  }
	      }
	  }
	nPa++;
      }
    hypre_points.clear();
    hypre_points.resize(SBoxes.size());
    nPa=0;
    for(vector<int> &SB : SBoxes)
      {
	for(int nz=SB[4];nz<=SB[5];nz+=spacing)
	  {
	    for(int ny=SB[2];ny<=SB[3];ny+=spacing)
	      {
		for(int nx=SB[0];nx<=SB[1];nx+=spacing)
		  {
		    ar3={{nx,ny,nz}};
		    std::map<array<int,3>,Point*,point_comp2>::iterator it=dupes.find(ar3);
		    assert(it != dupes.end());
		    hypre_points[nPa].push_back(it->second);
		  }
	      }
	  }
	nPa++;
      }
    if(no_overlaps)
      return;
    hypre_points_boxes(hypre_points,spacing,false,SBoxes,SPoints);
  }
}
