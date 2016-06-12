#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void clean_overlaps(int spacing,vector<bool>& STrouble,vector<vector<int>>& SBoxes,vector<vector<Point*>>& SPoints)
  {
    int RANK=-1;
    MPI_Comm_rank(MPI_COMM_WORLD,&RANK);
    bool RANKY=RANK==21;
    // if(RANKY)
    cerr << " ENTER DUPES " << RANK << " " << SBoxes.size() << " " << SPoints.size() << endl;
    vector<vector<int>>SBover;
    vector<vector<Point*>>SPover;
    int nBP=0;
    int Ngood=0;
    bool foundit=false;
    for(vector <int> &SB : SBoxes)
      {
	if(STrouble[nBP])
	  {
	    foundit=true;
	    SBover.resize(SBover.size()+1);
	    SPover.resize(SPover.size()+1);
	    SBover.back()=SB;
	    SPover.back().assign(SPoints[nBP].begin(),SPoints[nBP].end());
	  }
	else
	  {
	    if(foundit)
	      {
		SBoxes[Ngood]=SBoxes[nBP];
		SPoints[Ngood].assign(SPoints[nBP].begin(),SPoints[nBP].end());
	      }
	    Ngood++;
	  }
	nBP++;
      }
    assert(foundit);
    SBoxes.resize(nBP);
    SPoints.resize(nBP);
    std::map<array<int,4>,Point*,point_comp4> dupes;
    vector<int>pos(3);
    array<int,4>ar4;
    int nPa=0;
    for(vector<int> &SB : SBover)
      {
	int nPb=0;
	for(int nz=SB[4];nz<=SB[5];nz+=spacing)
	  {
	    for(int ny=SB[2];ny<=SB[3];ny+=spacing)
	      {
		for(int nx=SB[0];nx<=SB[1];nx+=spacing)
		  {
		    Point* pp=SPover[nPa][nPb];
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
		    if(!ret.second && pp != 0)
		      {
			dupes.erase(ar4);
			ret=dupes.insert(std::pair<array<int,4>,Point*>(ar4,pp));
			assert(ret.second);
		      }
		    nPb++;  
		  }
	      }
	  }
	nPa++;
      }
    vector<vector<Point*>>hypre_points(SBover.size());
    SBover.clear();
    SPover.clear();
    for(auto &d : dupes)
      hypre_points[d.first[3]].push_back(d.second);
    int counta=0;
    int countb=0;
    for(auto &hp : hypre_points)
      {
	if(countb > counta && hp.size() > 0)
	  hypre_points[counta].assign(hp.begin(),hp.end());
	if(hp.size() > 0)
	  counta++;
	countb++;
      }
    hypre_points.resize(counta);
    hypre_points_boxes(hypre_points,spacing,false,SBoxes,SPoints);
    return;
  }
}
