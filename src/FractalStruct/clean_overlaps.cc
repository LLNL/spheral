#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void clean_overlaps(Fractal_Memory& mem,int spacing,int VOLMIN,double FILLFACTOR,vector<bool>& STrouble,vector<vector<int>>& SBoxes,vector<vector<Point*>>& SPoints)
  {
    ofstream& FHT=mem.p_file->DUMPS;
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
		SBoxes[Ngood]=SB;
		SPoints[Ngood].assign(SPoints[nBP].begin(),SPoints[nBP].end());
	      }
	    Ngood++;
	  }
	nBP++;
      }
    assert(foundit);
    const int NB=SBoxes.size();
    SBoxes.resize(Ngood);
    SPoints.resize(Ngood);
    // int sp2=log((double)(spacing)+0.1)/log(2.0);
    // int level=mem.p_fractal->get_level_max()-sp2;
    FHT << " BOX PARAMS B " << VOLMIN << " " << FILLFACTOR << " " << NB << " " << Ngood <<"\n";
    // box_stats(mem,level,-10,SBoxes,SPoints);
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
    nPa=0;
    vector<Point*>Povers;
    for(vector<int> &SB : SBover)
      {
	for(int nz=SB[4];nz<=SB[5];nz+=spacing)
	  {
	    for(int ny=SB[2];ny<=SB[3];ny+=spacing)
	      {
		for(int nx=SB[0];nx<=SB[1];nx+=spacing)
		  {
		    array<int,4>ar4={{nx,ny,nz,nPa}};
		    std::map<array<int,4>,Point*>::iterator it=dupes.find(ar4);
		    assert(it != dupes.end());
		    if(it->second != 0)
		      continue;
		    Point* pp=new Point;
		    dupes[ar4]=pp;
		    pp->set_really_passive(true);
		    pp->set_pos_point(nx,ny,nz);
		  }
	      }
	  }
	nPa++;
      }
    vector<vector<Point*>>hypre_points(SBover.size());
    SBover.clear();
    SPover.clear();
    for(auto &d : dupes)
      {
	// if(d.second != 0)
	  hypre_points[d.first[3]].push_back(d.second);
      }
    dupes.clear();
    if(mem.p_mess->FractalRank == -21)
      cerr << " BOX PARAMS BB " << VOLMIN << " " << FILLFACTOR << "\n";
    hypre_points_boxes(mem,hypre_points,spacing,1,2.0,SBoxes,SPoints);
    hypre_points_zero(SPoints);
    if(mem.p_mess->FractalRank == -21)
      cerr << " BOX PARAMS CC " << VOLMIN << " " << FILLFACTOR << "\n";
    FHT << " BOX PARAMS C " << VOLMIN << " " << FILLFACTOR << "\n";
    // box_stats(mem,level,-20,SBoxes,SPoints);
  }
}
