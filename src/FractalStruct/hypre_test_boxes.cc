#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_test_boxes(Fractal_Memory& mem,int spacing,
			vector < vector<int> >& SBoxes,vector < vector<Point*> >& SPoints)
  {
    assert(SBoxes.size() == SPoints.size());
    int nBa=0;
    for(auto &SB : SBoxes)
      {
	int vola=(SB[1]-SB[0])/spacing+1;
	int volb=(SB[3]-SB[2])/spacing+1;
	int volc=(SB[5]-SB[4])/spacing+1;
	assert(!SPoints[nBa].empty());
	assert(vola > 0 && volb > 0 && volc > 0 && vola*volb*volc == SPoints[nBa].size());
	int nSa=0;
	for(int nz=SB[4];nz<=SB[5];nz+=spacing)
	  {
	    for(int ny=SB[2];ny<=SB[3];ny+=spacing)
	      {
		for(int nx=SB[0];nx<=SB[1];nx+=spacing)
		  {
		    Point* p=SPoints[nBa][nSa];
		    if(p != 0)
		      {
			assert(nx == p->get_pos_point_x());
			assert(ny == p->get_pos_point_y());
			assert(nz == p->get_pos_point_z());
			assert(p->get_inside());
			for(int ni=0;ni<6;ni++)
			  assert(p->get_point_ud(ni));
		      }
		    nSa++;
		  }
	      }
	  }
	nBa++;
      }
    nBa=0;
    bool baad=false;
    for(auto &SBa : SBoxes)
      {
	for(int nBb=nBa+1;nBb<SBoxes.size();nBb++)
	  if(overlap_boxes(SBa,SBoxes[nBb]))
	    {
	      baad=true;
	      cerr << " BOX OVERLAP " << mem.p_mess->FractalRank << " " << nBa << " " << nBb << endl;
	      cerr << SBa[0] <<  " " << SBa[1] <<  " " << SBa[2] <<  " " << SBa[3] <<  " " << SBa[4] <<  " " << SBa[5] << endl;
	      cerr << SBoxes[nBb][0] <<  " " << SBoxes[nBb][1] <<  " " << SBoxes[nBb][2] <<  " " << SBoxes[nBb][3] <<  " " << SBoxes[nBb][4] <<  " " << SBoxes[nBb][5] << endl;
	    }
	nBa++;
      }
    mem.p_mess->Full_Stop_Do_Not_Argue();
    assert(!baad);
  }
}
