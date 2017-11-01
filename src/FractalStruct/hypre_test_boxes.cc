#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void hypre_test_boxes(Fractal_Memory& mem,int level,
			vector < vector<int> >& SBoxes,vector < vector<Point*> >& SPoints)
  {
    assert(SBoxes.size() == SPoints.size());
    vector <int>pos(3);
    vector <int>pos1(3);
    ofstream& FHT=mem.p_file->DUMPS;
    int badtouches(0);
    int goodtouches(0);
    const int spacing=Misc::pow(2,mem.p_fractal->get_level_max()-level);
    const int FractalRank=mem.p_mess->FractalRank;
    vector <int> BOX=mem.BoxesLev[FractalRank][level];
    FHT << " MYBOX " << level << " " << spacing << " " << BOX[0] << " " << BOX[1] << " " << BOX[2];
    FHT << " " << BOX[3] << " " << BOX[4] << " " << BOX[5] << "\n";
    for(int FR : mem.Touchy)
      {
	vector <int>BOXT=mem.BoxesLev[FR][level];
	FHT << " TBOX " << level << " " << FR << " " << BOXT[0] << " " << BOXT[1] << " " << BOXT[2];
	FHT << " " << BOXT[3] << " " << BOXT[4] << " " << BOXT[5] << "\n";
	vector <int>BBOXT=mem.BBoxesLev[FR][level];
	FHT << " TBBOX " << level << " " << FR << " " << BBOXT[0] << " " << BBOXT[1] << " " << BBOXT[2];
	FHT << " " << BBOXT[3] << " " << BBOXT[4] << " " << BBOXT[5] << "\n";
      }
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
			  {
			    Point*p1=p->get_point_ud_0(ni,-17);
			    p1->get_pos_point(pos1);
			    if(vector_in_box(pos1,BOX))
			      continue;
			    int touches=0;
			    for(int FR : mem.Touchy)
			      {
				if(vector_in_box(pos1,mem.BoxesLev[FR][level]))
				  touches++;				
			      }
			    if(touches==1)
			      {
				goodtouches++;
				continue;
			      }
			    badtouches++;
			    FHT << " MYBOX " << level << " " << BOX[0] << " " << BOX[1] << " " << BOX[2];
			    FHT << " " << BOX[3] << " " << BOX[4] << " " << BOX[5] << "\n";
			    p->get_pos_point(pos);
			    FHT << " MYPOS0 " << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
			    FHT << " MYPOS1 " << pos1[0] << " " << pos1[1] << " " << pos1[2] << "\n";
			    if(touches > 1)
			      {
				for(int FR : mem.Touchy)
				  {
				    vector <int>BOXFR=mem.BoxesLev[FR][level];
				    if(vector_in_box(pos1,BOXFR))
				      {
					touches++;
					FHT << " OTHERBOX "  << " " << BOXFR[0] << " " << BOXFR[1] << " " << BOXFR[2];
					FHT << " " << BOXFR[3] << " " << BOXFR[4] << " " << BOXFR[5] << "\n";
				      }
				  }
			      }
			  }
			
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
	      cerr << " BOX OVERLAP " << mem.p_mess->FractalRank << " " << nBa << " " << nBb << "\n";
	      cerr << SBa[0] <<  " " << SBa[1] <<  " " << SBa[2] <<  " " << SBa[3] <<  " " << SBa[4] <<  " " << SBa[5] << "\n";
	      cerr << SBoxes[nBb][0] <<  " " << SBoxes[nBb][1] <<  " " << SBoxes[nBb][2] <<  " " << SBoxes[nBb][3] <<  " " << SBoxes[nBb][4] <<  " " << SBoxes[nBb][5] << "\n";
	    }
	nBa++;
      }
    FHT << " TOTALS " << level << " " << goodtouches << " " << badtouches << "\n";
    FHT.flush();
    mem.p_mess->Full_Stop_Do_Not_Argue();
    assert(!baad);
    assert(badtouches==0);
  }
}
