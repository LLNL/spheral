#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void receive_potf(Group& group,Fractal_Memory& mem,Fractal& frac)
  {
    int FractalRank=mem.p_mess->FractalRank;
    vector <int>BBox(6);
    BBox=mem.BBoxes[FractalRank];
    int SNodes=mem.Box_from_Slices[FractalRank].size();
    vector < vector <double> >potf;
    potf.resize(SNodes);
    vector <int>buffsize(SNodes);
    vector <int>fromSlices(SNodes);
    for(int S=0;S<SNodes;S++)
      {
	int Slice=mem.Box_from_Slices[FractalRank][S];
	fromSlices[S]=Slice;
	buffsize[S]=(BBox[5]-BBox[4]+1)*(BBox[3]-BBox[2]+1);
	int nx0=mem.Slice_to_Boxes_Boxes[Slice][FractalRank][0];
	int nx1=mem.Slice_to_Boxes_Boxes[Slice][FractalRank][1];
	buffsize[S]*=nx1-nx0+1;
	potf[S].resize(buffsize[S]);
      }
    mem.p_mess->recv_data(potf,fromSlices,buffsize);

    for(int S=0;S<SNodes;S++)
      {
	int total=0;
	int Slice=mem.Box_from_Slices[FractalRank][S];
	for(int nx=mem.Slice_to_Boxes_Boxes[Slice][FractalRank][0];nx<=mem.Slice_to_Boxes_Boxes[Slice][FractalRank][1];nx++)
	  {
	    for(int ny=BBox[2];ny<=BBox[3];ny++)
	      {
		for(int nz=BBox[4];nz<=BBox[5];nz++)
		  {
		    int number=frac.where_1(nx,ny,nz);
		    group.list_points[number]->set_potential_point(potf[S][total]);
		    total++;
		  }
	      }
	  }
      }
  }
}
