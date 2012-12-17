#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void slices_to_potf(Fractal_Memory& mem,Fractal& frac,const int& length)
  {
    int Slice=mem.p_mess->FractalRank;
    vector <int>BBox(6);
    vector <int> BoxS(6);
    BoxS[0]=mem.p_mess->Slices[Slice][0];
    BoxS[1]=mem.p_mess->Slices[Slice][1];
    BoxS[2]=0;
    BoxS[3]=length-1;
    BoxS[4]=0;
    BoxS[5]=length-1;
    vector <int> BoxSL(3);
    BoxSL[0]=BoxS[1]-BoxS[0]+1;
    BoxSL[1]=length;
    BoxSL[2]=length;
    int BNodes=mem.Slice_to_Boxes[Slice].size();
    vector < vector <double> >pot;
    pot.resize(BNodes);
    for(int B=0;B<BNodes;B++)
      {
	int FR=mem.Slice_to_Boxes[Slice][B];
	BBox=mem.BBoxes[FR];
	int buffsize=(BBox[5]-BBox[4]+1)*(BBox[3]-BBox[2]+1);
	int nx0=mem.Slice_to_Boxes_Boxes[Slice][FR][0];
	int nx1=mem.Slice_to_Boxes_Boxes[Slice][FR][1];
	buffsize*=nx1-nx0+1;
	pot[B].resize(buffsize);
	int total=0;
	int number=-1;
	for(int nx=nx0;nx<=nx1;nx++)
	  {
	    for(int ny=BBox[2];ny<=BBox[3];ny++)
	      {
		for(int nz=BBox[4];nz<=BBox[5];nz++)
		  {

		    number=frac.where(nx,ny,nz,BoxS,BoxSL);
		    pot[B].push_back(mem.p_mess->potR[number]);
		    total++;
		  }
	      }
	  }
	bool first=B==0;
	bool last=B==BNodes-1;
	mem.p_mess->send_data(pot[B],FR,total,first,last);
      }
  }
}
