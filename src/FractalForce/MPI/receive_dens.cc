#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void receive_dens(Fractal_Memory& mem,Fractal& frac,const int& length)
  {
    int Slice=mem.p_mess->FractalRank;
    vector <int>Box(6);
    int BNodes=mem.Slice_from_Boxes[Slice].size();
    vector < vector <double> >dens;
    dens.resize(BNodes);
    vector <int>buffsize(BNodes);
    vector <int>fromBoxes(BNodes);
    for(int B=0;B<BNodes;B++)
      {
	int FR=mem.Slice_from_Boxes[Slice][B];
	fromBoxes[B]=FR;
	Box=mem.Boxes[FR];
	buffsize[B]=(Box[5]-Box[4]+1)*(Box[3]-Box[2]+1);
	int nx0=mem.Box_to_Slices_Boxes[FR][Slice][0];
	int nx1=mem.Box_to_Slices_Boxes[FR][Slice][1];
	buffsize[B]*=nx1-nx0+1;
	dens[B].resize(buffsize[B]);
      }
    mem.p_mess->recv_data(dens,fromBoxes,buffsize);
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
    int number=-1;
    for(int B=0;B<BNodes;B++)
      {
	int FR=mem.Box_from_Slices[Slice][B];
	int nx0=mem.Box_to_Slices_Boxes[FR][Slice][0];
	int nx1=mem.Box_to_Slices_Boxes[FR][Slice][1];
	int total=0;
	for(int nx=nx0;nx<=nx1;nx++)
	  {
	    for(int ny=Box[2];ny<=Box[3];ny++)
	      {
		for(int nz=Box[4];nz<=Box[5];nz++)
		  {
		    number=frac.where(nx,ny,nz,BoxS,BoxSL);
		    mem.p_mess->potR[number]=dens[B][total];
		    total++;
		  }
	      }
	  }
      } 
  }
}
