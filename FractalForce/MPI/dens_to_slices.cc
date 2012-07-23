#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void dens_to_slices(Group& group,Fractal_Memory& mem,Fractal& frac)
  {
    int FractalRank=mem.p_mess->FractalRank;
    //    int FractalNodes=mem.p_mess->FractalNodes;
    vector <int>Box(6);
    frac.getBox(Box);
    int BSNodes=mem.Box_to_Slices[FractalRank].size();
    vector < vector <double> >dens;
    dens.resize(BSNodes);
    for(int BS=0;BS<BSNodes;BS++)
      {
	int total=0;
	int Slice=mem.Box_to_Slices[FractalRank][BS];
	for(int nx=mem.Box_to_Slices_Boxes[FractalRank][BS][0];nx<=mem.Box_to_Slices_Boxes[FractalRank][BS][1];nx++)
	  {
	    for(int ny=Box[2];ny<=Box[3];ny++)
	      {
		for(int nz=Box[4];nz<=Box[5];nz++)
		  {
		    dens[BS].push_back(group.list_points[frac.where_1(nx,ny,nz)]->get_density_point());
		    total++;
		  }
	      }
	  }
	bool first=BS==0;
	bool last=BS==BSNodes-1;
	mem.p_mess->send_data(dens[BS],Slice,total,first,last);
      }
  }
}
