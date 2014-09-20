#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void dens_to_slices(Group& group,Fractal_Memory& mem,Fractal& frac)
  {
    int FractalRank=mem.p_mess->FractalRank;
    int FractalNodes=mem.p_mess->FractalNodes;
    ofstream& FF=mem.p_file->DUMPS;
    //    ofstream& FF=mem.p_file->FileFractal;
    int zoom=Misc::pow(2,frac.get_level_max());
    int length_1=frac.get_grid_length();
    int length_S=length_1;
    bool period=frac.get_periodic();
    if(!period)
      {
	length_1++;
	length_S*=2;
      }

    vector <int> pos_point(3);
    vector <int>counts_in;
    vector <int>counts_out;
    vector <vector <int> > dataI_out;
    vector <vector <double> > dataR_out;
    vector <int> dataI_in;
    vector <double> dataR_in;
    int LOOPS=(FractalNodes-1)/256+1;
    for(int loop=0;loop<LOOPS;loop++)
      {
	dataI_out.clear();
	dataR_out.clear();
	dataI_in.clear();
	dataR_in.clear();
	dataI_out.resize(FractalNodes);
	dataR_out.resize(FractalNodes);
	counts_out.assign(FractalNodes,0);
	counts_in.assign(FractalNodes,0);
	int loopcount=0;
	for(vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr != group.list_points.end();++point_itr)
	  {
	    if(loopcount % LOOPS == loop)
	      {
		Point* p_point=*point_itr;
		bool pass=p_point->get_passive_point();
		p_point->get_pos_point(pos_point);
		int nx=(pos_point[0]/zoom+length_1) % length_1;
		int ny=(pos_point[1]/zoom+length_1) % length_1;
		int nz=(pos_point[2]/zoom+length_1) % length_1;
		int S=mem.p_mess->WhichSlice[nx];
		int slice_point=frac.where(nx,ny,nz,mem.p_mess->BoxS[S],mem.p_mess->BoxSL[S]);
		if(pass)
		  slice_point=-slice_point-1;
		double density=p_point->get_density_point();
		dataI_out[S].push_back(slice_point);
		dataR_out[S].push_back(density);
		counts_out[S]++;
	      }
	  }
	int how_manyI=-1;
	int how_manyR=-1;
	int integers=1;
	int doubles=1;
	mem.p_file->note(true," dens to slices a ");
	mem.p_mess->Send_Data_Some_How(0,counts_out,counts_in,integers,doubles,
				       dataI_out,dataI_in,how_manyI,
				       dataR_out,dataR_in,how_manyR);
	mem.p_file->note(true," dens to slices c ");
	dataI_out.clear();
	dataR_out.clear();      
	int counterIR=0;
	for(int FR=0;FR<FractalNodes;FR++)
	  {
	    for(int c=0;c<counts_in[FR];c++)
	      {
		int n=dataI_in[counterIR];
		bool pass=n < 0;
		if(pass)
		  n=-n-1;
		if(!pass)
		  mem.p_mess->potR[n]=dataR_in[counterIR];
		counterIR++;
	      }
	  }
      }
  }
}
