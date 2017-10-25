#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void slices_to_potf(Fractal_Memory& mem,Fractal& frac,int lev)
  {
    bool period=frac.get_periodic();
    int FractalNodes=mem.p_mess->FractalNodes;
    int zoom=Misc::pow(2,frac.get_level_max());
    int length_1=frac.get_grid_length();
    int length_11=length_1+1;
    int length_S=2*length_1;
    int length_S2=length_S+2;
    int division=Misc::pow(2,frac.get_level_max()-lev);
    int wrapping=length_1*division;
    int really_long=length_1*zoom;
    if(!period)
      {
	wrapping*=2;
	really_long=0;
      }
    vector <int> pos_point(3);
    bool notZERO= lev > 0;
    fprintf(mem.p_file->PFTime," slices to potf "); 
    int LOOPS=((length_1*length_1)/(512*512))+1;
    if(lev > 0)
      LOOPS*=3;
    for(int LOOP=0;LOOP<LOOPS;LOOP++)
      {
	double time1=-mem.p_mess->Clock();
	vector <int>counts_in;
	vector <int>counts_out(FractalNodes,0);
	vector <vector <int> > dataI_out(FractalNodes);
	vector <vector <double> > dataR_out(FractalNodes);
	vector <int> dataI_in;
	vector <double> dataR_in;
	vector<vector <Point*>> point_counter(FractalNodes);
	int number=0;
	for(auto &pgroup : mem.all_groups[lev])
	  {
	    Group& group=*pgroup;
	    for(auto &ppoint : group.list_points)
	      {
		if(number % LOOPS == LOOP)
		  {
		    ppoint->get_pos_point(pos_point);
		    int p_xi=((pos_point[0]+really_long) % wrapping)/division;
		    int p_yi=((pos_point[1]+really_long) % wrapping)/division;
		    int p_zi=((pos_point[2]+really_long) % wrapping)/division;
		    int S=mem.p_mess->WhichSlice[p_xi];
		    int slice_point=frac.where(p_xi,p_yi,p_zi,mem.p_mess->BoxS[S],mem.p_mess->BoxSL[S]);
		    dataI_out[S].push_back(slice_point);
		    point_counter[S].push_back(ppoint);
		    counts_out[S]++;
		  }
		number++;
	      }
	  }
	int how_manyI=-1;
	int how_manyR=-1;
	int integers=1;
	int doubles=0;
	mem.p_mess->Send_Data_Some_How(4,counts_out,counts_in,integers,doubles,
				       dataI_out,dataI_in,how_manyI,
				       dataR_out,dataR_in,how_manyR);
	dataI_out.clear();
	dataR_out.clear();    
	dataI_out.resize(FractalNodes);
	dataR_out.resize(FractalNodes);
	int counterI=0;
	for(int FR=0;FR<FractalNodes;FR++)
	  {
	    for(int c=0;c<counts_in[FR];c++)
	      {
		int NN=dataI_in[counterI];
		if(!period)
		  {
		    int nz=NN % length_S2;
		    int ny=(NN/length_S2) % length_S;
		    int nx=NN/(length_S2*length_S);
		    assert(nx < length_1);
		    assert(ny < length_1);
		    assert(nz < length_1);
		    NN=nz+(ny+nx*length_11)*length_11;
		  }
		dataR_out[FR].push_back(mem.p_mess->potRS[NN]);
		counterI+=integers;
	      }
	  }
	counts_out=counts_in;
	counts_in.assign(FractalNodes,0);
	clean_vector(dataI_in);
	clean_vector(dataR_in);
	how_manyI=-1;
	how_manyR=-1;
	integers=0;
	doubles=1;
	mem.p_mess->Send_Data_Some_How(7,counts_out,counts_in,integers,doubles,
				       dataI_out,dataI_in,how_manyI,
				       dataR_out,dataR_in,how_manyR);
	dataI_out.clear();
	dataR_out.clear();      
	int counterIR=0;
	for(int FR=0;FR<FractalNodes;FR++)
	  {
	    for(int c=0;c<counts_in[FR];c++)
	      {
		double potential=dataR_in[counterIR];
		Point* p_point=point_counter[FR][c];
		if(notZERO)
		  potential+=p_point->get_potential_point();
		p_point->set_potential_point(potential);
		counterIR++;
	      }
	  }
	time1+=mem.p_mess->Clock();
	fprintf(mem.p_file->PFTime," %3d %8.3E ",LOOP,time1);
      }
    mem.p_mess->free_potRS();
    fprintf(mem.p_file->PFTime,"\n");
  }
}
