#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void dens_to_slices(Group& group,Fractal_Memory& mem,Fractal& frac)
  {
    static bool printit=true;
    printit=false;
    int FractalNodes=mem.p_mess->FractalNodes;
    ofstream& FF=mem.p_file->DUMPS;
    int zoom=Misc::pow(2,frac.get_level_max());
    int length_1=frac.get_grid_length();
    int length_S=length_1;
    int length_S2=0;
    bool period=frac.get_periodic();
    if(!period)
      {
	length_1++;
	length_S*=2;
	length_S2=length_S+2;
      }
    vector <int> pos_point(3);
    fprintf(mem.p_file->PFTime," dens to slices ");
    int LOOPS=((length_1*length_1)/(512*512))+1;
    for(int LOOP=0;LOOP<LOOPS;LOOP++)
      {
	double time1=-mem.p_mess->Clock();
	vector <int>counts_in;
	vector <int>counts_out(FractalNodes,0);;
	vector <vector <int> > dataI_out(FractalNodes);
	vector <vector <double> > dataR_out(FractalNodes);
	vector <int> dataI_in;
	vector <double> dataR_in;
	int loop_count=0;
	for(auto &p_point : group.list_points)
	  {
	    bool doit=loop_count % LOOPS == LOOP && !p_point->get_passive_point() && p_point->get_mass_point();
	    if(doit)
	      {
		p_point->get_pos_point(pos_point);
		int nx=(pos_point[0]/zoom+length_1) % length_1;
		int ny=(pos_point[1]/zoom+length_1) % length_1;
		int nz=(pos_point[2]/zoom+length_1) % length_1;
		int S=mem.p_mess->WhichSlice[nx];
		dataI_out[S].push_back(frac.where(nx,ny,nz,mem.p_mess->BoxS[S],mem.p_mess->BoxSL[S]));
		dataR_out[S].push_back(p_point->get_density_point());
		if(printit)
		  FF << " DSA" << LOOP << " " << loop_count << " " << pos_point[0] << " "<<  pos_point[1] << " " << pos_point[2] << " " << S << " " << dataI_out[S].back() << " " << dataR_out[S].back() << "\n";
		counts_out[S]++;
	      }
	    loop_count++;
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
	clean_vector(counts_out);
	dataI_out.clear();
	dataR_out.clear();      
	if(LOOP == 0)
	  {
	    mem.p_mess->create_potRS();
	    mem.p_mess->zeroRS(-frac.get_density_0());
	  }
	int counterIR=0;
	for(int FR=0;FR<FractalNodes;FR++)
	  {
	    for(int c=0;c<counts_in[FR];c++)
	      {
		int NN=dataI_in[counterIR];
		if(!period)
		  {
		    int nz=NN % length_S2;
		    int ny=(NN/length_S2) % length_S;
		    int nx=NN/(length_S2*length_S);
		    assert(nx < length_1);
		    assert(ny < length_1);
		    assert(nz < length_1);
		    NN=nz+(ny+nx*length_1)*length_1;
		  }
		mem.p_mess->potRS[NN]=dataR_in[counterIR];
		if(printit)
		  FF << " DSB " << counterIR << " " << dataI_in[counterIR] << " " << NN << " " << dataR_in[counterIR] << "\n";
		counterIR++;
	      }
	  }
	time1+=mem.p_mess->Clock();
	fprintf(mem.p_file->PFTime," %3d %8.3E ",LOOP,time1);
      }
    fprintf(mem.p_file->PFTime,"\n");
    if(printit)
      FF << endl;
    printit=false;
  }
}
