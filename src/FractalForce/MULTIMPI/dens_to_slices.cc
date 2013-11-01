#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void dens_to_slices(Group& group,Fractal_Memory& mem,Fractal& frac)
  {
    int FractalRank=mem.p_mess->FractalRank;
    int FractalNodes=mem.p_mess->FractalNodes;
    ofstream& FF=mem.p_file->FileFractal;
    int zoom=Misc::pow(2,frac.get_level_max());
    int length_1=frac.get_grid_length();
    int length_S=length_1;
    bool period=frac.get_periodic();
    if(!period)
      {
	length_1++;
	length_S*=2;
      }
    int length_2=length_1*length_1;
    vector <int> counts_out(FractalNodes);
    counts_out.assign(FractalNodes,0);
    vector <vector <int> > dataI_out(FractalNodes);
    vector <vector <double> > dataR_out(FractalNodes);
    vector <int> pos_point(3);
    for(vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr != group.list_points.end();++point_itr)
      {
	Point* p_point=*point_itr;
	bool pass=p_point->get_passive_point();
	p_point->get_pos_point(pos_point);
	int nx=(pos_point[0]/zoom+length_1) % length_1;
	int ny=(pos_point[1]/zoom+length_1) % length_1;
	int nz=(pos_point[2]/zoom+length_1) % length_1;
	int S=mem.p_mess->WhichSlice[nx];
	int global_pos=nx+(ny+nz*length_1)*length_1;
	if(pass)
	  global_pos=-global_pos-1;
	int numberFR=p_point->get_number_in_list();
	double density=p_point->get_density_point();
	dataI_out[S].push_back(global_pos);
	dataI_out[S].push_back(numberFR);
	dataR_out[S].push_back(density);
	counts_out[S]++;
      }
    vector <int> counts_in(FractalNodes);
    vector <int> dataI_in;
    vector <double> dataR_in;
    int how_manyI=-1;
    int how_manyR=-1;
    int integers=2;
    int doubles=1;
    frac.timing(-1,34);
    mem.p_mess->Full_Stop();
    frac.timing(1,34);
    mem.p_file->note(true," dens to slices a ");
    FF << " countsa " << counts_out[0] << " " << counts_out[1] << endl;
    mem.p_mess->How_Many_Things_To_Send(counts_out,counts_in);
    mem.p_file->note(true," dens to slices b ");
    FF << " countsb " << counts_out[0] << " " << counts_out[1] << endl;
    FF << " countsc " << counts_in[0] << " " << counts_in[1] << endl;
    mem.p_mess->Send_Data_Somewhere_No_Block(counts_out,counts_in,integers,doubles,
				    dataI_out,dataI_in,how_manyI,
				    dataR_out,dataR_in,how_manyR);
    mem.p_file->note(true," dens to slices c ");
    dataI_out.clear();
    dataR_out.clear();      
    int Slice=FractalRank;
    vector <int> BoxS;
    BoxS=mem.p_mess->BoxS[Slice];
    vector <int> BoxSL;
    BoxSL=mem.p_mess->BoxSL[Slice];
    int nx,ny,nz;
    int counterR=0;
    int counterI=0;
    int number=-1;
    mem.p_mess->return_point.resize(how_manyR);
    mem.p_mess->return_node.resize(how_manyR);
    mem.p_mess->what_Slice_point.resize(how_manyR);
    for(int FR=0;FR<FractalNodes;FR++)
      {
	for(int c=0;c<counts_in[FR];c++)
	  {
	    int n=dataI_in[counterI];
	    bool pass=n < 0;
	    if(pass)
	      n=-n-1;
	    nx=n % length_1;
	    ny=(n/length_1) % length_1;
	    nz=n/length_2;
	    number=frac.where(nx,ny,nz,BoxS,BoxSL);
	    if(!pass)
	      mem.p_mess->potR[number]=dataR_in[counterR];
	    mem.p_mess->return_point[counterR]=dataI_in[counterI+1];
	    mem.p_mess->return_node[counterR]=FR;
	    mem.p_mess->what_Slice_point[counterR]=number;
	    //	    bool something=dataR_in[counterR] != 0.0;
	    //	    cout << " slice " << Slice << " " << something << " " << nx << " " << ny << " " << nz << " " << n << " " << number << " " << counterR << " " << dataR_in[counterR] << " " << mem.p_mess->potR[number] << endl;
	    counterR++;
	    counterI+=2;
	  }
      }
  }
}
