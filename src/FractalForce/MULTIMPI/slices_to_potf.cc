#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void slices_to_potf(Group& group,Fractal_Memory& mem,Fractal& frac)
  {
    int FractalNodes=mem.p_mess->FractalNodes;
    int length_1=frac.get_grid_length();
    int length_S=length_1;
    bool period=frac.get_periodic();
    if(!period)
      {
	length_1++;
	length_S*=2;
      }
    vector <int>counts_out(FractalNodes);
    counts_out.assign(FractalNodes,0);
    int number=-1;
    double potential=-1.0;
    vector < vector <int> > dataI_out(FractalNodes);
    vector < vector <double> > dataR_out(FractalNodes);
    int number_points=mem.p_mess->return_point.size();
    for(int number=0;number<number_points;number++)
      {
	int numberS=mem.p_mess->what_Slice_point[number];
	int number_list=mem.p_mess->return_point[number];
	int FR=mem.p_mess->return_node[number];
	dataR_out[FR].push_back(mem.p_mess->potR[numberS]);
	dataI_out[FR].push_back(number_list);
	counts_out[FR]++;
      }
    mem.p_mess->return_Slice_pos.clear();
    mem.p_mess->return_group.clear();
    mem.p_mess->return_point.clear();
    mem.p_mess->return_node.clear();
    vector <int>counts_in(FractalNodes);
    vector <int> dataI_in;
    vector <double> dataR_in;
    int how_manyI=-1;
    int how_manyR=-1;
    int integers=1;
    int doubles=1;
    mem.p_file->note(true," slices to potf a ");
    mem.p_mess->Send_Data_Some_How(counts_out,counts_in,integers,doubles,
				   dataI_out,dataI_in,how_manyI,
				   dataR_out,dataR_in,how_manyR);
    mem.p_file->note(true," slices to potf c ");
    dataI_out.clear();
    dataR_out.clear();      
    number=-1;
    int counter=0;
    potential=-1.0;
    for(int FR=0;FR<FractalNodes;FR++)
      {
	for(int c=0;c<counts_in[FR];c++)
	  {
	    number=dataI_in[counter];
	    potential=dataR_in[counter];
	    group.list_points[number]->set_potential_point(potential);
	    counter++;
	  }
      }
  }
}
