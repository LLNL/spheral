#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void slices_to_pot_init(Fractal_Memory& mem,Fractal& frac,const int& lev)
  {
    mem.p_file->note(true," slices to pot init enter ");
    int FractalNodes=mem.p_mess->FractalNodes;
    //    int length_1=frac.get_grid_length();
    //    int length_S=length_1;
    vector <int> counts_out(FractalNodes);
    counts_out.assign(FractalNodes,0);
    double potential=-1.0;
    vector < vector <int> > dataI_out(FractalNodes);
    vector < vector <double> > dataR_out(FractalNodes);
    int number_points=mem.p_mess->return_point.size();
    for(int number=0;number<number_points;number++)
      {
	int numberS=mem.p_mess->what_Slice_point[number];
	int number_list=mem.p_mess->return_point[number];
	int number_group=mem.p_mess->return_group[number];
	int FR=mem.p_mess->return_node[number];
	dataR_out[FR].push_back(mem.p_mess->potR[numberS]);
	dataI_out[FR].push_back(number_group);
	dataI_out[FR].push_back(number_list);
	counts_out[FR]++;
      }
    mem.p_mess->what_Slice_point.clear();
    mem.p_mess->return_group.clear();
    mem.p_mess->return_point.clear();
    mem.p_mess->return_node.clear();
    vector <int> counts_in(FractalNodes);
    vector <int> dataI_in;
    vector <double> dataR_in;
    int how_manyI=-1;
    int how_manyR=-1;
    int integers=2;
    int doubles=1;
    frac.timing(-1,40);
    mem.p_mess->Full_Stop();
    frac.timing(1,40);
    mem.p_file->note(true," slices to pot init a ");
    mem.p_mess->How_Many_Things_To_Send(counts_out,counts_in);
    mem.p_file->note(true," slices to pot init b ");
    mem.p_mess->Send_Data_Somewhere_No_Block(counts_out,counts_in,integers,doubles,
				    dataI_out,dataI_in,how_manyI,
				    dataR_out,dataR_in,how_manyR);
    mem.p_file->note(true," slices to pot init c ");
    dataI_out.clear();
    dataR_out.clear();      

    int number_group=-1;
    int number_point=1;
    int counterI=0;
    int counterR=0;
    potential=-1.0;
    Point* p_point=0;
    for(int FR=0;FR<FractalNodes;FR++)
      {
	for(int c=0;c<counts_in[FR];c++)
	  {
	    number_group=dataI_in[counterI];
	    number_point=dataI_in[counterI+1];
	    potential=dataR_in[counterR];
	    //	    mem.p_file->FileFractal << " check it " << lev << " " << FR << " " << c << " " << number_group << " " << number_point << endl;
	    p_point=mem.all_groups[lev][number_group]->list_points[number_point];
	    if(lev != 0)
	      potential+=p_point->get_potential_point();
	    p_point->set_potential_point(potential);
	    counterI+=2;
	    counterR++;
	  }
      }
    mem.p_file->note(true," slices to pot init exit ");
  }
}
