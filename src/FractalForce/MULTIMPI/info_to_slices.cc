#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void info_to_slices(Fractal_Memory& mem,Fractal& frac,const int& lev)
  {
    int FractalNodes=mem.p_mess->FractalNodes;
    int FractalRank=mem.p_mess->FractalRank;
    int zoom=Misc::pow(2,frac.get_level_max());
    int length_1=frac.get_grid_length();
    vector <int>counts_out(FractalNodes);
    counts_out.assign(FractalNodes,0);
    vector <vector <int> > dataI_out(FractalNodes);
    vector <vector <double> > dataR_out(FractalNodes);
    int division=Misc::pow(2,frac.get_level_max()-lev);
    int wrapping=length_1*division;
    int really_long=length_1*zoom;
    vector <int> pos_point(3);
    int group_number=0;  
    for(vector <Group*>::const_iterator group_itr=mem.all_groups[lev].begin();
	group_itr!=mem.all_groups[lev].end();group_itr++)
      {
	Group& group=**group_itr;
	int point_number=0;
	for(vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr != group.list_points.end();++point_itr)
	  {
	    Point& point=**point_itr;
	    point.get_pos_point(pos_point);
	    int p_xi=((pos_point[0]+really_long) % wrapping)/division;
	    int p_yi=((pos_point[1]+really_long) % wrapping)/division;
	    int p_zi=((pos_point[2]+really_long) % wrapping)/division;
	    int S=mem.p_mess->WhichSlice[p_xi];
	    int slice_point=frac.where(p_xi,p_yi,p_zi,mem.p_mess->BoxS[S],mem.p_mess->BoxSL[S]);
	    dataI_out[S].push_back(slice_point);
	    dataI_out[S].push_back(group_number);
	    dataI_out[S].push_back(point_number);
	    counts_out[S]++;
	    point_number++;
	  }
	group_number++;
      }
    vector <int>counts_in(FractalNodes);
    vector <int> dataI_in;
    vector <double> dataR_in;
    int how_manyI=-1;
    int how_manyR=-1;
    int integers=3;
    int doubles=0;
    mem.p_file->note(true," info to slices a ");
    mem.p_mess->Send_Data_Some_How(counts_out,counts_in,integers,doubles,
				   dataI_out,dataI_in,how_manyI,
				   dataR_out,dataR_in,how_manyR);
    mem.p_file->note(true," info to slices c ");
    dataI_out.clear();
    dataR_out.clear();      
    int Slice=FractalRank;
    vector <int> BoxS;
    BoxS=mem.p_mess->BoxS[Slice];
    vector <int> BoxSL;
    BoxSL=mem.p_mess->BoxSL[Slice];
    int counterI=0;
    int number=0;
    int how_many=how_manyI/3;
    mem.p_mess->what_Slice_point.resize(how_many);
    mem.p_mess->return_group.resize(how_many);
    mem.p_mess->return_point.resize(how_many);
    mem.p_mess->return_node.resize(how_many);
    for(int FR=0;FR<FractalNodes;FR++)
      {
	for(int c=0;c<counts_in[FR];c++)
	  {
	    mem.p_mess->what_Slice_point[number]=dataI_in[counterI];
	    mem.p_mess->return_group[number]=dataI_in[counterI+1];
	    mem.p_mess->return_point[number]=dataI_in[counterI+2];
	    mem.p_mess->return_node[number]=FR;
	    counterI+=3;
	    number++;
	  }
      }
  }
}
