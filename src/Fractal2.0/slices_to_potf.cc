#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void slices_to_potf(Fractal_Memory& mem,Fractal& frac,int lev)
  {
    ofstream& FF=mem.p_file->DUMPS;
    int FractalNodes=mem.p_mess->FractalNodes;
    int zoom=Misc::pow(2,frac.get_level_max());
    int length_1=frac.get_grid_length();
    int division=Misc::pow(2,frac.get_level_max()-lev);
    int wrapping=length_1*division;
    int really_long=length_1*zoom;
    vector <int> pos_point(3);
    bool notZERO= lev > 0;
    vector <int>counts_in;
    vector <int>counts_out;
    vector <vector <int> > dataI_out;
    vector <vector <double> > dataR_out;
    vector <int> dataI_in;
    vector <double> dataR_in;
    int LOOPS=1;
    vector <double>times;
    times.push_back(mem.p_mess->Clock());
    dataI_out.clear();
    dataR_out.clear();
    dataI_in.clear();
    dataR_in.clear();
    dataI_out.resize(FractalNodes);
    dataR_out.resize(FractalNodes);
    counts_out.assign(FractalNodes,0);
    counts_in.assign(FractalNodes,0);
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
    int how_manyI=-1;
    int how_manyR=-1;
    int integers=3;
    int doubles=0;
    mem.p_file->note(true," info to slices a ");
    mem.p_mess->Send_Data_Some_How(4,counts_out,counts_in,integers,doubles,
				   dataI_out,dataI_in,how_manyI,
				   dataR_out,dataR_in,how_manyR);
    mem.p_file->note(true," info to slices c ");
    times.push_back(mem.p_mess->Clock());
    dataI_out.clear();
    dataR_out.clear();    
    dataI_out.resize(FractalNodes);
    dataR_out.resize(FractalNodes);
    int counterI=0;
    for(int FR=0;FR<FractalNodes;FR++)
      {
	for(int c=0;c<counts_in[FR];c++)
	  {
	    dataR_out[FR].push_back(mem.p_mess->potR[dataI_in[counterI]]);
	    dataI_out[FR].push_back(dataI_in[counterI+1]);
	    dataI_out[FR].push_back(dataI_in[counterI+2]);
	    //	    FF << " SPA " << counterI << " " << dataI_in[counterI] << " " << mem.p_mess->potR[dataI_in[counterI]] << " " << dataI_in[counterI+1] << "\n";
	    counterI+=integers;
	  }
      }
    
    counts_out=counts_in;
    counts_in.assign(FractalNodes,0);
    dataI_in.clear();
    dataR_in.clear();
    //
    how_manyI=-1;
    how_manyR=-1;
    integers=2;
    doubles=1;
    mem.p_file->note(true," slices to potf a ");
    mem.p_mess->Send_Data_Some_How(7,counts_out,counts_in,integers,doubles,
				   dataI_out,dataI_in,how_manyI,
				   dataR_out,dataR_in,how_manyR);
    mem.p_file->note(true," slices to potf c ");
    dataI_out.clear();
    dataR_out.clear();      
    
    counterI=0;
    int counterR=0;
    for(int FR=0;FR<FractalNodes;FR++)
      {
	for(int c=0;c<counts_in[FR];c++)
	  {
	    int number_group=dataI_in[counterI];
	    int number_point=dataI_in[counterI+1];
	    double potential=dataR_in[counterR];
	    Point* p_point=mem.all_groups[lev][number_group]->list_points[number_point];
	    if(notZERO)
	      potential+=p_point->get_potential_point();
	    p_point->set_potential_point(potential);
	    //	    FF << " SPB " << counterI << " " << dataI_in[counterI] << " " << dataR_in[counterR] << "\n";
	    counterI+=integers;
	    counterR++;
	  }
      }
    times.push_back(mem.p_mess->Clock());
    fprintf(mem.p_file->PFTime," slices to potf "); 
    for(int ni=0;ni<2*LOOPS;ni++)
      fprintf(mem.p_file->PFTime,"%3d %9.3E ",ni/2,times[ni+1]-times[ni]);
    fprintf(mem.p_file->PFTime,"\n");
    mem.p_file->note(true," slices to potf exit ");
  }
}
