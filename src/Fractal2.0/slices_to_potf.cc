#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void slices_to_potf(Fractal_Memory& mem,Fractal& frac,int lev)
  {
//     ofstream& FF=mem.p_file->DUMPS;
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
    bool simpleGROUPS= mem.all_groups[lev].size() == 1;
    vector <int>counts_in;
    vector <int>counts_out;
    vector <vector <int> > dataI_out;
    vector <vector <double> > dataR_out;
    vector <int> dataI_in;
    vector <double> dataR_in;
    vector <int> point_counter;
    //
//     counts_out.assign(FractalNodes,0);
//     for(vector <Group*>::const_iterator group_itr=mem.all_groups[lev].begin();
// 	group_itr!=mem.all_groups[lev].end();group_itr++)
//       {
// 	Group& group=**group_itr;
// 	for(vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr != group.list_points.end();++point_itr)
// 	  {
// 	    Point& point=**point_itr;
// 	    point.get_pos_point(pos_point);
// 	    int p_xi=((pos_point[0]+really_long) % wrapping)/division;
// 	    int S=mem.p_mess->WhichSlice[p_xi];
// 	    counts_out[S]++;
// 	  }
//       }
    fprintf(mem.p_file->PFTime," slices to potf "); 
//     double time1=-mem.p_mess->Clock();
//     vector <int>maxSR;
//     mem.p_mess->MAX_Things_To_Send_Receive_I(counts_out,counts_in,maxSR);
//     time1+=mem.p_mess->Clock();
//     fprintf(mem.p_file->PFTime," %3d %8.3E ",-1,time1);
//     int maxOUT=maxSR[0]*2;
//     int maxIN=maxSR[1]*2;
//     int maxINOUT=max(maxOUT,maxIN);
//     int maxIO=6000000;
//     counts_in.clear();
//     counts_out.clear();
//     //
//     int LOOPS=(maxINOUT-1)/maxIO+1;
//     FF << " LOOPS " << maxIO << " " << maxOUT << " " << maxIN << " " << LOOPS << "\n";
    int LOOPS=((length_1*length_1)/(512*512))+1;
    if(lev > 0)
      LOOPS*=3;
    //
    //    LOOPS=1;
    //
    for(int LOOP=0;LOOP<LOOPS;LOOP++)
      {
	double time1=-mem.p_mess->Clock();
	dataI_out.clear();
	dataR_out.clear();
	dataI_in.clear();
	dataR_in.clear();
	dataI_out.resize(FractalNodes);
	dataR_out.resize(FractalNodes);
	counts_out.assign(FractalNodes,0);
	counts_in.assign(FractalNodes,0);
	int number=0;  
	point_counter.clear();
	for(vector <Group*>::const_iterator group_itr=mem.all_groups[lev].begin();
	    group_itr!=mem.all_groups[lev].end();group_itr++)
	  {
	    Group& group=**group_itr;
	    point_counter.push_back(number);
	    for(vector <Point*>::const_iterator point_itr=group.list_points.begin();point_itr != group.list_points.end();++point_itr)
	      {
		if(number % LOOPS == LOOP)
		  {
		    Point& point=**point_itr;
		    point.get_pos_point(pos_point);
		    int p_xi=((pos_point[0]+really_long) % wrapping)/division;
		    int p_yi=((pos_point[1]+really_long) % wrapping)/division;
		    int p_zi=((pos_point[2]+really_long) % wrapping)/division;
		    int S=mem.p_mess->WhichSlice[p_xi];
		    int slice_point=frac.where(p_xi,p_yi,p_zi,mem.p_mess->BoxS[S],mem.p_mess->BoxSL[S]);
		    dataI_out[S].push_back(slice_point);
		    dataI_out[S].push_back(number);
		    counts_out[S]++;
		  }
		number++;
	      }
	  }
	int how_manyI=-1;
	int how_manyR=-1;
	int integers=2;
	int doubles=0;
	//	mem.p_file->note(true," info to slices a ");
	mem.p_mess->Send_Data_Some_How(4,counts_out,counts_in,integers,doubles,
				       dataI_out,dataI_in,how_manyI,
				       dataR_out,dataR_in,how_manyR);
	//	mem.p_file->note(true," info to slices c ");
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
		dataI_out[FR].push_back(dataI_in[counterI+1]);
		//	    FF << " SPA " << lev << " " << counterI << " " << dataI_in[counterI] << " " << NN << " " << mem.p_mess->potRS[NN] << " " << dataI_in[counterI+1] << " " << dataI_in[counterI+2] << "\n";
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
	integers=1;
	doubles=1;
	//	mem.p_file->note(true," slices to potf a ");
	mem.p_mess->Send_Data_Some_How(7,counts_out,counts_in,integers,doubles,
				       dataI_out,dataI_in,how_manyI,
				       dataR_out,dataR_in,how_manyR);
	//	mem.p_file->note(true," slices to potf c ");
	dataI_out.clear();
	dataR_out.clear();      
    
	int number_group=-1;
	int number_point=-1;
	int counterIR=0;
	for(int FR=0;FR<FractalNodes;FR++)
	  {
	    for(int c=0;c<counts_in[FR];c++)
	      {
		if(simpleGROUPS)
		  {
		    number_group=0;
		    number_point=dataI_in[counterIR];
		  }
		else
		  {
		    number_group=std::upper_bound(point_counter.begin(),point_counter.end(),dataI_in[counterIR])-point_counter.begin()-1;
		    number_point=dataI_in[counterIR]-point_counter[number_group];
		  }
		double potential=dataR_in[counterIR];
		Point* p_point=mem.all_groups[lev][number_group]->list_points[number_point];
		if(notZERO)
		  potential+=p_point->get_potential_point();
		p_point->set_potential_point(potential);
		//	    FF << " SPB " << lev << " " << counterI << " " << dataI_in[counterI] << " " << dataI_in[counterI+1] << " " << dataR_in[counterR] << "\n";
		counterIR++;
	      }
	  }
	time1+=mem.p_mess->Clock();
	fprintf(mem.p_file->PFTime," %3d %8.3E ",LOOP,time1);
      }
    mem.p_mess->free_potRS();
    fprintf(mem.p_file->PFTime,"\n");
    //    mem.p_file->note(true," slices to potf exit ");
  }
}
