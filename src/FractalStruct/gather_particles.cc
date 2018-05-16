#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void gather_particles(Fractal_Memory& mem,Fractal& frac)
  {
    Group* p_fake_group=new Group;
    //    ofstream& FF=frac.p_file->FileFractal;
    if(!mem.MPIrun)
      return;
    int FractalNodes=mem.p_mess->FractalNodes;
    bool sendrad=mem.calc_shear;
    vector <double> pf(4);
    int totalI=0;
    int LOOPS=3;
    for(int LOOP=0;LOOP<LOOPS;LOOP++)
      {
	vector <int> counts_out(FractalNodes,0);
	vector <vector <int> > dataI_out(FractalNodes);
	vector <vector <double> > dataR_out(FractalNodes);
	vector <int> counts_in;
	vector <double> dataR_in;
	vector <int> dataI_in;
	int number=-1;
	int FR=-1;
	for(int particle=0; particle < frac.get_number_particles(); ++particle)
	  {
	    if(particle % LOOPS == LOOP)
	      {
		Particle* P=frac.particle_list[particle];
		P->get_world(FR,number);
		int lev=P->get_highest_level();
		if(P->get_p_highest_level_group()==0)
		  continue;
		if(number < 0)
		  continue;
		P->get_field_pf(pf);
		dataR_out[FR].push_back(pf[0]);
		dataR_out[FR].push_back(pf[1]);
		dataR_out[FR].push_back(pf[2]);
		dataR_out[FR].push_back(pf[3]);
		if(sendrad)
		  dataR_out[FR].push_back(P->get_rad_max());
		dataI_out[FR].push_back(number);
		dataI_out[FR].push_back(lev);
		counts_out[FR]++;
	      }
	  }
	int how_manyI=-1;
	int how_manyR=-1;
	int integers=2;
	int doubles=4;
	if(sendrad)
	  doubles=5;
	mem.p_mess->Send_Data_Some_How(1,counts_out,counts_in,integers,doubles,
				       dataI_out,dataI_in,how_manyI,
				       dataR_out,dataR_in,how_manyR);
	totalI+=how_manyI;
	dataR_out.clear();
	dataI_out.clear();
	int p2=0;
	int p4=0;
	Particle* P=0;
	for(int FR=0;FR<FractalNodes;FR++)
	  {
	    for(int c=0;c<counts_in[FR];c++)
	      {
		if(dataI_in[p2] >= 0)
		  {
		    P=frac.particle_list_world[dataI_in[p2]];
		    P->set_highest_level(dataI_in[p2+1]);
		    P->set_field_pf(dataR_in[p4],dataR_in[p4+1],dataR_in[p4+2],dataR_in[p4+3]);
		    if(sendrad)
		      P->set_rad_max(dataR_in[p4+4]);
		    P->set_p_highest_level_group(p_fake_group);
		  }
		p2+=integers;
		p4+=doubles;
	      }
	  }
      }
    delete p_fake_group;
    //    FF << " gather " << how_manyI << " " << frac.particle_list_world.size() << "\n";
    // frac.set_number_particles((totalI/2));
    frac.set_number_particles(frac.particle_list_world.size());
    frac.particle_list=frac.particle_list_world;
    // frac.particle_list_world.clear();
    clean_deque(frac.particle_list_world);
    size_t totaltmp=mem.p_mess->parts_tmp.size();
    size_t totaltmpp=mem.p_mess->parts_tmpp.size();
    for(int FR=0;FR<totaltmp;FR++)
      {
	Particle* pt=mem.p_mess->parts_tmp[FR];
	delete [] pt;
      }
    // mem.p_mess->parts_tmp.clear();
    clean_deque(mem.p_mess->parts_tmp);
    for(int FR=0;FR<totaltmpp;FR++)
      {
	Particle* ptp=mem.p_mess->parts_tmpp[FR];
	delete [] ptp;
      }
    // mem.p_mess->parts_tmpp.clear();
    clean_deque(mem.p_mess->parts_tmpp);
    remove_pseudo_particles(mem,frac);
    // frac.particle_list.resize((totalI/2));
    // frac.set_number_particles((totalI/2));
  }
}

