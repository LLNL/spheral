#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void gather_particles(Fractal_Memory& mem,Fractal& frac)
  {
    Group* p_fake_group=new Group;
    ofstream& FF=frac.p_file->FileFractal;
    //    ofstream& FP=frac.p_file->FileParticle;
    FF << " entered gather_particles " << "\n";
    if(!mem.MPIrun)
      return;
    //    int FractalRank=mem.p_mess->FractalRank;
    int FractalNodes=mem.p_mess->FractalNodes;
    bool sendrad=mem.calc_shear;
    vector <double> pf(4);
    vector <int>counts_out(FractalNodes);
    counts_out.assign(FractalNodes,0);
    vector <vector <int> > dataI_out(FractalNodes);
    vector <vector <double> > dataR_out(FractalNodes);
    int number=-1;
    int FR=-1;
    for(int particle=0; particle < frac.get_number_particles(); ++particle)
      {
	Particle* P=frac.particle_list[particle];
	P->get_world(FR,number);
	int lev=P->get_highest_level();
	//	FF << " gather00 " << P->get_p_highest_level_group() << " " << number << " " << FR << "\n";
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
    vector <int>counts_in(FractalNodes);
    counts_in.assign(FractalNodes,0);
    vector <double> dataR_in;
    vector <int> dataI_in;
    int how_manyI=-1;
    int how_manyR=-1;
    int integers=2;
    int doubles=4;
    if(sendrad)
      doubles=5;
    mem.p_file->note(true," gather particles a ");
    mem.p_mess->Send_Data_Some_How(counts_out,counts_in,integers,doubles,
				   dataI_out,dataI_in,how_manyI,
				   dataR_out,dataR_in,how_manyR);
    mem.p_file->note(true," gather particles c ");
    dataR_out.clear();
    dataI_out.clear();
    //    really_clear(dataR_out);
    //    really_clear(dataI_out);
    FF << " world size " << frac.particle_list_world.size() << "\n";
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
    delete p_fake_group;
    FF << " gather b " << how_manyI << " " << frac.particle_list_world.size() << "\n";
    frac.set_number_particles((how_manyI/2));
    frac.particle_list=frac.particle_list_world;
    frac.particle_list_world.clear();
    //    really_clear(frac.particle_list_world);
    FF << " gather e " << mem.p_mess->parts_tmp << "\n";
    Particle* pt=mem.p_mess->parts_tmp;
    delete [] pt;
    pt=0;
    Particle* ptp=mem.p_mess->parts_tmpp;
    delete [] ptp;
    ptp=0;
    remove_pseudo_particles(mem,frac);
    frac.particle_list.resize((how_manyI/2));
    frac.set_number_particles((how_manyI/2));
    FF << " finished remove_particles " << "\n";
  }
}

