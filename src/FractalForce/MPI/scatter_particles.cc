#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void scatter_particles(Fractal_Memory& mem,Fractal& frac)
  {
    if(frac.get_periodic())
      frac.wrap();
    if(!mem.MPIrun)
      return;
    int FractalRank=mem.p_mess->FractalRank;
    int FractalNodes=mem.p_mess->FractalNodes;
    vector <double> pos_left(3);
    vector <double> pos_right(3);
    vector <double> pos(3);
    vector <double> RealPBox(6);
    int* countsI_out=new int[FractalNodes];
    vector <vector <int> > posmI_out(FractalNodes);
    vector <vector <double> > posmR_out(FractalNodes);
    left_right(frac,pos_left,pos_right);
    for(int particle=0; particle < frac.get_number_particles(); ++particle)
      frac.particle_list[particle]->set_world(FractalRank,particle);
    frac.particle_list_world=frac.particle_list;
    frac.set_number_particles_world(frac.get_number_particles());
    frac.particle_list.clear();
    int count_stay=0;
    for(int FR=0;FR<FractalNodes;FR++)
      {
	RealPBox=mem.RealPBoxes[FR];
	int counts=0;
	if(overlap(pos_left,pos_right,RealPBox))
	  {
	    for(int particle=0; particle < frac.get_number_particles_world(); ++particle)
	      {
		Particle* P=frac.particle_list_world[particle];
		P->get_pos(pos);
		if(pos[0] <  RealPBox[0] ||
		   pos[0] >= RealPBox[1] ||
		   pos[1] <  RealPBox[2] ||
		   pos[1] >= RealPBox[3] ||
		   pos[2] <  RealPBox[4] ||
		   pos[2] >= RealPBox[5]) 
		  continue;
		if(FR == FractalRank)
		  {
		    frac.particle_list.push_back(P);
		    count_stay++;
		  }
		else
		  {
		    posmR_out[FR].push_back(pos[0]);
		    posmR_out[FR].push_back(pos[1]);
		    posmR_out[FR].push_back(pos[2]);
		    posmR_out[FR].push_back(P->get_mass());
		    posmI_out[FR].push_back(particle);
		    counts++;
		  }
	      }
	  }
	countsI_out[FR]=counts;
      }
    int* countsI_in=new int[FractalNodes];
    double* posmR_in=0;
    int* posmI_in=0;
    int how_many=-1;
    mem.p_mess->How_Many_Particles_To_Send(countsI_out,countsI_in);
    mem.p_mess->Send_Particles_Somewhere(countsI_out,countsI_in,
					 posmR_out,posmI_out,
					 how_many,posmR_in,posmI_in);
    frac.particle_list.resize(count_stay+how_many);
    Particle* particles_tmp=new Particle[how_many];
    mem.p_mess->parts_tmp=particles_tmp;
    int particle=0;
    int p4=0;
    Particle* P=0;
    for(int FR=0;FR<FractalNodes;FR++)
      {
	for(int c=0;c<countsI_in[FR];c++)
	  {
	    P=&particles_tmp[particle];
	    frac.particle_list[particle+count_stay]=P;
	    p4=particle*4;
	    P->set_posmIFR(posmR_in[p4],posmR_in[p4+1],posmR_in[p4+2],posmR_in[p4+3],posmI_in[particle],FR);
	    particle++;
	  }
      }
    frac.set_number_particles(particle+count_stay);
  }
}
