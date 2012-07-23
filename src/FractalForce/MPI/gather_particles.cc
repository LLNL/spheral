#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void gather_particles(Fractal_Memory& mem,Fractal& frac)
  {
    if(!mem.MPIrun)
      return;
    int FractalRank=mem.p_mess->FractalRank;
    int FractalNodes=mem.p_mess->FractalNodes;
    vector <double> pos(3);
    vector <double> RealBox(6);
    RealBox=mem.RealBoxes[FractalRank];
    vector <double> pf(4);
    int* countsI_out=new int[FractalNodes];
    for(int FR=0;FR<FractalNodes;FR++)
      countsI_out[FR]=0;
    vector <vector <int> > potfI_out(FractalNodes);
    vector <vector <double> > potfR_out(FractalNodes);
    int number=-1;
    int FR=-1;
    int count_stay=0;
    for(int particle=0; particle < frac.get_number_particles(); ++particle)
      {
	Particle* P=frac.particle_list[particle];
	P->get_world(number,FR);
	P->get_pos(pos);
	if(pos[0] <  RealBox[0] ||
	   pos[0] >= RealBox[1] ||
	   pos[1] <  RealBox[2] ||
	   pos[1] >= RealBox[3] ||
	   pos[2] <  RealBox[4] ||
	   pos[2] >= RealBox[5]) 
	  continue;
	if(FR == FractalRank)
	  {
	    frac.particle_list_world[number]=P;
	    count_stay++;
	  }
	else
	  {
	    P->get_field_pf(pf);
	    potfR_out[FR].push_back(pf[0]);
	    potfR_out[FR].push_back(pf[1]);
	    potfR_out[FR].push_back(pf[2]);
	    potfR_out[FR].push_back(pf[3]);
	    potfI_out[FR].push_back(number);
	    countsI_out[FR]++;
	  }
      }
    int* countsI_in=new int[FractalNodes];
    double* potfR_in=0;
    int* potfI_in=0;
    int how_many=-1;
    mem.p_mess->How_Many_Particles_To_Send(countsI_out,countsI_in);
    mem.p_mess->Send_Particles_Somewhere(countsI_out,countsI_in,
					 potfR_out,potfI_out,
					 how_many,potfR_in,potfI_in);
    int totals=frac.get_number_particles_world();
    assert(count_stay+how_many==totals);
    int particle=0;
    int p4=-1;
    Particle* P=0;
    int pp=-1;
    for(int FR=0;FR<FractalNodes;FR++)
      {
	for(int c=0;c<countsI_in[FR];c++)
	  {
	    pp=potfI_in[particle];
	    assert(pp>=0);
	    assert(pp<totals);
	    P=frac.particle_list_world[pp];
	    p4=particle*4;
	    P->set_field_pf(potfR_in[p4],potfR_in[p4+1],potfR_in[p4+2],potfR_in[p4+3]);
	    particle++;
	  }
      }
    frac.set_number_particles(totals);
    frac.particle_list=frac.particle_list_world;
    frac.particle_list_world.clear();
    delete [] mem.p_mess->parts_tmp;
  }
  void left_right(Fractal& frac,vector <double>& pos_left,vector <double>& pos_right)
  {
    pos_left.assign(3,1.0e30);
    pos_right.assign(3,-1.0e30);
    vector <double> pos(3);
    for(int particle=0; particle < frac.get_number_particles(); ++particle)
      {
	(frac.particle_list[particle])->get_pos(pos);
	for(int ni=0;ni<3;ni++)
	  {
	    pos_left[ni]=min(pos_left[ni],pos[ni]);
	    pos_right[ni]=max(pos_right[ni],pos[ni]);
	  }
      }
  }	
  template <class T> bool overlap(vector <T>& xleft,vector <T>& xright,vector <T>& yleft,vector <T>& yright)
  {
    return (xleft[0] <= yright[0] && xright[0] >= yleft[0] && 
	    xleft[1] <= yright[1] && xright[1] >= yleft[1] &&
	    xleft[2] <= yright[2] && xright[2] >= yleft[2]);
  }
  template <class T> bool overlap(vector <T>& xleft,vector <T>& xright,vector <T>& yleftright)
  {
    return (xleft[0] <= yleftright[3] && xright[0] >= yleftright[0] && 
	    xleft[1] <= yleftright[4] && xright[1] >= yleftright[1] &&
	    xleft[2] <= yleftright[5] && xright[2] >= yleftright[2]);
  }
}

