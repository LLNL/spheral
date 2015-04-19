#include "libs.hh"
#include "classes.hh"
#include "headers.hh"
namespace FractalSpace
{
  void scatter_particles(Fractal_Memory& mem,Fractal& frac)
  {
    assert(mem.FractalNodes==mem.p_mess->FractalNodes);
    ofstream& FF=frac.p_file->FileFractal;
    frac.p_file->FileFractal << " entered scatter particles " << endl;
    if(frac.get_periodic())
      {
	FF << " entered scatter particlelist size " << frac.particle_list.size() << endl;
	frac.wrap();
	FF << " enter add pseudo particles " << endl;
	add_pseudo_particles(mem,frac);
	FF << " exited scatter particlelist size " << frac.particle_list.size() << endl;
      }
    if(!mem.MPIrun)
      return;
    int FractalRank=mem.p_mess->FractalRank;
    int FractalNodes=mem.p_mess->FractalNodes;
    vector <double> pos_left(3);
    vector <double> pos_right(3);
    vector <double> pos(3);
    vector <double> RealBox(6);
    vector <double> RealPBox(6);
    vector <int> counts_out(FractalNodes);
    vector <vector <int> > dataI_out(FractalNodes);
    vector <vector <double> > dataR_out(FractalNodes);
    left_right(frac,pos_left,pos_right);
    FF << " poslr " << FractalRank;
    FF << " " << pos_left[0] << " " << pos_left[1] << " " << pos_left[2];
    FF << " " << pos_right[0] << " " << pos_right[1] << " " << pos_right[2] << endl;
    for(int particle=0; particle < frac.get_number_particles(); ++particle)
      frac.particle_list[particle]->set_world(FractalRank,particle);
    frac.particle_list_world=frac.particle_list;
    frac.set_number_particles_world(frac.get_number_particles());
    for(int FR=0;FR<FractalNodes;FR++)
      {
	counts_out[FR]=0;
	RealPBox=mem.RealPBoxes[FR];
	RealBox=mem.RealBoxes[FR];
	FF << " rpb " << FR << " " << FractalRank;
	FF << " " << RealPBox[0] << " " << RealPBox[1] << " " << RealPBox[2];
	FF << " " << RealPBox[3] << " " << RealPBox[4] << " " << RealPBox[5] << endl;
	if(!overlap(pos_left,pos_right,RealPBox))
	  continue;
	FF << " found an overlap " << FR << " " << FractalRank << endl;
	for(int particle=0; particle < frac.get_number_particles_world(); ++particle)
	  {
	    Particle* P=frac.particle_list_world[particle];
	    P->get_pos(pos);
	    if(pos[0] <  RealPBox[0] || pos[0] >= RealPBox[1] || 
	       pos[1] <  RealPBox[2] || pos[1] >= RealPBox[3] || 
	       pos[2] <  RealPBox[4] || pos[2] >= RealPBox[5]) 
	      continue;
	    int part=particle;
	    if(pos[0] <  RealBox[0] || pos[0] >= RealBox[1] || 
	       pos[1] <  RealBox[2] || pos[1] >= RealBox[3] || 
	       pos[2] <  RealBox[4] || pos[2] >= RealBox[5]) 
	      part=-particle-1;
	    dataR_out[FR].push_back(pos[0]);
	    dataR_out[FR].push_back(pos[1]);
	    dataR_out[FR].push_back(pos[2]);
	    dataR_out[FR].push_back(P->get_mass());
	    dataI_out[FR].push_back(part);
	    counts_out[FR]++;
	    //	    FF << " scata " << FR << " " << counts_out[FR] << " " << FractalRank << endl;
	  }
      }
    vector <int> counts_in(FractalNodes);
    vector <double> dataR_in;
    vector <int> dataI_in;
    int how_manyI=-1;
    int how_manyR=-1;
    int integers=1;
    int doubles=4;
    frac.timing(-1,33);
    mem.p_mess->Full_Stop();
    frac.timing(1,33);
    mem.p_file->note(true," scatter particles a ");
    mem.p_mess->How_Many_Things_To_Send(counts_out,counts_in);
    mem.p_file->note(true," scatter particles b ");
    mem.p_mess->Send_Data_Somewhere(counts_out,counts_in,integers,doubles,
				    dataI_out,dataI_in,how_manyI,
				    dataR_out,dataR_in,how_manyR);
    mem.p_file->note(true," scatter particles c ");
    dataR_out.clear();
    dataI_out.clear();
    frac.particle_list.resize(how_manyI);
    Particle* particles_tmp=new Particle[how_manyI];
    mem.p_mess->parts_tmp=particles_tmp;
    int field_length=4;
    if(mem.calc_density_particle)
      field_length=5;
    if(mem.calc_shear)
      field_length=6;
    const bool change_it= field_length > 4;
    int particle=0;
    int p4=-1;
    Particle* P=0;
    for(int FR=0;FR<FractalNodes;FR++)
      {
	for(int c=0;c<counts_in[FR];c++)
	  {
	    P=&particles_tmp[particle];
	    assert(P);
	    frac.particle_list[particle]=P;
	    p4=particle*4;
	    //	    FF << " datain " << particle << " " << p4 << " " << FR << " " << counts_in[FR] << endl;
	    //	    FF << "dataina " << P << " " << dataR_in[p4] << " " << dataR_in[p4+1]<< " " << dataR_in[p4+2] << " " << dataR_in[p4+3];
	    //	    FF << " " << dataI_in[particle] << endl;
	    P->set_posmIFR(dataR_in[p4],dataR_in[p4+1],dataR_in[p4+2],dataR_in[p4+3],dataI_in[particle],FR);
	    if(change_it)
	      P->field_resize(field_length);
	    particle++;
	  }
      }
    frac.set_number_particles(particle);
  }
}
